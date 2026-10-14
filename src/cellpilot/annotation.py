"""Strict, auditable CellTypist annotation.

The adapter deliberately supports one frozen model artifact.  It never asks
CellTypist to resolve a model name (which would use ``~/.celltypist``), never
enables majority voting for the primary prediction, and never invents a
fallback label when annotation fails.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse


MODEL_NAME = "Healthy_COVID19_PBMC.pkl"
MODEL_URL = (
    "https://celltypist.cog.sanger.ac.uk/models/COVID19_PBMC_Haniffa/v1/"
    "Healthy_COVID19_PBMC.pkl"
)
MODEL_SHA256 = "1eee68d23b4ac0f6999a0a92d6c950d671f31fcdb003ce0760c580b0bd6a3b79"
MODEL_VERSION = "v1"
MODEL_SOURCE = "https://doi.org/10.1038/s41591-021-01329-2"
MODEL_DATE = "2022-03-10 05:08:08.224597"
MODEL_N_CELL_TYPES = 51
MODEL_N_FEATURES = 3_443
_MAX_MODEL_BYTES = 128 * 1024 * 1024


class AnnotationInputError(ValueError):
    """The expression matrix is not a valid CellTypist query."""


class ModelArtifactError(RuntimeError):
    """The fixed model could not be acquired or failed integrity checks."""


class AnnotationExecutionError(RuntimeError):
    """CellTypist failed or returned a malformed result."""


@dataclass(frozen=True)
class ModelArtifact:
    """Verified local model plus provenance recorded beside it."""

    path: Path
    metadata_path: Path
    name: str
    url: str
    sha256: str
    size_bytes: int
    version: str
    source: str
    model_date: str
    declared_cell_types: int
    declared_features: int
    downloaded: bool

    def to_metadata(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["path"] = str(self.path)
        payload["metadata_path"] = str(self.metadata_path)
        return payload


@dataclass(frozen=True)
class AnnotationResult:
    """Per-cell primary labels/probabilities and the evidence needed to audit them."""

    predictions: pd.DataFrame
    probabilities: pd.DataFrame
    model_metadata: dict[str, Any]
    input_validation: dict[str, Any]


def _matrix_values(matrix: Any) -> np.ndarray:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    return np.asarray(values)


def _inverse_log1p_row_sums(matrix: Any, *, chunk_size: int = 2_048) -> np.ndarray:
    """Compute per-cell totals on the count scale without densifying the matrix."""

    n_rows = int(matrix.shape[0])
    totals = np.empty(n_rows, dtype=np.float64)
    for start in range(0, n_rows, chunk_size):
        stop = min(start + chunk_size, n_rows)
        block = matrix[start:stop]
        if sparse.issparse(block):
            block = block.tocsr(copy=True)
            block.data = np.expm1(block.data.astype(np.float64, copy=False))
            totals[start:stop] = np.asarray(block.sum(axis=1)).ravel()
        else:
            totals[start:stop] = np.expm1(
                np.asarray(block, dtype=np.float64)
            ).sum(axis=1)
    return totals


def validate_celltypist_input(
    adata: Any,
    *,
    target_sum: float = 10_000.0,
    min_genes: int = 1_000,
    normalization_rtol: float = 2e-3,
    normalization_atol: float = 2.0,
) -> dict[str, Any]:
    """Validate non-negative, finite, full-gene log1p-normalized input.

    "Full-gene" cannot be proven from a matrix alone.  This function therefore
    enforces three observable safeguards: a minimum feature count, rejection
    when ``.raw`` exposes a larger feature universe than ``.X``, and a per-cell
    inverse-log total close to 10,000.  :func:`annotate_cells` adds a fourth
    safeguard by checking coverage of the frozen model's 3,443 features.
    """

    if target_sum <= 0:
        raise ValueError("target_sum must be positive")
    if min_genes < 1:
        raise ValueError("min_genes must be positive")
    for attribute in ("X", "n_obs", "n_vars", "obs_names", "var_names"):
        if not hasattr(adata, attribute):
            raise AnnotationInputError(f"input is missing AnnData attribute {attribute!r}")
    if int(adata.n_obs) == 0 or int(adata.n_vars) == 0:
        raise AnnotationInputError("CellTypist input is empty")
    if int(adata.n_vars) < min_genes:
        raise AnnotationInputError(
            f"CellTypist requires full-gene input; found {adata.n_vars} genes, "
            f"below the configured minimum of {min_genes}"
        )
    if not getattr(adata.obs_names, "is_unique", False):
        raise AnnotationInputError("cell identifiers must be unique")
    if not getattr(adata.var_names, "is_unique", False):
        raise AnnotationInputError("gene identifiers must be unique")
    if any(not str(name).strip() for name in adata.var_names):
        raise AnnotationInputError("gene identifiers must be non-empty")

    raw = getattr(adata, "raw", None)
    if raw is not None and int(raw.n_vars) > int(adata.n_vars):
        raise AnnotationInputError(
            "CellTypist input appears feature-subsetted: .raw contains "
            f"{raw.n_vars} genes but .X contains {adata.n_vars}; pass the full-gene "
            "log1p-normalized matrix explicitly"
        )

    values = _matrix_values(adata.X)
    if not np.issubdtype(values.dtype, np.number):
        raise AnnotationInputError("expression values must be numeric")
    if not np.isfinite(values).all():
        raise AnnotationInputError("expression matrix contains non-finite values")
    minimum = float(values.min()) if values.size else 0.0
    maximum = float(values.max()) if values.size else 0.0
    if minimum < 0:
        raise AnnotationInputError(
            "expression matrix contains negative values; scaled/regressed data are invalid"
        )
    # No normalized count for one gene can exceed the complete target library.
    theoretical_max = float(np.log1p(target_sum))
    if maximum > theoretical_max + 0.02:
        raise AnnotationInputError(
            "expression matrix is not plausible log1p-normalized data: observed "
            f"maximum {maximum:.4g} exceeds log1p(target_sum) {theoretical_max:.4g}"
        )

    totals = _inverse_log1p_row_sums(adata.X)
    if not np.isfinite(totals).all():
        raise AnnotationInputError("inverse-log expression totals are non-finite")
    if (totals <= 0).any():
        raise AnnotationInputError("expression matrix contains an empty cell")
    valid_totals = np.isclose(
        totals,
        target_sum,
        rtol=normalization_rtol,
        atol=normalization_atol,
    )
    if not valid_totals.all():
        bad = totals[~valid_totals]
        raise AnnotationInputError(
            "expression matrix is not log1p-normalized to the required per-cell "
            f"target {target_sum:g}; {bad.size}/{totals.size} cells fail "
            f"(observed range {totals.min():.3f}-{totals.max():.3f})"
        )

    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "matrix_format": "sparse" if sparse.issparse(adata.X) else "dense",
        "minimum": minimum,
        "maximum": maximum,
        "target_sum": float(target_sum),
        "inverse_log_total_min": float(totals.min()),
        "inverse_log_total_max": float(totals.max()),
        "full_gene_checks": [
            f"n_genes >= {min_genes}",
            ".raw does not contain more genes than .X",
            "all inverse-log per-cell totals match target_sum",
            "model-feature coverage checked before prediction",
        ],
    }


def _sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_fixed_model(url: str, timeout_seconds: int) -> Any:
    """Open the fixed HTTPS artifact using Requests' certifi trust bundle."""

    import requests

    response = requests.get(
        url,
        stream=True,
        timeout=timeout_seconds,
        headers={"User-Agent": "CellPilot-AI/0.1"},
    )
    response.raise_for_status()
    return response


def _stream_model(response: Any, destination: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with destination.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if not chunk:
                continue
            size += len(chunk)
            if size > _MAX_MODEL_BYTES:
                raise ModelArtifactError(
                    f"model response exceeds safety limit of {_MAX_MODEL_BYTES} bytes"
                )
            digest.update(chunk)
            handle.write(chunk)
    return digest.hexdigest(), size


def ensure_model(
    model_dir: str | Path,
    *,
    model_url: str | None = None,
    expected_sha256: str | None = None,
    timeout_seconds: int = 120,
    force_download: bool = False,
) -> ModelArtifact:
    """Acquire and verify the one fixed model inside ``model_dir`` only.

    ``model_url`` and ``expected_sha256`` exist for configuration compatibility,
    but any non-matching value is rejected.  The known digest makes a change at
    the remote URL fail closed instead of silently changing benchmark results.
    """

    requested_url = MODEL_URL if model_url is None else model_url
    requested_sha = MODEL_SHA256 if expected_sha256 is None else expected_sha256.lower()
    if requested_url != MODEL_URL:
        raise ModelArtifactError(
            f"only the frozen model URL is permitted; received {requested_url!r}"
        )
    if requested_sha != MODEL_SHA256:
        raise ModelArtifactError(
            "configured model SHA-256 does not match the frozen Healthy_COVID19_PBMC "
            f"artifact ({MODEL_SHA256})"
        )
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")

    directory = Path(model_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    model_path = directory / MODEL_NAME
    metadata_path = directory / f"{MODEL_NAME}.metadata.json"
    if model_path.parent.resolve() != directory:
        raise ModelArtifactError("model destination escaped model_dir")

    downloaded = False
    if model_path.exists() and not force_download:
        if not model_path.is_file() or model_path.is_symlink():
            raise ModelArtifactError(f"cached model is not a regular file: {model_path}")
        observed_sha = _sha256_file(model_path)
        size_bytes = model_path.stat().st_size
        if observed_sha != MODEL_SHA256:
            raise ModelArtifactError(
                f"cached model checksum mismatch: expected {MODEL_SHA256}, "
                f"observed {observed_sha}; remove it or use force_download=True"
            )
    else:
        temporary_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=directory, prefix=f".{MODEL_NAME}.", suffix=".part", delete=False
            ) as temporary:
                temporary_path = Path(temporary.name)
            try:
                response_context = _open_fixed_model(MODEL_URL, timeout_seconds)
                with response_context as response:
                    observed_sha, size_bytes = _stream_model(response, temporary_path)
            except Exception as exc:
                if isinstance(exc, ModelArtifactError):
                    raise
                raise ModelArtifactError(f"failed to download fixed CellTypist model: {exc}") from exc
            if observed_sha != MODEL_SHA256:
                raise ModelArtifactError(
                    f"downloaded model checksum mismatch: expected {MODEL_SHA256}, "
                    f"observed {observed_sha}"
                )
            temporary_path.replace(model_path)
            temporary_path = None
            downloaded = True
        finally:
            if temporary_path is not None and temporary_path.exists():
                temporary_path.unlink()

    artifact = ModelArtifact(
        path=model_path,
        metadata_path=metadata_path,
        name=MODEL_NAME,
        url=MODEL_URL,
        sha256=observed_sha,
        size_bytes=int(size_bytes),
        version=MODEL_VERSION,
        source=MODEL_SOURCE,
        model_date=MODEL_DATE,
        declared_cell_types=MODEL_N_CELL_TYPES,
        declared_features=MODEL_N_FEATURES,
        downloaded=downloaded,
    )
    persisted = artifact.to_metadata()
    persisted["checksum_algorithm"] = "sha256"
    persisted["verified_at_utc"] = datetime.now(timezone.utc).isoformat()
    persisted["policy"] = "fixed URL and fixed digest; no CellTypist model-name resolution"
    temporary_metadata = metadata_path.with_suffix(metadata_path.suffix + ".part")
    temporary_metadata.write_text(
        json.dumps(persisted, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary_metadata.replace(metadata_path)
    return artifact


def _load_celltypist(model_dir: Path) -> Any:
    # CellTypist creates its configured folder at import time.  Keep that write
    # below the caller-provided workspace model directory rather than $HOME.
    runtime_dir = model_dir / ".celltypist_runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CELLTYPIST_FOLDER"] = str(runtime_dir)
    try:
        import celltypist
    except Exception as exc:  # pragma: no cover - environment-specific import failures
        raise AnnotationExecutionError(f"CellTypist import failed: {exc}") from exc
    return celltypist


def _model_features(model: Any) -> np.ndarray:
    if hasattr(model, "features"):
        return np.asarray(model.features, dtype=str)
    classifier = getattr(model, "classifier", None)
    if classifier is not None and hasattr(classifier, "features"):
        return np.asarray(classifier.features, dtype=str)
    raise AnnotationExecutionError("loaded CellTypist model does not expose feature names")


def _model_cell_types(model: Any) -> np.ndarray:
    if hasattr(model, "cell_types"):
        return np.asarray(model.cell_types, dtype=str)
    classifier = getattr(model, "classifier", None)
    if classifier is not None and hasattr(classifier, "classes_"):
        return np.asarray(classifier.classes_, dtype=str)
    raise AnnotationExecutionError("loaded CellTypist model does not expose cell types")


def annotate_cells(
    adata: Any,
    model_dir: str | Path,
    *,
    model_url: str | None = None,
    expected_sha256: str | None = None,
    target_sum: float = 10_000.0,
    min_genes: int = 1_000,
    min_model_feature_fraction: float = 0.5,
) -> AnnotationResult:
    """Run fixed-model, per-cell CellTypist annotation without majority voting.

    Exceptions are propagated as :class:`AnnotationExecutionError`; callers
    receive no fabricated labels and no alternate annotator is invoked.
    """

    if not 0 < min_model_feature_fraction <= 1:
        raise ValueError("min_model_feature_fraction must be in (0, 1]")
    validation = validate_celltypist_input(
        adata, target_sum=target_sum, min_genes=min_genes
    )
    artifact = ensure_model(
        model_dir,
        model_url=model_url,
        expected_sha256=expected_sha256,
    )
    celltypist = _load_celltypist(artifact.path.parent)
    try:
        model = celltypist.models.Model.load(model=str(artifact.path))
    except Exception as exc:
        raise AnnotationExecutionError(
            f"verified model could not be loaded by CellTypist: {exc}"
        ) from exc

    features = _model_features(model)
    model_cell_types = _model_cell_types(model)
    query_genes = np.asarray(adata.var_names, dtype=str)
    overlap = int(np.isin(features, query_genes).sum())
    overlap_fraction = overlap / len(features) if len(features) else 0.0
    validation["model_features"] = int(len(features))
    validation["model_feature_overlap"] = overlap
    validation["model_feature_overlap_fraction"] = float(overlap_fraction)
    if overlap_fraction < min_model_feature_fraction:
        raise AnnotationInputError(
            "query does not provide sufficient frozen-model feature coverage for "
            f"full-gene annotation: {overlap}/{len(features)} "
            f"({overlap_fraction:.1%}) < {min_model_feature_fraction:.1%}"
        )

    try:
        raw_result = celltypist.annotate(
            adata.copy(),
            model=model,
            mode="best match",
            majority_voting=False,
        )
    except Exception as exc:
        raise AnnotationExecutionError(
            f"CellTypist annotation failed; no labels were substituted: {exc}"
        ) from exc

    labels = getattr(raw_result, "predicted_labels", None)
    probabilities = getattr(raw_result, "probability_matrix", None)
    if not isinstance(labels, pd.DataFrame) or "predicted_labels" not in labels.columns:
        raise AnnotationExecutionError(
            "CellTypist result lacks the required predicted_labels column"
        )
    if not isinstance(probabilities, pd.DataFrame) or probabilities.empty:
        raise AnnotationExecutionError("CellTypist result lacks a probability matrix")
    expected_index = pd.Index(adata.obs_names.astype(str))
    labels = labels.copy()
    probabilities = probabilities.copy()
    labels.index = labels.index.astype(str)
    probabilities.index = probabilities.index.astype(str)
    if not labels.index.equals(expected_index) or not probabilities.index.equals(expected_index):
        raise AnnotationExecutionError(
            "CellTypist output cell order/identifiers do not exactly match the input"
        )
    probability_values = probabilities.to_numpy(dtype=float)
    if not np.isfinite(probability_values).all():
        raise AnnotationExecutionError("CellTypist returned non-finite probabilities")
    if ((probability_values < 0) | (probability_values > 1)).any():
        raise AnnotationExecutionError("CellTypist returned probabilities outside [0, 1]")

    predicted = labels["predicted_labels"].astype("string")
    if predicted.isna().any() or (predicted.str.len() == 0).any():
        raise AnnotationExecutionError("CellTypist returned missing/empty primary labels")
    unknown_outputs = sorted(set(predicted.astype(str)) - set(model_cell_types))
    if unknown_outputs:
        raise AnnotationExecutionError(
            f"CellTypist returned labels absent from the loaded model: {unknown_outputs}"
        )
    per_cell = pd.DataFrame(
        {
            "predicted_label": predicted.to_numpy(dtype=str),
            "confidence": probabilities.max(axis=1).to_numpy(dtype=float),
        },
        index=expected_index,
    )
    per_cell.index.name = "cell_id"

    metadata = artifact.to_metadata()
    metadata.update(
        {
            "celltypist_version": str(getattr(celltypist, "__version__", "unknown")),
            "observed_model_features": int(len(features)),
            "observed_model_cell_types": [str(value) for value in model_cell_types],
            "prediction_mode": "best match",
            "majority_voting": False,
            "fallback_labels_used": False,
        }
    )
    return AnnotationResult(
        predictions=per_cell,
        probabilities=probabilities,
        model_metadata=metadata,
        input_validation=validation,
    )


__all__ = [
    "AnnotationExecutionError",
    "AnnotationInputError",
    "AnnotationResult",
    "MODEL_NAME",
    "MODEL_SHA256",
    "MODEL_URL",
    "ModelArtifact",
    "ModelArtifactError",
    "annotate_cells",
    "ensure_model",
    "validate_celltypist_input",
]
