"""Independent five-class FACS-enriched PBMC annotation benchmark.

Reference labels are stored in :class:`BenchmarkCohort`, never in the AnnData
object handed to an annotator.  This prevents accidental use of ``.obs`` labels
as model inputs or post-hoc label substitution.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import tarfile
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Mapping

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from sklearn.metrics import (
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_fscore_support,
    recall_score,
)

from .data import validate_raw_counts, verify_extracted_tar


BROAD_CLASSES = (
    "B cell",
    "Monocyte",
    "CD4+ T cell",
    "CD8+ T cell",
    "NK cell",
)
UNMAPPED_LABEL = "__unmapped__"
DEFAULT_LABEL_MAP_PATH = Path(__file__).with_name("resources") / "benchmark_label_map.yaml"
_DOWNLOAD_LIMIT_BYTES = 4 * 1024**3
_EXTRACT_LIMIT_BYTES = 16 * 1024**3
_EXTRACT_MEMBER_LIMIT = 512


@dataclass(frozen=True)
class FacsDatasetSpec:
    """Frozen source and biological caveats for one official 10x panel."""

    dataset_id: str
    reference_label: str
    landing_page: str
    download_url: str
    reported_cells: int
    reported_purity_percent: float
    enrichment: str
    archive_sha256: str
    license: str = "CC BY 4.0"


FACS_DATASETS: dict[str, FacsDatasetSpec] = {
    "b_cells": FacsDatasetSpec(
        dataset_id="b_cells",
        reference_label="B cell",
        landing_page="https://www.10xgenomics.com/datasets/cd-19-plus-b-cells-1-standard-1-1-0",
        download_url=(
            "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/"
            "b_cells_filtered_gene_bc_matrices.tar.gz"
        ),
        reported_cells=10_085,
        reported_purity_percent=100.0,
        enrichment="FACS-enriched CD19+ B cells",
        archive_sha256="fdc8a77c9bfc0e78c378850df5caa77cf54352c3e6f82eea8acc35032089a622",
    ),
    "cd14_monocytes": FacsDatasetSpec(
        dataset_id="cd14_monocytes",
        reference_label="Monocyte",
        landing_page="https://www.10xgenomics.com/datasets/cd-14-plus-monocytes-1-standard-1-1-0",
        download_url=(
            "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/"
            "cd14_monocytes_filtered_gene_bc_matrices.tar.gz"
        ),
        reported_cells=2_612,
        reported_purity_percent=98.0,
        enrichment="FACS-enriched CD14+ monocytes",
        archive_sha256="dcdc45e252a838b46c4d3672647218602f817904fc933ae12178ba9a3fa48b99",
    ),
    "cd4_t_helper": FacsDatasetSpec(
        dataset_id="cd4_t_helper",
        reference_label="CD4+ T cell",
        landing_page="https://www.10xgenomics.com/datasets/cd-4-plus-helper-t-cells-1-standard-1-1-0",
        download_url=(
            "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd4_t_helper/"
            "cd4_t_helper_filtered_gene_bc_matrices.tar.gz"
        ),
        reported_cells=11_213,
        reported_purity_percent=99.0,
        enrichment="FACS-enriched CD4+ helper T cells",
        archive_sha256="3dd88aa2e33f5f833eef265fd32f564cdb81703c2e391703b0648f808d06ee93",
    ),
    "cd8_cytotoxic_t": FacsDatasetSpec(
        dataset_id="cd8_cytotoxic_t",
        reference_label="CD8+ T cell",
        landing_page="https://www.10xgenomics.com/datasets/cd-8-plus-cytotoxic-t-cells-1-standard-1-1-0",
        download_url=(
            "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cytotoxic_t/"
            "cytotoxic_t_filtered_gene_bc_matrices.tar.gz"
        ),
        reported_cells=10_209,
        reported_purity_percent=98.0,
        enrichment="FACS-enriched CD8+ cytotoxic T cells",
        archive_sha256="e6557491c09c1abdd4bf0efd0afa6139c14e39c11402d022b1601b4a07132d02",
    ),
    "cd56_nk": FacsDatasetSpec(
        dataset_id="cd56_nk",
        reference_label="NK cell",
        landing_page="https://www.10xgenomics.com/datasets/cd-56-plus-natural-killer-cells-1-standard-1-1-0",
        download_url=(
            "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd56_nk/"
            "cd56_nk_filtered_gene_bc_matrices.tar.gz"
        ),
        reported_cells=8_385,
        reported_purity_percent=92.0,
        enrichment="FACS-enriched CD56+ natural killer cells",
        archive_sha256="17765693a1a67ce4e3b069c415d348702029a13a9dbc19885433397140cb5e00",
    ),
}


BENCHMARK_LIMITATIONS = (
    "FACS enrichment labels are experimental population labels, not perfect "
    "single-cell ground truth; reported purity ranges from 92% to 100%.",
    "Each broad class is a separate 10x library (and may be a separate donor), "
    "so biological class is confounded with library/donor and technical effects.",
    "Inputs are Cell Ranger filtered gene-barcode matrices after cell calling, "
    "not FASTQ files or unfiltered droplets.",
    "The benchmark covers only five broad blood lineages and does not establish "
    "fine-subtype performance or out-of-reference rejection performance.",
    "The frozen CellTypist model cites Stephenson et al. 2021 rather than these "
    "10x panels; this reduces known study-level leakage but does not prove the "
    "absence of every indirect training-data relationship.",
    "Cell-level metrics do not estimate donor-to-donor generalization because "
    "there are no independent biological replicates within each class panel.",
)


class BenchmarkDataError(RuntimeError):
    """Dataset acquisition, integrity, or cohort construction failed."""


class UnsafeArchiveError(BenchmarkDataError):
    """A tar member attempted traversal, linking, or unsafe extraction."""


@dataclass(frozen=True)
class FetchedFacsDataset:
    spec: FacsDatasetSpec
    adata: Any
    archive_path: Path
    matrix_dir: Path
    manifest_path: Path
    sha256: str


@dataclass(frozen=True)
class BenchmarkCohort:
    """Sanitized query matrix and separately held reference information."""

    adata: Any
    reference_labels: pd.Series
    sample_metadata: pd.DataFrame
    metadata: dict[str, Any]


@dataclass(frozen=True)
class LabelMap:
    schema_version: str
    model_name: str
    model_sha256: str
    target_classes: tuple[str, ...]
    mappings: dict[str, str]
    explicitly_unmapped: dict[str, str]
    path: Path


@dataclass(frozen=True)
class BenchmarkResult:
    metrics: dict[str, Any]
    per_class: pd.DataFrame
    confusion_matrix: pd.DataFrame
    confusion_matrix_normalized: pd.DataFrame
    per_cell: pd.DataFrame
    metadata: dict[str, Any]


def _sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_fixed_url(url: str, destination: Path, timeout_seconds: int) -> str:
    """Atomically stream one allow-listed source with a bounded response size."""

    import requests

    if url not in {spec.download_url for spec in FACS_DATASETS.values()}:
        raise BenchmarkDataError(f"URL is not an allow-listed FACS dataset source: {url}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    digest = hashlib.sha256()
    total = 0
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".part",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
        with requests.get(
            url,
            stream=True,
            timeout=timeout_seconds,
            headers={"User-Agent": "CellPilot-AI/0.1"},
        ) as response:
            response.raise_for_status()
            with temporary.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > _DOWNLOAD_LIMIT_BYTES:
                        raise BenchmarkDataError(
                            f"download exceeds {_DOWNLOAD_LIMIT_BYTES} byte safety limit"
                        )
                    digest.update(chunk)
                    output.write(chunk)
        temporary.replace(destination)
        temporary = None
        return digest.hexdigest()
    except Exception as exc:
        if isinstance(exc, BenchmarkDataError):
            raise
        raise BenchmarkDataError(f"failed to download {url}: {exc}") from exc
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _validate_member_name(name: str) -> PurePosixPath:
    normalized = name.replace("\\", "/")
    member_path = PurePosixPath(normalized)
    if (
        normalized in {"", "."}
        or member_path.is_absolute()
        or ".." in member_path.parts
        or not member_path.parts
        or member_path.parts[0].endswith(":")
    ):
        raise UnsafeArchiveError(f"unsafe path in tar archive: {name!r}")
    return member_path


def safe_extract_tar(archive: str | Path, destination: str | Path) -> None:
    """Extract only directories/regular files, rejecting traversal and links."""

    archive_path = Path(archive)
    output_root = Path(destination).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode="r:*") as bundle:
        members = bundle.getmembers()
        if len(members) > _EXTRACT_MEMBER_LIMIT:
            raise UnsafeArchiveError(
                f"archive contains {len(members)} members; limit is {_EXTRACT_MEMBER_LIMIT}"
            )
        declared_size = 0
        for member in members:
            member_path = _validate_member_name(member.name)
            if not (member.isdir() or member.isfile()):
                raise UnsafeArchiveError(
                    f"archive member is not a regular file/directory: {member.name!r}"
                )
            if member.size < 0:
                raise UnsafeArchiveError(f"archive member has negative size: {member.name!r}")
            declared_size += member.size
            if declared_size > _EXTRACT_LIMIT_BYTES:
                raise UnsafeArchiveError(
                    f"declared extracted size exceeds {_EXTRACT_LIMIT_BYTES} bytes"
                )
            target = (output_root / Path(*member_path.parts)).resolve()
            if target != output_root and output_root not in target.parents:
                raise UnsafeArchiveError(f"archive member escapes destination: {member.name!r}")

        actual_size = 0
        for member in members:
            member_path = _validate_member_name(member.name)
            target = output_root / Path(*member_path.parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            source = bundle.extractfile(member)
            if source is None:
                raise UnsafeArchiveError(f"could not read archive member: {member.name!r}")
            with source, target.open("xb") as output:
                while True:
                    chunk = source.read(1024 * 1024)
                    if not chunk:
                        break
                    actual_size += len(chunk)
                    if actual_size > _EXTRACT_LIMIT_BYTES:
                        raise UnsafeArchiveError("actual extracted size exceeds safety limit")
                    output.write(chunk)
            if target.stat().st_size != member.size:
                raise UnsafeArchiveError(
                    f"extracted size mismatch for archive member {member.name!r}"
                )


def _find_10x_matrix_dir(root: Path) -> Path:
    candidates = sorted(
        {
            path.parent
            for pattern in ("matrix.mtx", "matrix.mtx.gz")
            for path in root.rglob(pattern)
        }
    )
    valid = []
    for candidate in candidates:
        genes = any((candidate / name).is_file() for name in ("genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz"))
        barcodes = any((candidate / name).is_file() for name in ("barcodes.tsv", "barcodes.tsv.gz"))
        if genes and barcodes:
            valid.append(candidate)
    if len(valid) != 1:
        raise BenchmarkDataError(
            f"expected exactly one complete 10x matrix directory below {root}, found {len(valid)}"
        )
    return valid[0]


def fetch_facs_dataset(
    dataset_id: str,
    data_dir: str | Path,
    *,
    timeout_seconds: int = 300,
    expected_sha256: str | None = None,
) -> FetchedFacsDataset:
    """Download, safely extract, load, and validate one official 10x panel."""

    if dataset_id not in FACS_DATASETS:
        raise KeyError(
            f"unknown FACS dataset {dataset_id!r}; choose from {sorted(FACS_DATASETS)}"
        )
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    pinned_sha256 = expected_sha256 or FACS_DATASETS[dataset_id].archive_sha256
    pinned_sha256 = pinned_sha256.strip().lower()
    if len(pinned_sha256) != 64 or any(
        value not in "0123456789abcdef" for value in pinned_sha256
    ):
        raise ValueError("expected_sha256 must be 64 lowercase hexadecimal characters")

    import scanpy as sc

    spec = FACS_DATASETS[dataset_id]
    root = Path(data_dir).expanduser().resolve()
    archive_dir = root / "archives"
    archive = archive_dir / f"{dataset_id}.tar.gz"
    extracted = root / "extracted" / dataset_id
    manifest_path = root / "manifests" / f"{dataset_id}.json"
    root.mkdir(parents=True, exist_ok=True)

    if archive.exists():
        if not archive.is_file() or archive.is_symlink():
            raise BenchmarkDataError(f"cached archive is not a regular file: {archive}")
        observed_sha = _sha256_file(archive)
    else:
        observed_sha = _download_fixed_url(spec.download_url, archive, timeout_seconds)
    if observed_sha != pinned_sha256:
        raise BenchmarkDataError(
            f"archive checksum mismatch for {dataset_id}: expected {pinned_sha256}, "
            f"observed {observed_sha}"
        )

    if extracted.exists():
        matrix_dir = _find_10x_matrix_dir(extracted)
    else:
        extracted.parent.mkdir(parents=True, exist_ok=True)
        temporary_root = Path(
            tempfile.mkdtemp(prefix=f".{dataset_id}.", dir=extracted.parent)
        )
        try:
            safe_extract_tar(archive, temporary_root)
            _find_10x_matrix_dir(temporary_root)
            temporary_root.replace(extracted)
        except Exception:
            shutil.rmtree(temporary_root, ignore_errors=True)
            raise
        matrix_dir = _find_10x_matrix_dir(extracted)

    extracted_digests = verify_extracted_tar(archive, extracted)
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    validation = validate_raw_counts(adata)
    manifest = {
        "dataset_id": dataset_id,
        "reference_label": spec.reference_label,
        "landing_page": spec.landing_page,
        "download_url": spec.download_url,
        "archive_sha256": observed_sha,
        "archive_sha256_policy": (
            "Checksum observed during the first verified project retrieval on "
            "2026-10-14 and subsequently frozen in source; not represented as an "
            "upstream-published checksum."
        ),
        "archive_size_bytes": archive.stat().st_size,
        "license": spec.license,
        "reported_cells": spec.reported_cells,
        "reported_purity_percent": spec.reported_purity_percent,
        "enrichment": spec.enrichment,
        "input_kind": "Cell Ranger filtered gene-barcode count matrix",
        "extracted_member_sha256": extracted_digests,
        "validation": validation,
        "limitations": list(BENCHMARK_LIMITATIONS),
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_manifest = manifest_path.with_suffix(".json.part")
    temporary_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary_manifest.replace(manifest_path)
    return FetchedFacsDataset(
        spec=spec,
        adata=adata,
        archive_path=archive,
        matrix_dir=matrix_dir,
        manifest_path=manifest_path,
        sha256=observed_sha,
    )


def fetch_all_facs_datasets(
    data_dir: str | Path,
    *,
    timeout_seconds: int = 300,
    expected_sha256: Mapping[str, str] | None = None,
) -> dict[str, FetchedFacsDataset]:
    """Fetch all five panels in a fixed order."""

    checksums = dict(expected_sha256 or {})
    unknown = set(checksums) - set(FACS_DATASETS)
    if unknown:
        raise KeyError(f"checksums supplied for unknown datasets: {sorted(unknown)}")
    return {
        dataset_id: fetch_facs_dataset(
            dataset_id,
            data_dir,
            timeout_seconds=timeout_seconds,
            expected_sha256=checksums.get(dataset_id),
        )
        for dataset_id in FACS_DATASETS
    }


def _dataset_seed(random_seed: int, dataset_id: str) -> int:
    digest = hashlib.sha256(f"{random_seed}:{dataset_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False)


def build_facs_cohort(
    datasets: Mapping[str, Any],
    *,
    max_cells_per_class: int | None = None,
    random_seed: int = 42,
) -> BenchmarkCohort:
    """Subsample each class independently, then concatenate sanitized matrices."""

    import anndata as ad

    if max_cells_per_class is not None and max_cells_per_class <= 0:
        raise ValueError("max_cells_per_class must be positive or None")
    if random_seed < 0:
        raise ValueError("random_seed must be non-negative")
    missing = set(FACS_DATASETS) - set(datasets)
    extra = set(datasets) - set(FACS_DATASETS)
    if missing or extra:
        raise BenchmarkDataError(
            f"datasets must contain exactly the five frozen panels; missing={sorted(missing)}, "
            f"extra={sorted(extra)}"
        )

    sanitized: list[Any] = []
    metadata_records: list[dict[str, Any]] = []
    selected_counts: dict[str, int] = {}
    for dataset_id, spec in FACS_DATASETS.items():
        source = datasets[dataset_id]
        validate_raw_counts(source)
        if not source.obs_names.is_unique or not source.var_names.is_unique:
            raise BenchmarkDataError(
                f"{dataset_id} must have unique cell and gene identifiers"
            )
        n_available = int(source.n_obs)
        n_selected = (
            n_available
            if max_cells_per_class is None
            else min(n_available, max_cells_per_class)
        )
        if n_selected == 0:
            raise BenchmarkDataError(f"{dataset_id} has no cells")
        if n_selected < n_available:
            rng = np.random.default_rng(_dataset_seed(random_seed, dataset_id))
            positions = np.sort(
                rng.choice(n_available, size=n_selected, replace=False)
            )
        else:
            positions = np.arange(n_available)
        subset = source[positions, :]
        original_cell_ids = source.obs_names[positions].astype(str)
        # Reconstruct AnnData from X + the gene index only.  No source .obs,
        # .uns, embeddings, or layers can carry a reference label downstream.
        clean = ad.AnnData(
            X=subset.X.copy(),
            obs=pd.DataFrame(
                index=[f"temporary_{len(sanitized)}_{i}" for i in range(n_selected)]
            ),
            var=pd.DataFrame(index=subset.var_names.astype(str)),
        )
        sanitized.append(clean)
        selected_counts[spec.reference_label] = n_selected
        metadata_records.extend(
            {
                "dataset_id": dataset_id,
                "original_cell_id": original_cell_id,
                "reference_label": spec.reference_label,
                "reported_purity_percent": spec.reported_purity_percent,
            }
            for original_cell_id in original_cell_ids
        )

    combined = ad.concat(
        sanitized,
        axis=0,
        join="inner",
        merge="same",
        index_unique=None,
    )
    if combined.n_vars == 0:
        raise BenchmarkDataError("the five panels have no common genes")
    anonymous_ids = pd.Index(
        [f"benchmark_cell_{index:07d}" for index in range(combined.n_obs)],
        name="cell_id",
    )
    combined.obs = pd.DataFrame(index=anonymous_ids)
    combined.uns.clear()
    combined.obsm.clear()
    combined.obsp.clear()

    sample_metadata = pd.DataFrame(metadata_records, index=anonymous_ids)
    reference_labels = sample_metadata["reference_label"].copy()
    reference_labels.name = "reference_label"
    metadata = {
        "benchmark": "10x FACS-enriched five-class PBMC",
        "target_classes": list(BROAD_CLASSES),
        "n_cells": int(combined.n_obs),
        "n_common_genes": int(combined.n_vars),
        "max_cells_per_class": max_cells_per_class,
        "random_seed": int(random_seed),
        "sampling": "independent deterministic subsampling per dataset before concatenation",
        "selected_cells_per_class": selected_counts,
        "reference_label_storage": "external pandas Series; absent from annotator AnnData",
        "limitations": list(BENCHMARK_LIMITATIONS),
    }
    return BenchmarkCohort(
        adata=combined,
        reference_labels=reference_labels,
        sample_metadata=sample_metadata,
        metadata=metadata,
    )


def normalize_log1p_for_celltypist(
    adata: Any, *, target_sum: float = 10_000.0
) -> Any:
    """Return a sanitized full-gene query normalized exactly once."""

    if target_sum <= 0:
        raise ValueError("target_sum must be positive")
    validate_raw_counts(adata)
    query = adata.copy()
    totals = np.asarray(query.X.sum(axis=1)).ravel().astype(np.float64)
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise BenchmarkDataError("every benchmark cell must have a positive finite count total")
    scale = target_sum / totals
    if sparse.issparse(query.X):
        normalized = query.X.tocsr().multiply(scale[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        query.X = normalized
    else:
        query.X = np.log1p(np.asarray(query.X, dtype=np.float64) * scale[:, None])
    query.obs = pd.DataFrame(index=query.obs_names.copy())
    query.uns.clear()
    query.obsm.clear()
    query.obsp.clear()
    return query


def load_label_map(path: str | Path = DEFAULT_LABEL_MAP_PATH) -> LabelMap:
    """Load and strictly validate the frozen broad-lineage mapping."""

    map_path = Path(path).resolve()
    try:
        payload = yaml.safe_load(map_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"could not load benchmark label map {map_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("label map root must be a mapping")
    allowed_keys = {
        "schema_version",
        "model",
        "target_classes",
        "mappings",
        "explicitly_unmapped",
        "policy",
    }
    unknown_keys = set(payload) - allowed_keys
    if unknown_keys:
        raise ValueError(f"unknown label-map keys: {sorted(unknown_keys)}")
    model = payload.get("model")
    if not isinstance(model, dict):
        raise ValueError("label map model metadata must be a mapping")
    from .annotation import MODEL_NAME, MODEL_SHA256

    if model.get("name") != MODEL_NAME or model.get("sha256") != MODEL_SHA256:
        raise ValueError("label map is not pinned to the frozen CellTypist model artifact")
    targets = tuple(payload.get("target_classes", ()))
    if targets != BROAD_CLASSES:
        raise ValueError(
            f"target_classes must exactly equal the frozen order {list(BROAD_CLASSES)}"
        )
    raw_mappings = payload.get("mappings")
    raw_unmapped = payload.get("explicitly_unmapped")
    if not isinstance(raw_mappings, dict) or not isinstance(raw_unmapped, dict):
        raise ValueError("mappings and explicitly_unmapped must be mappings")
    mappings = {str(key): str(value) for key, value in raw_mappings.items()}
    explicitly_unmapped = {
        str(key): str(reason) for key, reason in raw_unmapped.items()
    }
    if any(value not in BROAD_CLASSES for value in mappings.values()):
        raise ValueError("all mapped values must be one of the five frozen target classes")
    overlap = set(mappings) & set(explicitly_unmapped)
    if overlap:
        raise ValueError(f"labels cannot be both mapped and unmapped: {sorted(overlap)}")
    if len(mappings) + len(explicitly_unmapped) != 51:
        raise ValueError(
            "frozen label map must account for all 51 Healthy_COVID19_PBMC labels"
        )
    return LabelMap(
        schema_version=str(payload.get("schema_version")),
        model_name=str(model["name"]),
        model_sha256=str(model["sha256"]),
        target_classes=targets,
        mappings=mappings,
        explicitly_unmapped=explicitly_unmapped,
        path=map_path,
    )


def _coerce_series(values: Any, name: str) -> pd.Series:
    if isinstance(values, pd.Series):
        result = values.copy()
    else:
        result = pd.Series(values)
    result.name = name
    result.index = result.index.astype(str)
    if not result.index.is_unique:
        raise ValueError(f"{name} index must be unique")
    if result.empty:
        raise ValueError(f"{name} must not be empty")
    if result.isna().any():
        raise ValueError(f"{name} contains missing values")
    return result.astype(str)


def calculate_benchmark_metrics(
    reference_labels: pd.Series,
    predicted_labels: pd.Series,
    *,
    label_map: LabelMap | str | Path | None = None,
    confidence: pd.Series | None = None,
    runtime_seconds: float = 0.0,
) -> BenchmarkResult:
    """Calculate strict metrics where every unmapped prediction is incorrect."""

    mapping = load_label_map() if label_map is None else (
        load_label_map(label_map) if isinstance(label_map, (str, Path)) else label_map
    )
    if not isinstance(mapping, LabelMap):
        raise TypeError("label_map must be a LabelMap or path")
    truth = _coerce_series(reference_labels, "reference_label")
    raw_prediction = _coerce_series(predicted_labels, "predicted_label")
    if not truth.index.equals(raw_prediction.index):
        raise ValueError(
            "reference and prediction cell identifiers/order must match exactly; "
            "implicit reindexing is forbidden"
        )
    invalid_truth = sorted(set(truth) - set(BROAD_CLASSES))
    if invalid_truth:
        raise ValueError(f"reference labels outside frozen classes: {invalid_truth}")
    if not math.isfinite(runtime_seconds) or runtime_seconds < 0:
        raise ValueError("runtime_seconds must be finite and non-negative")

    mapped_values: list[str] = []
    statuses: list[str] = []
    reasons: list[str] = []
    for label in raw_prediction:
        if label in mapping.mappings:
            mapped_values.append(mapping.mappings[label])
            statuses.append("mapped")
            reasons.append("")
        elif label in mapping.explicitly_unmapped:
            mapped_values.append(UNMAPPED_LABEL)
            statuses.append("explicitly_unmapped")
            reasons.append(mapping.explicitly_unmapped[label])
        else:
            mapped_values.append(UNMAPPED_LABEL)
            statuses.append("unknown_model_label")
            reasons.append("label absent from the frozen version-controlled map")
    mapped = pd.Series(mapped_values, index=truth.index, name="mapped_prediction")
    covered = mapped.ne(UNMAPPED_LABEL)

    y_true = truth.to_numpy(dtype=str)
    y_pred = mapped.to_numpy(dtype=str)
    labels = list(BROAD_CLASSES)
    precision, recall, per_f1, support = precision_recall_fscore_support(
        y_true,
        y_pred,
        labels=labels,
        zero_division=0,
    )
    per_class = pd.DataFrame(
        {
            "precision": precision,
            "recall": recall,
            "f1": per_f1,
            "support": support.astype(int),
        },
        index=pd.Index(labels, name="reference_class"),
    )
    matrix_labels = labels + [UNMAPPED_LABEL]
    cm_all = confusion_matrix(y_true, y_pred, labels=matrix_labels)
    cm = pd.DataFrame(
        cm_all[: len(labels), :],
        index=pd.Index(labels, name="reference_class"),
        columns=pd.Index(matrix_labels, name="predicted_class"),
    )
    row_totals = cm.sum(axis=1).replace(0, np.nan)
    cm_normalized = cm.div(row_totals, axis=0).fillna(0.0)

    per_cell = pd.DataFrame(
        {
            "reference_label": truth,
            "predicted_label": raw_prediction,
            "mapped_prediction": mapped,
            "mapping_status": statuses,
            "mapping_reason": reasons,
            "covered": covered,
            "strict_correct": mapped.eq(truth),
        },
        index=truth.index,
    )
    if confidence is not None:
        scores = confidence.copy()
        scores.index = scores.index.astype(str)
        if not truth.index.equals(scores.index):
            raise ValueError("confidence cell identifiers/order must match exactly")
        score_values = scores.to_numpy(dtype=float)
        if not np.isfinite(score_values).all() or ((score_values < 0) | (score_values > 1)).any():
            raise ValueError("confidence values must be finite and in [0, 1]")
        per_cell["confidence"] = score_values

    metrics = {
        "n_cells": int(len(truth)),
        "coverage": float(covered.mean()),
        "n_mapped": int(covered.sum()),
        "n_unmapped": int((~covered).sum()),
        "strict_accuracy": float(mapped.eq(truth).mean()),
        "strict_macro_f1": float(
            f1_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)
        ),
        "strict_balanced_accuracy": float(
            recall_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)
        ),
        "strict_mcc": float(matthews_corrcoef(y_true, y_pred)),
        "runtime_seconds": float(runtime_seconds),
    }
    metadata = {
        "metric_policy": (
            "strict: predictions absent from or explicitly excluded by the frozen "
            "mapping count as errors; metrics are never conditioned on covered cells"
        ),
        "target_classes": labels,
        "unmapped_sentinel": UNMAPPED_LABEL,
        "label_map_path": str(mapping.path),
        "label_map_schema_version": mapping.schema_version,
        "model_name": mapping.model_name,
        "model_sha256": mapping.model_sha256,
        "limitations": list(BENCHMARK_LIMITATIONS),
    }
    return BenchmarkResult(
        metrics=metrics,
        per_class=per_class,
        confusion_matrix=cm,
        confusion_matrix_normalized=cm_normalized,
        per_cell=per_cell,
        metadata=metadata,
    )


def _extract_annotation_tables(result: Any) -> tuple[pd.Series, pd.Series | None, dict[str, Any]]:
    predictions = getattr(result, "predictions", result)
    if not isinstance(predictions, pd.DataFrame):
        raise TypeError("annotator must return AnnotationResult or a pandas DataFrame")
    if "predicted_label" not in predictions.columns:
        raise ValueError("annotator output must contain predicted_label")
    confidence = predictions["confidence"] if "confidence" in predictions.columns else None
    metadata = {
        "model": dict(getattr(result, "model_metadata", {})),
        "input_validation": dict(getattr(result, "input_validation", {})),
    }
    return predictions["predicted_label"], confidence, metadata


def run_benchmark(
    cohort: BenchmarkCohort,
    annotator: Callable[[Any], Any],
    *,
    label_map: LabelMap | str | Path | None = None,
    target_sum: float = 10_000.0,
) -> BenchmarkResult:
    """Normalize a label-free query, time annotation, and score strict outputs."""

    if list(cohort.adata.obs.columns):
        raise BenchmarkDataError(
            "benchmark AnnData contains obs columns; reference leakage protection requires none"
        )
    if cohort.adata.uns:
        raise BenchmarkDataError(
            "benchmark AnnData contains .uns metadata; rebuild the sanitized cohort"
        )
    if not pd.Index(cohort.adata.obs_names.astype(str)).equals(
        cohort.reference_labels.index
    ):
        raise BenchmarkDataError(
            "external reference labels do not exactly align with the sanitized query"
        )
    query = normalize_log1p_for_celltypist(cohort.adata, target_sum=target_sum)
    if list(query.obs.columns) or query.uns:
        raise AssertionError("reference leakage protection failed while normalizing query")

    started = time.perf_counter()
    annotation_result = annotator(query)
    elapsed = time.perf_counter() - started
    predicted, confidence, annotation_metadata = _extract_annotation_tables(
        annotation_result
    )
    result = calculate_benchmark_metrics(
        cohort.reference_labels,
        predicted,
        label_map=label_map,
        confidence=confidence,
        runtime_seconds=elapsed,
    )
    metadata = dict(result.metadata)
    metadata["cohort"] = cohort.metadata
    metadata["annotation"] = annotation_metadata
    return BenchmarkResult(
        metrics=result.metrics,
        per_class=result.per_class,
        confusion_matrix=result.confusion_matrix,
        confusion_matrix_normalized=result.confusion_matrix_normalized,
        per_cell=result.per_cell,
        metadata=metadata,
    )


def run_facs_benchmark(
    cohort: BenchmarkCohort,
    model_dir: str | Path,
    *,
    label_map: LabelMap | str | Path | None = None,
    model_url: str | None = None,
    expected_model_sha256: str | None = None,
    target_sum: float = 10_000.0,
) -> BenchmarkResult:
    """Convenience wrapper around the fixed CellTypist annotation adapter."""

    from .annotation import annotate_cells

    def fixed_annotator(query: Any) -> Any:
        return annotate_cells(
            query,
            model_dir,
            model_url=model_url,
            expected_sha256=expected_model_sha256,
            target_sum=target_sum,
        )

    return run_benchmark(
        cohort,
        fixed_annotator,
        label_map=label_map,
        target_sum=target_sum,
    )


def write_benchmark_artifacts(
    result: BenchmarkResult,
    output_dir: str | Path,
    *,
    cohort: BenchmarkCohort | None = None,
    fetched_datasets: Mapping[str, FetchedFacsDataset] | None = None,
    label_map: LabelMap | str | Path | None = None,
) -> list[Path]:
    """Persist recomputable metrics, predictions, provenance, and a heatmap."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    mapping = load_label_map() if label_map is None else (
        load_label_map(label_map) if isinstance(label_map, (str, Path)) else label_map
    )
    if not isinstance(mapping, LabelMap):
        raise TypeError("label_map must be a LabelMap or path")

    paths: list[Path] = []

    def write_json(name: str, payload: Any) -> Path:
        path = destination / name
        temporary = path.with_suffix(path.suffix + ".part")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
        temporary.replace(path)
        paths.append(path)
        return path

    write_json(
        "metrics.json",
        {
            "schema_version": "1.0",
            **result.metrics,
            "metric_policy": result.metadata["metric_policy"],
        },
    )
    write_json("benchmark_metadata.json", result.metadata)

    per_class_path = destination / "per_class_metrics.csv"
    result.per_class.to_csv(per_class_path)
    paths.append(per_class_path)
    confusion_path = destination / "confusion_matrix.csv"
    result.confusion_matrix.to_csv(confusion_path)
    paths.append(confusion_path)
    normalized_path = destination / "confusion_matrix_normalized.csv"
    result.confusion_matrix_normalized.to_csv(normalized_path)
    paths.append(normalized_path)
    predictions_path = destination / "cell_predictions.csv.gz"
    result.per_cell.to_csv(predictions_path, compression="gzip")
    paths.append(predictions_path)

    if cohort is not None:
        sample_path = destination / "sample_manifest.csv.gz"
        cohort.sample_metadata.to_csv(sample_path, compression="gzip")
        paths.append(sample_path)

    source_payload: dict[str, Any] = {
        "reference_label_design": (
            "Each class is the official FACS-enriched population named by the "
            "corresponding 10x dataset; labels were held outside annotator AnnData."
        ),
        "datasets": [],
        "limitations": list(BENCHMARK_LIMITATIONS),
    }
    if fetched_datasets:
        for dataset_id, fetched in fetched_datasets.items():
            source_manifest = json.loads(
                fetched.manifest_path.read_text(encoding="utf-8")
            )
            source_payload["datasets"].append(
                {
                    "dataset_id": dataset_id,
                    "reference_label": fetched.spec.reference_label,
                    "landing_page": fetched.spec.landing_page,
                    "download_url": fetched.spec.download_url,
                    "archive_sha256": fetched.sha256,
                    "archive_size_bytes": fetched.archive_path.stat().st_size,
                    "reported_cells": fetched.spec.reported_cells,
                    "reported_purity_percent": fetched.spec.reported_purity_percent,
                    "license": fetched.spec.license,
                    "extracted_member_sha256": source_manifest.get(
                        "extracted_member_sha256", {}
                    ),
                }
            )
    write_json("data_manifest.json", source_payload)

    frozen_map_path = destination / "label_map.yaml"
    shutil.copy2(mapping.path, frozen_map_path)
    paths.append(frozen_map_path)

    figure, axis = plt.subplots(figsize=(9, 5.5))
    sns.heatmap(
        result.confusion_matrix,
        annot=True,
        fmt="d",
        cmap="Blues",
        cbar=False,
        ax=axis,
    )
    axis.set_title("CellTypist strict five-class FACS benchmark")
    figure.tight_layout()
    figure_path = destination / "confusion_matrix.png"
    figure.savefig(figure_path, dpi=170, bbox_inches="tight")
    plt.close(figure)
    paths.append(figure_path)

    report_path = destination / "BENCHMARK_REPORT.md"
    report_lines = [
        "# Independent annotation benchmark",
        "",
        "**Status:** completed on real public FACS-enriched 10x datasets",
        "",
        "## Measured strict metrics",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for key in (
        "n_cells",
        "coverage",
        "strict_accuracy",
        "strict_macro_f1",
        "strict_balanced_accuracy",
        "strict_mcc",
        "runtime_seconds",
    ):
        value = result.metrics[key]
        formatted = str(value) if isinstance(value, int) else f"{float(value):.6f}"
        report_lines.append(f"| `{key}` | {formatted} |")
    report_lines.extend(
        [
            "",
            "## Reference-label design",
            "",
            "Each reference class is the named population in an official 10x FACS-enriched dataset. Labels are held in an external Series and are absent from the AnnData object supplied to CellTypist.",
            "",
            "## Limitations",
            "",
            *[f"- {item}" for item in BENCHMARK_LIMITATIONS],
            "",
            "Unmapped predictions count as errors in every strict metric. These results must not be described as donor-level generalization or fine-subtype accuracy.",
            "",
        ]
    )
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    paths.append(report_path)
    return paths


__all__ = [
    "BENCHMARK_LIMITATIONS",
    "BROAD_CLASSES",
    "BenchmarkCohort",
    "BenchmarkDataError",
    "BenchmarkResult",
    "DEFAULT_LABEL_MAP_PATH",
    "FACS_DATASETS",
    "FacsDatasetSpec",
    "FetchedFacsDataset",
    "LabelMap",
    "UNMAPPED_LABEL",
    "UnsafeArchiveError",
    "build_facs_cohort",
    "calculate_benchmark_metrics",
    "fetch_all_facs_datasets",
    "fetch_facs_dataset",
    "load_label_map",
    "normalize_log1p_for_celltypist",
    "run_benchmark",
    "run_facs_benchmark",
    "safe_extract_tar",
    "write_benchmark_artifacts",
]
