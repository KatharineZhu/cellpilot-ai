from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import cellpilot.annotation as annotation


def _normalized_adata(n_cells: int = 3, n_genes: int = 1_200) -> ad.AnnData:
    counts = np.zeros((n_cells, n_genes), dtype=np.float64)
    span = min(20, n_genes)
    for cell in range(n_cells):
        counts[cell, :span] = np.arange(1, span + 1)
    counts *= 10_000.0 / counts.sum(axis=1, keepdims=True)
    result = ad.AnnData(sparse.csr_matrix(np.log1p(counts)))
    result.obs_names = [f"cell_{index}" for index in range(n_cells)]
    result.var_names = [f"GENE_{index}" for index in range(n_genes)]
    return result


def test_validate_celltypist_input_accepts_full_gene_log1p_matrix() -> None:
    result = annotation.validate_celltypist_input(_normalized_adata())
    assert result["n_cells"] == 3
    assert result["n_genes"] == 1_200
    assert result["inverse_log_total_min"] == pytest.approx(10_000.0)


@pytest.mark.parametrize("failure", ["negative", "raw_counts", "subsetted"])
def test_validate_celltypist_input_fails_closed(failure: str) -> None:
    query = _normalized_adata()
    if failure == "negative":
        query.X = query.X.toarray()
        query.X[0, 0] = -0.1
        expected = "negative"
    elif failure == "raw_counts":
        query.X = sparse.csr_matrix(np.ones(query.shape))
        expected = "not log1p-normalized"
    else:
        query.raw = query
        query = query[:, :1_000].copy()
        expected = "feature-subsetted"
    with pytest.raises(annotation.AnnotationInputError, match=expected):
        annotation.validate_celltypist_input(query, min_genes=500)


def test_ensure_model_uses_only_fixed_workspace_path_and_records_digest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"synthetic fixed model bytes"
    digest = hashlib.sha256(payload).hexdigest()
    monkeypatch.setattr(annotation, "MODEL_SHA256", digest)
    seen: list[tuple[str, int]] = []

    class FakeResponse:
        def __enter__(self) -> "FakeResponse":
            return self

        def __exit__(self, *args: object) -> None:
            return None

        def iter_content(self, chunk_size: int):
            assert chunk_size == 1024 * 1024
            yield payload

    def fake_open(url: str, timeout: int) -> FakeResponse:
        seen.append((url, timeout))
        return FakeResponse()

    monkeypatch.setattr(annotation, "_open_fixed_model", fake_open)
    model_dir = tmp_path / "models"
    artifact = annotation.ensure_model(model_dir)

    assert seen == [(annotation.MODEL_URL, 120)]
    assert artifact.path == model_dir.resolve() / annotation.MODEL_NAME
    assert artifact.path.read_bytes() == payload
    assert artifact.sha256 == digest
    metadata = json.loads(artifact.metadata_path.read_text(encoding="utf-8"))
    assert metadata["sha256"] == digest
    assert metadata["url"] == annotation.MODEL_URL
    assert metadata["policy"].startswith("fixed URL")

    with pytest.raises(annotation.ModelArtifactError, match="only the frozen model URL"):
        annotation.ensure_model(model_dir, model_url="https://example.invalid/model.pkl")


def test_annotate_cells_forces_primary_per_cell_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = _normalized_adata(n_cells=2, n_genes=20)
    artifact = annotation.ModelArtifact(
        path=tmp_path / annotation.MODEL_NAME,
        metadata_path=tmp_path / "model.metadata.json",
        name=annotation.MODEL_NAME,
        url=annotation.MODEL_URL,
        sha256=annotation.MODEL_SHA256,
        size_bytes=123,
        version="v1",
        source=annotation.MODEL_SOURCE,
        model_date=annotation.MODEL_DATE,
        declared_cell_types=51,
        declared_features=20,
        downloaded=False,
    )
    monkeypatch.setattr(annotation, "ensure_model", lambda *args, **kwargs: artifact)

    class FakeModel:
        features = np.asarray(query.var_names)
        cell_types = np.array(["B_naive", "CD14_mono"])

    calls: list[dict[str, object]] = []

    def fake_annotate(adata: ad.AnnData, **kwargs: object) -> SimpleNamespace:
        calls.append(kwargs)
        index = pd.Index(adata.obs_names.astype(str))
        return SimpleNamespace(
            predicted_labels=pd.DataFrame(
                {"predicted_labels": ["B_naive", "CD14_mono"]}, index=index
            ),
            probability_matrix=pd.DataFrame(
                [[0.9, 0.1], [0.2, 0.8]],
                index=index,
                columns=["B_naive", "CD14_mono"],
            ),
        )

    fake_celltypist = SimpleNamespace(
        __version__="test",
        models=SimpleNamespace(Model=SimpleNamespace(load=lambda model: FakeModel())),
        annotate=fake_annotate,
    )
    monkeypatch.setattr(annotation, "_load_celltypist", lambda model_dir: fake_celltypist)

    result = annotation.annotate_cells(
        query,
        tmp_path,
        min_genes=10,
        min_model_feature_fraction=0.5,
    )
    assert calls[0]["mode"] == "best match"
    assert calls[0]["majority_voting"] is False
    assert result.predictions["predicted_label"].tolist() == [
        "B_naive",
        "CD14_mono",
    ]
    assert result.predictions["confidence"].tolist() == pytest.approx([0.9, 0.8])
    assert result.model_metadata["fallback_labels_used"] is False


def test_annotate_cells_never_substitutes_labels_on_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = _normalized_adata(n_cells=1, n_genes=10)
    artifact = SimpleNamespace(path=tmp_path / annotation.MODEL_NAME)
    monkeypatch.setattr(annotation, "ensure_model", lambda *args, **kwargs: artifact)
    loaded_model = SimpleNamespace(
        features=np.asarray(query.var_names), cell_types=np.array(["B_naive"])
    )

    def explode(*args: object, **kwargs: object) -> None:
        raise RuntimeError("classifier exploded")

    fake_celltypist = SimpleNamespace(
        models=SimpleNamespace(
            Model=SimpleNamespace(load=lambda **kwargs: loaded_model)
        ),
        annotate=explode,
    )
    monkeypatch.setattr(annotation, "_load_celltypist", lambda model_dir: fake_celltypist)

    with pytest.raises(annotation.AnnotationExecutionError, match="no labels were substituted"):
        annotation.annotate_cells(
            query,
            tmp_path,
            min_genes=5,
            min_model_feature_fraction=0.5,
        )
