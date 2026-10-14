from __future__ import annotations

import io
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellpilot.benchmark import (
    BROAD_CLASSES,
    BenchmarkDataError,
    FACS_DATASETS,
    UNMAPPED_LABEL,
    UnsafeArchiveError,
    build_facs_cohort,
    calculate_benchmark_metrics,
    fetch_facs_dataset,
    load_label_map,
    run_benchmark,
    safe_extract_tar,
)


def _synthetic_datasets(n_cells: int = 4, n_genes: int = 8) -> dict[str, ad.AnnData]:
    datasets: dict[str, ad.AnnData] = {}
    for dataset_position, dataset_id in enumerate(FACS_DATASETS):
        matrix = np.zeros((n_cells, n_genes), dtype=np.float32)
        for cell in range(n_cells):
            matrix[cell, (cell + dataset_position) % n_genes] = cell + 1
            matrix[cell, (cell + dataset_position + 1) % n_genes] = 2
        value = ad.AnnData(sparse.csr_matrix(matrix))
        value.obs_names = [f"source_{dataset_position}_{cell}" for cell in range(n_cells)]
        value.var_names = [f"GENE_{gene}" for gene in range(n_genes)]
        # Deliberately include source metadata that must not survive sanitization.
        value.obs["source_reference_label"] = FACS_DATASETS[dataset_id].reference_label
        value.uns["source_reference_label"] = FACS_DATASETS[dataset_id].reference_label
        datasets[dataset_id] = value
    return datasets


def test_build_facs_cohort_is_deterministic_and_keeps_labels_external() -> None:
    datasets = _synthetic_datasets()
    first = build_facs_cohort(datasets, max_cells_per_class=2, random_seed=7)
    second = build_facs_cohort(datasets, max_cells_per_class=2, random_seed=7)

    assert first.adata.shape == (10, 8)
    assert first.sample_metadata["original_cell_id"].tolist() == second.sample_metadata[
        "original_cell_id"
    ].tolist()
    assert list(first.adata.obs.columns) == []
    assert first.adata.uns == {}
    assert first.reference_labels.index.equals(first.adata.obs_names)
    assert set(first.reference_labels) == set(BROAD_CLASSES)
    assert all(
        label not in cell_id
        for label in BROAD_CLASSES
        for cell_id in first.adata.obs_names
    )


def test_run_benchmark_does_not_expose_reference_labels_to_annotator() -> None:
    cohort = build_facs_cohort(
        _synthetic_datasets(n_cells=2), max_cells_per_class=1, random_seed=11
    )
    raw_predictions = {
        "B cell": "B_naive",
        "Monocyte": "CD14_mono",
        "CD4+ T cell": "CD4.Naive",
        "CD8+ T cell": "CD8.Naive",
        "NK cell": "NK_56hi",
    }

    def label_blind_annotator(query: ad.AnnData) -> pd.DataFrame:
        assert list(query.obs.columns) == []
        assert query.uns == {}
        assert query.obs_names.tolist() == cohort.adata.obs_names.tolist()
        # The test supplies outcomes by row solely to exercise scoring.  The
        # annotator callable receives no reference Series or sample metadata.
        return pd.DataFrame(
            {
                "predicted_label": [raw_predictions[label] for label in BROAD_CLASSES],
                "confidence": [0.9] * len(BROAD_CLASSES),
            },
            index=query.obs_names,
        )

    result = run_benchmark(cohort, label_blind_annotator)
    assert result.metrics["coverage"] == 1.0
    assert result.metrics["strict_macro_f1"] == 1.0
    assert result.metrics["strict_balanced_accuracy"] == 1.0
    assert result.metrics["strict_mcc"] == 1.0
    assert result.metadata["annotation"] == {
        "model": {},
        "input_validation": {},
    }


def test_strict_metrics_count_unmapped_predictions_as_errors() -> None:
    index = pd.Index([f"cell_{value}" for value in range(5)])
    truth = pd.Series(BROAD_CLASSES, index=index)
    prediction = pd.Series(
        ["B_naive", "CD14_mono", "CD4.Naive", "CD8.Naive", "pDC"],
        index=index,
    )
    result = calculate_benchmark_metrics(
        truth,
        prediction,
        label_map=load_label_map(),
        runtime_seconds=1.25,
    )

    assert result.metrics["coverage"] == pytest.approx(0.8)
    assert result.metrics["strict_accuracy"] == pytest.approx(0.8)
    assert result.metrics["strict_macro_f1"] == pytest.approx(0.8)
    assert result.metrics["strict_balanced_accuracy"] == pytest.approx(0.8)
    assert result.metrics["runtime_seconds"] == pytest.approx(1.25)
    assert result.confusion_matrix.loc["NK cell", UNMAPPED_LABEL] == 1
    assert result.per_class.loc["NK cell", "recall"] == 0
    assert result.per_cell.loc["cell_4", "mapping_status"] == "explicitly_unmapped"


def test_metrics_forbid_implicit_cell_reindexing() -> None:
    truth = pd.Series(BROAD_CLASSES, index=[f"cell_{i}" for i in range(5)])
    prediction = pd.Series(
        ["B_naive", "CD14_mono", "CD4.Naive", "CD8.Naive", "NK_56hi"],
        index=list(reversed(truth.index)),
    )
    with pytest.raises(ValueError, match="implicit reindexing is forbidden"):
        calculate_benchmark_metrics(truth, prediction)


def test_safe_extract_tar_rejects_path_traversal(tmp_path: Path) -> None:
    archive = tmp_path / "malicious.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        payload = b"do not write outside"
        member = tarfile.TarInfo("../escaped.txt")
        member.size = len(payload)
        bundle.addfile(member, io.BytesIO(payload))

    destination = tmp_path / "extract"
    with pytest.raises(UnsafeArchiveError, match="unsafe path"):
        safe_extract_tar(archive, destination)
    assert not (tmp_path / "escaped.txt").exists()


def test_facs_fetch_rejects_cached_archive_that_does_not_match_frozen_hash(
    tmp_path: Path,
) -> None:
    archive = tmp_path / "archives" / "b_cells.tar.gz"
    archive.parent.mkdir(parents=True)
    archive.write_bytes(b"not the frozen official archive")

    with pytest.raises(BenchmarkDataError, match="archive checksum mismatch"):
        fetch_facs_dataset("b_cells", tmp_path)
