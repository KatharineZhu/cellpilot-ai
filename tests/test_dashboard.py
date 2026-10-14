from __future__ import annotations

import base64
import json
from pathlib import Path

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from cellpilot.provenance import write_manifest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_PATH = PROJECT_ROOT / "src" / "cellpilot" / "dashboard.py"
ONE_PIXEL_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUB"
    "AScY42YAAAAASUVORK5CYII="
)


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_image(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(ONE_PIXEL_PNG)


def _complete_run(tmp_path: Path) -> Path:
    run_dir = tmp_path / "complete-run"
    tables = run_dir / "tables"
    benchmark = run_dir / "benchmark"
    tables.mkdir(parents=True)
    benchmark.mkdir(parents=True)

    _write_json(
        run_dir / "analysis_summary.json",
        {
            "cells_after_qc": 8,
            "genes_after_qc": 1200,
            "n_clusters": 2,
            "n_highly_variable_genes": 500,
            "limitations": ["Synthetic dashboard fixture; not a biological result."],
        },
    )
    _write_json(
        run_dir / "stage_status.json",
        {"stages": {"analysis": {"status": "completed"}}},
    )
    pd.DataFrame(
        {
            "cell_id": ["cell_0", "cell_1"],
            "n_genes_by_counts": [700, 800],
            "total_counts": [1500, 1800],
            "pct_counts_mt": [2.0, 3.0],
            "pass_all": [True, False],
        }
    ).to_csv(tables / "qc_cells_before.csv", index=False)
    pd.DataFrame(
        {
            "cluster": [0, 0, 1, 1],
            "rank": [1, 2, 1, 2],
            "names": ["MS4A1", "CD79A", "NKG7", "GNLY"],
            "scores": [8.0, 7.0, 9.0, 8.0],
        }
    ).to_csv(tables / "markers.csv", index=False)
    pd.DataFrame(
        {
            "cluster": [0, 1],
            "predicted_label": ["B_naive", "NK_56hi"],
            "n_cells": [4, 4],
            "mean_confidence": [0.91, 0.89],
            "fraction": [1.0, 1.0],
        }
    ).to_csv(tables / "celltypist_cluster_consensus.csv", index=False)
    pd.DataFrame(
        {
            "cell_id": [f"cell_{index}" for index in range(8)],
            "predicted_label": ["B_naive"] * 4 + ["NK_56hi"] * 4,
            "confidence": [0.9] * 8,
        }
    ).to_csv(tables / "celltypist_cells.csv", index=False)
    _write_json(
        run_dir / "annotation" / "celltypist_metadata.json",
        {"status": "completed", "n_predictions": 8},
    )
    _write_json(
        run_dir / "interpretation" / "interpretation.json",
        {
            "metadata": {
                "backend": "deterministic",
                "validated": True,
                "external_api_called": False,
            },
            "interpretations": [],
        },
    )
    _write_json(
        benchmark / "metrics.json",
        {
            "schema_version": "1.0",
            "n_cells": 10,
            "coverage": 0.9,
            "strict_accuracy": 0.8,
            "strict_macro_f1": 0.75,
            "runtime_seconds": 1.25,
        },
    )
    pd.DataFrame(
        {
            "reference_class": ["B cell", "NK cell"],
            "precision": [0.8, 0.7],
            "recall": [0.8, 0.7],
            "f1": [0.8, 0.7],
            "support": [5, 5],
        }
    ).to_csv(benchmark / "per_class_metrics.csv", index=False)

    for relative in (
        "plots/umap_clusters.png",
        "plots/marker_dotplot.png",
        "plots/qc_violin_before_filtering.png",
        "plots/umap_celltypist.png",
        "benchmark/confusion_matrix.png",
    ):
        _write_image(run_dir / relative)

    write_manifest(
        run_dir,
        project_root=tmp_path,
        packages=(),
        command=["cellpilot", "all"],
    )
    return run_dir


def _load_app(run_dir: Path, monkeypatch: pytest.MonkeyPatch) -> AppTest:
    monkeypatch.setenv("CELLPILOT_RUN_DIR", str(run_dir))
    return AppTest.from_file(str(DASHBOARD_PATH), default_timeout=20).run()


def test_dashboard_renders_complete_run_across_all_pages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    app = _load_app(_complete_run(tmp_path), monkeypatch)

    assert not app.exception
    assert [(metric.label, metric.value) for metric in app.metric] == [
        ("Cells after QC", "8"),
        ("Genes after QC", "1200"),
        ("Leiden clusters", "2"),
        ("HVGs", "500"),
    ]

    app.sidebar.radio[0].set_value("QC").run()
    assert not app.exception
    assert [(metric.label, metric.value) for metric in app.metric] == [
        ("Cells passing all QC rules", "50.0%")
    ]
    assert app.dataframe[0].value["metric"].tolist() == [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
    ]

    app.sidebar.radio[0].set_value("Clusters & markers").run()
    assert not app.exception
    assert app.selectbox[0].options == ["0", "1"]
    assert set(app.dataframe[0].value["cluster"].astype(str)) == {"0"}

    app.sidebar.radio[0].set_value("CellTypist annotation").run()
    assert not app.exception
    assert len(app.dataframe) == 2
    assert app.dataframe[1].value.shape == (8, 3)
    assert any("8 predictions" in caption.value for caption in app.caption)

    app.sidebar.radio[0].set_value("Evidence-grounded interpretation").run()
    assert not app.exception
    assert any("not experimental validation" in item.value for item in app.warning)

    app.sidebar.radio[0].set_value("Benchmark").run()
    assert not app.exception
    assert [(metric.label, metric.value) for metric in app.metric] == [
        ("Cells benchmarked", "10"),
        ("Mapped coverage", "90.0%"),
        ("Strict accuracy", "80.0%"),
        ("Strict macro-F1", "0.750"),
    ]
    assert len(app.dataframe) == 1

    app.sidebar.radio[0].set_value("Provenance & limitations").run()
    assert not app.exception
    assert any("Manifest verification passed" in item.value for item in app.success)
    assert any("Synthetic dashboard fixture" in item.value for item in app.markdown)


def test_dashboard_reports_malformed_json_without_crashing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "malformed-run"
    run_dir.mkdir()
    (run_dir / "analysis_summary.json").write_text("{not-json", encoding="utf-8")

    app = _load_app(run_dir, monkeypatch)

    assert not app.exception
    assert any("Could not read analysis_summary.json" in item.value for item in app.error)
    assert any("this run is incomplete" in item.value for item in app.warning)


def test_dashboard_reports_invalid_manifest_schema_without_crashing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "invalid-manifest-run"
    _write_json(run_dir / "manifest.json", {"schema_version": "1.0"})

    app = _load_app(run_dir, monkeypatch)
    app.sidebar.radio[0].set_value("Provenance & limitations").run()

    assert not app.exception
    assert any("Manifest verification failed" in item.value for item in app.error)
    assert any("verification_error" in item.value for item in app.json)


@pytest.mark.integration
def test_dashboard_reads_the_verified_pbmc3k_run(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run_dir = PROJECT_ROOT / "results" / "pbmc3k" / "verified_run"
    if not (run_dir / "analysis_summary.json").is_file():
        pytest.skip("verified PBMC3k artifacts are not present in this checkout")

    app = _load_app(run_dir, monkeypatch)
    assert not app.exception
    assert [(metric.label, metric.value) for metric in app.metric] == [
        ("Cells after QC", "2638"),
        ("Genes after QC", "13656"),
        ("Leiden clusters", "6"),
        ("HVGs", "2000"),
    ]

    for page in (
        "QC",
        "Clusters & markers",
        "CellTypist annotation",
        "Evidence-grounded interpretation",
        "Benchmark",
        "Provenance & limitations",
    ):
        app.sidebar.radio[0].set_value(page).run()
        assert not app.exception, page

    app.sidebar.radio[0].set_value("Benchmark").run()
    benchmark_metrics = {metric.label: metric.value for metric in app.metric}
    assert benchmark_metrics == {
        "Cells benchmarked": "2500",
        "Mapped coverage": "96.6%",
        "Strict accuracy": "93.3%",
        "Strict macro-F1": "0.948",
    }

    app.sidebar.radio[0].set_value("Provenance & limitations").run()
    assert any("Manifest verification passed" in item.value for item in app.success)
