"""Read-only Streamlit dashboard for completed CellPilot run artifacts."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from cellpilot.provenance import verify_manifest


def _argument_value(flag: str) -> str | None:
    if flag in sys.argv:
        index = sys.argv.index(flag)
        if index + 1 < len(sys.argv):
            return sys.argv[index + 1]
    return None


def _default_run_dir() -> str:
    explicit = os.getenv("CELLPILOT_RUN_DIR") or _argument_value("--run-dir")
    if explicit:
        return explicit
    pointer = Path("results/LATEST_RUN.txt")
    if pointer.exists():
        return pointer.read_text(encoding="utf-8").strip()
    return "results"


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        st.error(f"Could not read {path.name}: {exc}")
        return None
    if not isinstance(payload, dict):
        st.error(f"Could not read {path.name}: the JSON root must be an object.")
        return None
    return payload


@st.cache_data(show_spinner=False)
def _verify_run_manifest(run_dir: str, manifest_mtime_ns: int) -> dict[str, Any]:
    """Verify a manifest once per manifest revision for a responsive UI."""

    del manifest_mtime_ns  # Included in the cache key to invalidate changed manifests.
    try:
        report = verify_manifest(Path(run_dir))
    except (OSError, ValueError) as exc:
        return {
            "valid": False,
            "checked_inputs": 0,
            "checked_outputs": 0,
            "issues": [
                {
                    "kind": "verification_error",
                    "path": "manifest.json",
                    "error": str(exc),
                }
            ],
        }
    return {
        "valid": report.valid,
        "checked_inputs": report.checked_inputs,
        "checked_outputs": report.checked_outputs,
        "issues": [issue.model_dump(mode="json") for issue in report.issues],
    }


def _show_image(path: Path, caption: str) -> None:
    if path.exists():
        st.image(str(path), caption=caption, width="stretch")
    else:
        st.info(f"Artifact not available: {path.name}")


def render_dashboard() -> None:
    st.set_page_config(page_title="CellPilot AI", page_icon="🧬", layout="wide")
    st.title("CellPilot AI")
    st.caption(
        "Read-only view of generated artifacts. The dashboard does not rerun analysis or invent missing values."
    )

    run_input = st.sidebar.text_input("Run directory", value=_default_run_dir())
    run_dir = Path(run_input).expanduser().resolve()
    pages = [
        "Overview",
        "QC",
        "Clusters & markers",
        "CellTypist annotation",
        "Evidence-grounded interpretation",
        "Benchmark",
        "Provenance & limitations",
    ]
    page = st.sidebar.radio("Page", pages)

    if not run_dir.is_dir():
        st.error(f"Run directory does not exist: {run_dir}")
        st.stop()

    summary = _read_json(run_dir / "analysis_summary.json")
    stage_status = _read_json(run_dir / "stage_status.json")
    manifest = _read_json(run_dir / "manifest.json")

    if page == "Overview":
        if summary is None:
            st.warning("No analysis_summary.json is present; this run is incomplete.")
        else:
            columns = st.columns(4)
            columns[0].metric("Cells after QC", summary.get("cells_after_qc", "unknown"))
            columns[1].metric("Genes after QC", summary.get("genes_after_qc", "unknown"))
            columns[2].metric("Leiden clusters", summary.get("n_clusters", "unknown"))
            columns[3].metric("HVGs", summary.get("n_highly_variable_genes", "unknown"))
        left, right = st.columns(2)
        with left:
            _show_image(run_dir / "plots" / "umap_clusters.png", "Leiden clusters")
        with right:
            _show_image(run_dir / "plots" / "marker_dotplot.png", "Top marker genes")
        if stage_status:
            st.subheader("Stage status")
            st.json(stage_status)

    elif page == "QC":
        _show_image(
            run_dir / "plots" / "qc_violin_before_filtering.png",
            "QC distributions before filtering",
        )
        qc_path = run_dir / "tables" / "qc_cells_before.csv"
        if qc_path.exists():
            qc = pd.read_csv(qc_path)
            pass_rate = float(qc["pass_all"].mean()) if "pass_all" in qc else None
            if pass_rate is not None:
                st.metric("Cells passing all QC rules", f"{pass_rate:.1%}")
            numeric_summary = (
                qc.select_dtypes(include="number")
                .describe()
                .transpose()
                .rename_axis("metric")
                .reset_index()
            )
            st.dataframe(numeric_summary, width="stretch")
        else:
            st.info("QC table is unavailable.")

    elif page == "Clusters & markers":
        _show_image(run_dir / "plots" / "umap_clusters.png", "UMAP by Leiden cluster")
        _show_image(run_dir / "plots" / "marker_dotplot.png", "Top-ranked markers")
        markers_path = run_dir / "tables" / "markers.csv"
        if markers_path.exists():
            markers = pd.read_csv(markers_path)
            clusters = sorted(markers["cluster"].astype(str).unique())
            selected = st.selectbox("Cluster", clusters)
            st.dataframe(
                markers.loc[markers["cluster"].astype(str) == selected],
                width="stretch",
            )
        else:
            st.info("Marker table is unavailable.")

    elif page == "CellTypist annotation":
        metadata = _read_json(run_dir / "annotation" / "celltypist_metadata.json")
        if metadata:
            st.json(metadata)
        cells_path = run_dir / "tables" / "celltypist_cells.csv"
        consensus_path = run_dir / "tables" / "celltypist_cluster_consensus.csv"
        if consensus_path.exists():
            st.subheader("Cluster consensus")
            st.dataframe(pd.read_csv(consensus_path), width="stretch")
        if cells_path.exists():
            cells = pd.read_csv(cells_path)
            st.subheader("Per-cell predictions")
            st.dataframe(cells.head(500), width="stretch")
            st.caption(f"Showing up to 500 of {len(cells):,} predictions.")
            _show_image(
                run_dir / "plots" / "umap_celltypist.png", "UMAP by CellTypist label"
            )
        elif not metadata:
            st.warning("CellTypist output is absent; no annotation claim is shown.")

    elif page == "Evidence-grounded interpretation":
        interpretation = _read_json(
            run_dir / "interpretation" / "interpretation.json"
        )
        if interpretation:
            st.warning(
                "These labels are evidence-bounded hypotheses, not experimental validation."
            )
            st.json(interpretation)
        else:
            st.info("Interpretation was not generated for this run.")

    elif page == "Benchmark":
        benchmark_dir = run_dir / "benchmark"
        metrics = _read_json(benchmark_dir / "metrics.json")
        if metrics:
            columns = st.columns(4)
            columns[0].metric("Cells benchmarked", metrics.get("n_cells", "unknown"))
            coverage = metrics.get("coverage")
            accuracy = metrics.get("strict_accuracy")
            macro_f1 = metrics.get("strict_macro_f1")
            columns[1].metric(
                "Mapped coverage",
                f"{float(coverage):.1%}" if isinstance(coverage, (int, float)) else "unknown",
            )
            columns[2].metric(
                "Strict accuracy",
                f"{float(accuracy):.1%}" if isinstance(accuracy, (int, float)) else "unknown",
            )
            columns[3].metric(
                "Strict macro-F1",
                f"{float(macro_f1):.3f}"
                if isinstance(macro_f1, (int, float))
                else "unknown",
            )
            st.json(metrics)
            _show_image(
                benchmark_dir / "confusion_matrix.png", "Strict broad-class confusion matrix"
            )
            per_class = benchmark_dir / "per_class_metrics.csv"
            if per_class.exists():
                st.dataframe(pd.read_csv(per_class), width="stretch")
        else:
            st.info("No completed benchmark is attached to this run.")

    else:
        if manifest:
            verification = _verify_run_manifest(
                str(run_dir), (run_dir / "manifest.json").stat().st_mtime_ns
            )
            if verification["valid"]:
                st.success(
                    "Manifest verification passed for "
                    f"{verification['checked_inputs']} recorded inputs and "
                    f"{verification['checked_outputs']} recorded outputs."
                )
            else:
                st.error("Manifest verification failed; do not trust this run as intact.")
                st.json(verification["issues"])
            st.subheader("Artifact manifest")
            st.json(manifest)
        else:
            st.warning("No manifest is present; artifact integrity has not been recorded.")
        if summary and summary.get("limitations"):
            st.subheader("Declared limitations")
            for limitation in summary["limitations"]:
                st.markdown(f"- {limitation}")


def main() -> None:
    render_dashboard()


if __name__ == "__main__":
    main()
