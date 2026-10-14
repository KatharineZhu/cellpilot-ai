"""Deterministic Scanpy analysis pipeline for PBMC 3k."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .data import validate_raw_counts
from .markers import rank_markers, top_marker_panel


def _save_figure(path: Path) -> None:
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.gcf()
    figure.savefig(path, dpi=160, bbox_inches="tight")
    plt.close("all")


def _copy_matrix(matrix: Any) -> Any:
    return matrix.copy() if hasattr(matrix, "copy") else np.array(matrix, copy=True)


def _matrix_min(matrix: Any) -> float:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    return float(values.min()) if values.size else 0.0


def _get(config: Any, name: str, default: Any) -> Any:
    return getattr(config, name, default)


def run_analysis(
    adata: Any,
    *,
    qc_config: Any,
    analysis_config: Any,
    run_dir: Path,
) -> tuple[Any, pd.DataFrame, dict[str, Any]]:
    """Run QC through marker analysis and persist real scientific artifacts."""

    import scanpy as sc

    run_dir = run_dir.resolve()
    table_dir = run_dir / "tables"
    plot_dir = run_dir / "plots"
    data_dir = run_dir / "data"
    for directory in (table_dir, plot_dir, data_dir):
        directory.mkdir(parents=True, exist_ok=True)

    raw_validation = validate_raw_counts(adata)
    adata = adata.copy()
    adata.var_names_make_unique()
    mito_prefixes = tuple(
        str(prefix).upper() for prefix in _get(qc_config, "mito_prefix", ("MT-",))
    )
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefixes)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    min_genes = int(_get(qc_config, "min_genes", 200))
    configured_max_genes = _get(qc_config, "max_genes", 2500)
    max_genes = (
        int(configured_max_genes) if configured_max_genes is not None else None
    )
    configured_max_pct_mt = _get(qc_config, "max_pct_mt", 5.0)
    max_pct_mt = (
        float(configured_max_pct_mt) if configured_max_pct_mt is not None else None
    )
    min_counts = _get(qc_config, "min_counts", None)
    max_counts = _get(qc_config, "max_counts", None)
    min_cells_per_gene = int(_get(qc_config, "min_cells_per_gene", 3))

    qc_before = adata.obs[
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    ].copy()
    qc_before["pass_min_genes"] = qc_before["n_genes_by_counts"] >= min_genes
    qc_before["pass_max_genes"] = (
        True
        if max_genes is None
        else qc_before["n_genes_by_counts"] <= max_genes
    )
    qc_before["pass_mito"] = (
        True if max_pct_mt is None else qc_before["pct_counts_mt"] < max_pct_mt
    )
    qc_before["pass_min_counts"] = (
        True if min_counts is None else qc_before["total_counts"] >= int(min_counts)
    )
    qc_before["pass_max_counts"] = (
        True if max_counts is None else qc_before["total_counts"] <= int(max_counts)
    )
    qc_before["pass_all"] = qc_before[
        [
            "pass_min_genes",
            "pass_max_genes",
            "pass_mito",
            "pass_min_counts",
            "pass_max_counts",
        ]
    ].all(axis=1)
    qc_before.index.name = "cell_id"
    qc_before.to_csv(table_dir / "qc_cells_before.csv")

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.25,
        multi_panel=True,
        show=False,
    )
    _save_figure(plot_dir / "qc_violin_before_filtering.png")

    n_cells_before, n_genes_before = int(adata.n_obs), int(adata.n_vars)
    adata = adata[qc_before["pass_all"].to_numpy()].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    if adata.n_obs < 50 or adata.n_vars < 200:
        raise ValueError(
            f"QC retained an implausibly small matrix: {adata.n_obs} cells x {adata.n_vars} genes"
        )

    adata.layers["counts"] = _copy_matrix(adata.X)
    counts_snapshot = _copy_matrix(adata.layers["counts"])
    target_sum = float(_get(analysis_config, "target_sum", 10_000.0))
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    if _matrix_min(adata.X) < 0:
        raise AssertionError("Log-normalized expression unexpectedly contains negative values")
    adata.layers["log1p"] = _copy_matrix(adata.X)
    adata.raw = adata

    # Assert that normalization did not mutate the preserved count layer.
    if sparse.issparse(counts_snapshot):
        delta = counts_snapshot != adata.layers["counts"]
        if delta.nnz:
            raise AssertionError("Preserved count layer changed during normalization")
    elif not np.array_equal(counts_snapshot, adata.layers["counts"]):
        raise AssertionError("Preserved count layer changed during normalization")

    requested_hvgs = int(_get(analysis_config, "n_top_genes", 2000))
    n_top_genes = min(requested_hvgs, max(1, adata.n_vars - 1))
    hvg_flavor = str(_get(analysis_config, "hvg_flavor", "seurat"))
    hvg_kwargs: dict[str, Any] = {
        "n_top_genes": n_top_genes,
        "flavor": hvg_flavor,
        "subset": False,
    }
    if hvg_flavor == "seurat_v3":
        hvg_kwargs["layer"] = "counts"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    sc.pl.highly_variable_genes(adata, show=False)
    _save_figure(plot_dir / "highly_variable_genes.png")

    hvg_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
    pca_input = adata[:, hvg_mask].copy()
    sc.pp.scale(pca_input, max_value=10)
    requested_pcs = int(_get(analysis_config, "n_pcs", 40))
    n_comps = min(requested_pcs, pca_input.n_obs - 1, pca_input.n_vars - 1)
    if n_comps < 2:
        raise ValueError("Insufficient dimensions for PCA")
    seed = int(_get(analysis_config, "random_seed", 0))
    sc.tl.pca(
        pca_input,
        n_comps=n_comps,
        svd_solver="arpack",
        random_state=seed,
    )
    adata.obsm["X_pca"] = np.array(pca_input.obsm["X_pca"], copy=True)
    adata.uns["pca"] = dict(pca_input.uns["pca"])
    sc.pl.pca_variance_ratio(adata, n_pcs=n_comps, log=True, show=False)
    _save_figure(plot_dir / "pca_variance_ratio.png")

    n_neighbors = min(
        int(_get(analysis_config, "n_neighbors", 10)), max(2, adata.n_obs - 1)
    )
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_comps,
        use_rep="X_pca",
        random_state=seed,
    )
    sc.tl.umap(
        adata,
        min_dist=float(_get(analysis_config, "umap_min_dist", 0.5)),
        random_state=seed,
    )
    sc.tl.leiden(
        adata,
        resolution=float(_get(analysis_config, "leiden_resolution", 0.5)),
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden",
    )
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")

    marker_top_n = int(_get(analysis_config, "marker_top_n", 25))
    markers = rank_markers(
        adata,
        groupby="leiden",
        n_genes=marker_top_n,
        method=str(_get(analysis_config, "marker_method", "wilcoxon")),
    )
    markers.to_csv(table_dir / "markers.csv", index=False)

    cluster_counts = (
        adata.obs["leiden"]
        .value_counts(sort=False)
        .rename_axis("cluster")
        .rename("n_cells")
        .reset_index()
    )
    cluster_counts["fraction"] = cluster_counts["n_cells"] / adata.n_obs
    cluster_counts.to_csv(table_dir / "cluster_counts.csv", index=False)

    umap_table = pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["umap_1", "umap_2"],
    )
    umap_table.insert(0, "cluster", adata.obs["leiden"].astype(str).to_numpy())
    umap_table.index.name = "cell_id"
    umap_table.to_csv(table_dir / "umap_coordinates.csv")

    sc.pl.umap(
        adata,
        color="leiden",
        legend_loc="on data",
        frameon=False,
        title="PBMC 3k — Leiden clusters",
        show=False,
    )
    _save_figure(plot_dir / "umap_clusters.png")

    panel = top_marker_panel(markers, per_cluster=3)
    if panel:
        sc.pl.dotplot(
            adata,
            panel,
            groupby="leiden",
            use_raw=True,
            standard_scale="var",
            show=False,
        )
        _save_figure(plot_dir / "marker_dotplot.png")

    processed_path = data_dir / "pbmc3k_processed.h5ad"
    adata.write_h5ad(processed_path, compression="gzip")

    summary: dict[str, Any] = {
        "raw_validation": raw_validation,
        "cells_before_qc": n_cells_before,
        "genes_before_qc": n_genes_before,
        "cells_after_qc": int(adata.n_obs),
        "genes_after_qc": int(adata.n_vars),
        "cells_removed": int(n_cells_before - adata.n_obs),
        "genes_removed": int(n_genes_before - adata.n_vars),
        "n_highly_variable_genes": int(hvg_mask.sum()),
        "n_pcs": int(n_comps),
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "random_seed": seed,
        "normalization": f"normalize_total(target_sum={target_sum}) then log1p exactly once",
        "marker_test": "Wilcoxon rank-sum on full-gene log-normalized adata.raw; BH-adjusted p-values",
        "limitations": [
            "PBMC 3k is a Cell Ranger filtered matrix, not raw droplets or FASTQ.",
            "The dataset represents one donor; marker tests are descriptive cell-level analyses, not donor-level inference.",
            "This workflow does not model ambient RNA, doublets, batches, or condition effects.",
        ],
    }
    (run_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return adata, markers, summary
