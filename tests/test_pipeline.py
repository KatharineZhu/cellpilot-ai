from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
from scipy import sparse

from cellpilot.pipeline import run_analysis


def test_pipeline_preserves_counts_and_produces_artifacts(tmp_path) -> None:
    rng = np.random.default_rng(42)
    matrix = rng.poisson(1.2, size=(90, 320)).astype(np.float32)
    matrix[:45, :20] += rng.poisson(8, size=(45, 20))
    matrix[45:, 20:40] += rng.poisson(8, size=(45, 20))
    adata = ad.AnnData(sparse.csr_matrix(matrix))
    adata.obs_names = [f"cell_{index}" for index in range(adata.n_obs)]
    adata.var_names = [f"GENE_{index}" for index in range(adata.n_vars)]
    original = adata.X.copy()

    qc = SimpleNamespace(
        min_genes=1,
        max_genes=None,
        max_pct_mt=100.0,
        min_cells_per_gene=1,
    )
    analysis = SimpleNamespace(
        target_sum=10_000,
        n_top_genes=120,
        n_pcs=12,
        n_neighbors=8,
        leiden_resolution=0.8,
        random_seed=11,
        marker_top_n=5,
        marker_method="wilcoxon",
    )

    result, markers, summary = run_analysis(
        adata,
        qc_config=qc,
        analysis_config=analysis,
        run_dir=tmp_path,
    )

    assert result.n_obs == adata.n_obs
    assert (result.layers["counts"] != original).nnz == 0
    assert result.raw is not None
    assert float(result.X.min()) >= 0
    assert summary["n_clusters"] >= 2
    assert not markers.empty
    for relative in (
        "data/pbmc3k_processed.h5ad",
        "tables/markers.csv",
        "tables/umap_coordinates.csv",
        "plots/qc_violin_before_filtering.png",
        "plots/umap_clusters.png",
        "plots/marker_dotplot.png",
        "analysis_summary.json",
    ):
        path = tmp_path / relative
        assert path.exists() and path.stat().st_size > 0
