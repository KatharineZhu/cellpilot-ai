from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from cellpilot.markers import rank_markers


def test_rank_markers_uses_full_log_matrix() -> None:
    rng = np.random.default_rng(7)
    matrix = rng.poisson(1, size=(60, 8)).astype(float)
    matrix[:30, 0] += 12
    matrix[30:, 1] += 12
    adata = ad.AnnData(matrix)
    adata.var_names = ["GROUP_A", "GROUP_B", "G2", "G3", "G4", "G5", "G6", "G7"]
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)
    adata.raw = adata

    markers = rank_markers(adata, n_genes=3)

    top = markers.loc[markers["rank"] == 1].set_index("cluster")["names"].to_dict()
    assert top == {"0": "GROUP_A", "1": "GROUP_B"}
    assert (markers["pvals_adj"] >= 0).all()
    assert (markers["pvals_adj"] <= 1).all()
