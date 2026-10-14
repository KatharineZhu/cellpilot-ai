"""Marker-gene analysis on the preserved log-normalized expression matrix."""

from __future__ import annotations

from typing import Any

import pandas as pd


def rank_markers(
    adata: Any,
    *,
    groupby: str = "leiden",
    n_genes: int = 25,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """Rank markers with BH correction, using ``adata.raw`` by construction."""

    import scanpy as sc

    if adata.raw is None:
        raise ValueError("adata.raw must contain the log-normalized full-gene matrix")
    if groupby not in adata.obs:
        raise KeyError(f"Missing grouping column: {groupby}")
    if not isinstance(adata.obs[groupby].dtype, pd.CategoricalDtype):
        adata.obs[groupby] = adata.obs[groupby].astype("category")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        use_raw=True,
        corr_method="benjamini-hochberg",
        pts=True,
        key_added="rank_genes_groups",
    )
    frames: list[pd.DataFrame] = []
    for group in adata.obs[groupby].cat.categories:
        frame = sc.get.rank_genes_groups_df(
            adata,
            group=str(group),
            key="rank_genes_groups",
        ).head(n_genes)
        frame.insert(0, "cluster", str(group))
        frame.insert(1, "rank", range(1, len(frame) + 1))
        frames.append(frame)
    markers = pd.concat(frames, ignore_index=True)
    preferred = [
        "cluster",
        "rank",
        "names",
        "scores",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "pct_nz_group",
        "pct_nz_reference",
    ]
    return markers[[column for column in preferred if column in markers.columns]]


def top_marker_panel(markers: pd.DataFrame, *, per_cluster: int = 3) -> list[str]:
    """Return an ordered, de-duplicated plotting panel from a marker table."""

    selected = markers.loc[markers["rank"] <= per_cluster, "names"].astype(str)
    return list(dict.fromkeys(selected.tolist()))
