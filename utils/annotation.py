#!/usr/bin/env python3
"""
Cell type annotation utilities for single-cell RNA-seq analysis
Handles marker gene analysis and cell type assignment
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from utils.processing import choose_leiden_resolution

# Module-level constants: single sources of truth
MARKER_GENES = {  # Try to find top 30 markers
    # General neuron/excitatory
    "Neuron": ["Snap25", "Rbfox3", "Syp"],
    "Excit": ["Slc17a7", "Camk2a", "Satb2"],
    # Excitatory layer-specific markers
    "ExN_L2-4": ["Cux1", "Cux2", "Satb2"],
    # Include both names for Ctip2/Bcl11b
    "ExN_L5": ["Bcl11b", "Ctip2", "Fezf2"],
    "ExN_L6": ["Tbr1", "Sox5"],
    "ExN_L6b": ["Ctgf"],
    # Inhibitory (generic + subclasses)
    "Inhib": ["Gad1", "Gad2", "Slc6a1"],
    "InN_SST": ["Sst", "Npy", "Chodl"],
    "InN_VIP": ["Vip", "Cck", "Calb2"],
    "InN_PVALB": ["Pvalb", "Gabra1", "Reln"],
    # Glia and vascular
    "Astro": ["Slc1a2", "Slc1a3", "Aqp4", "Aldh1l1", "Gfap"],
    "Oligo": ["Plp1", "Mog", "Mobp", "Mbp"],
    "OPC": ["Pdgfra", "Cspg4", "Sox10"],
    "Micro": ["P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Sall1", "Aif1"],
    "Endo": ["Pecam1", "Kdr", "Flt1", "Klf2", "Slco1a4"],
    "Peri": ["Pdgfrb", "Rgs5", "Kcnj8", "Abcc9"],
    "VLMC": ["Col1a1", "Col1a2", "Lum", "Dcn"],
    "SMC": ["Acta2", "Myh11", "Tagln"],
}

MAJOR_LABELS = [
    "Excit",
    "Inhib",
    "Astro",
    "Oligo",
    "OPC",
    "Micro",
    "Endo",
    "Peri",
    "VLMC",
    "SMC",
]


def plot_marker_genes(adata, save_dir=None):
    """Plot marker genes across clusters

    Args:
        adata: AnnData object with clustering results
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.
    """
    # Plot marker genes
    available_markers = []
    for cell_type, genes in MARKER_GENES.items():
        available = [g for g in genes if g in adata.var_names]
        available_markers.extend(available)

    # Dedupe while preserving order
    if available_markers:
        seen = set()
        available_markers = [
            g for g in available_markers if not (g in seen or seen.add(g))
        ]

    if available_markers:
        sc.pl.dotplot(
            adata,
            available_markers,
            groupby="leiden",
            standard_scale="var",
            figsize=(15, 8),
            show=False,
        )
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        if save_dir:
            plt.savefig(
                save_dir / "marker_genes_dotplot.png", dpi=300, bbox_inches="tight"
            )
            print(f"  Saved: {save_dir}/marker_genes_dotplot.png")
            plt.close()
        else:
            plt.show()


def annotate_cell_types(
    adata,
    save_dir=None,
    label_mode="cell",  # "cell" or "cluster"
    margin=0.05,
    cluster_agg="median",
    fallback_cluster_map=None,
):
    """Annotate cell types based on marker genes and clustering.

    Args:
        adata: AnnData object with clustering results.
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.
        label_mode: "cell" for per-cell labeling, "cluster" for cluster-level labeling.
        margin: Confidence margin between top and second-best scores.
        cluster_agg: Aggregation statistic for cluster-level scores ("median" or "mean").
        fallback_cluster_map: Optional mapping from leiden cluster id (str) to label; applied only
            to cells that remain unlabeled after scoring.

    Returns:
        AnnData object with cell type annotations
    """
    print("Annotating cell types...")

    # Plot marker genes first, saved to save_dir/marker_genes_dotplot.png
    plot_marker_genes(adata, save_dir=save_dir)

    # Score-based major cell type assignment
    if label_mode == "cluster":
        assign_major_celltypes_by_cluster_scores(adata, margin=margin, agg=cluster_agg)
    else:
        assign_major_celltypes_by_scores(adata, margin=margin)

    # Optional fallback: map any remaining unlabeled cells by a provided cluster map
    if fallback_cluster_map is not None:
        missing_mask = adata.obs.get("celltype").isna()
        if missing_mask.any():
            adata.obs.loc[missing_mask, "celltype"] = adata.obs.loc[
                missing_mask, "leiden"
            ].map(fallback_cluster_map)

    # Refine excitatory neurons into cortical layers where possible
    adata = refine_by_subtype_scores(
        adata,
        subtype_labels=["ExN_L2-4", "ExN_L5", "ExN_L6", "ExN_L6b"],
        eligible_celltypes=["ExN", "Excit", "Neuron"],
        margin=None,
    )

    # Refine inhibitory neuron subtypes where possible
    adata = refine_by_subtype_scores(
        adata,
        subtype_labels=["InN_SST", "InN_VIP", "InN_PVALB"],
        eligible_celltypes=["Inhib", "Neuron"],
        margin=None,
    )

    # Plot annotated cell types
    sc.pl.umap(
        adata,
        color="celltype",
        legend_loc="right margin",
        title="Cell type annotation",
        show=False,
    )

    if save_dir:
        plt.savefig(save_dir / "celltype_umap.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/celltype_umap.png")
        plt.close()
    else:
        plt.show()

    return adata


def _recluster_subset(
    adata,
    mask,
    subset_name,
    save_dir=None,
    n_top_genes=3000,
    n_pcs=30,
    n_neighbors=15,
    resolution=0.6,
    auto_resolution=True,
    resolution_grid=None,
    random_state=0,
    write_h5ad=False,
):
    """Re-cluster and UMAP a subset of cells and map labels back to parent AnnData.

    Args:
        adata: Parent AnnData (expects .raw to be set to full gene space).
        mask: Boolean array-like for cells to include in the subset.
        subset_name: Short name used in output keys/files (e.g., "excit", "inhib").
        save_dir: Optional Path for saving plots/files.
        n_top_genes: Number of HVGs to select within the subset.
        n_pcs: Number of PCs to compute/use for neighbors.
        n_neighbors: k for kNN graph.
        resolution: Leiden resolution if auto_resolution is False.
        auto_resolution: If True, sweep resolutions and pick robust choice.
        resolution_grid: Optional list of resolutions for sweeping.
        random_state: Random seed for UMAP.
        write_h5ad: If True, write subset AnnData to disk.
    """
    if save_dir is not None:
        save_dir = Path(save_dir)

    # Build subset from raw (full gene space) if available
    use_base = (
        adata.raw.to_adata() if getattr(adata, "raw", None) is not None else adata
    )
    sub = use_base[np.asarray(mask)].copy()
    if sub.n_obs == 0:
        return adata  # nothing to do

    # Subset-specific HVGs and processing
    sc.pp.highly_variable_genes(sub, flavor="seurat_v3", n_top_genes=int(n_top_genes))
    if "highly_variable" in sub.var:
        sub = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=min(int(n_pcs), max(2, sub.n_vars - 1)), svd_solver="arpack")
    sc.pp.neighbors(
        sub,
        n_neighbors=int(n_neighbors),
        n_pcs=min(int(n_pcs), sub.obsm["X_pca"].shape[1]),
    )

    # Choose resolution if requested
    chosen_res = float(resolution)
    if auto_resolution:
        chosen_res = choose_leiden_resolution(
            sub,
            resolution_grid=resolution_grid,
            min_cluster_size=20,
            save_dir="./plots/leiden_resolution_sweep/",
        )
        print(f"Chosen Leiden resolution: {chosen_res}")
    sc.tl.leiden(sub, resolution=float(chosen_res), key_added=f"leiden_{subset_name}")

    # UMAP for visualization
    sc.tl.umap(sub, random_state=int(random_state))

    # Plot and save
    if save_dir is not None:
        colors = [f"leiden_{subset_name}"]
        if "celltype" in sub.obs:
            colors.append("celltype")
        sc.pl.umap(
            sub,
            color=colors,
            legend_loc="right margin",
            show=False,
        )
        plt.savefig(save_dir / f"umap_{subset_name}.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/umap_{subset_name}.png")
        plt.close()

        if write_h5ad:
            out_path = save_dir / f"subset_{subset_name}.h5ad"
            sub.write(out_path)
            print(f"  Saved: {out_path}")

    # Map labels back to parent AnnData
    col = f"leiden_{subset_name}"
    if col not in adata.obs:
        adata.obs[col] = np.nan
    adata.obs.loc[sub.obs_names, col] = sub.obs[col].astype(str).values

    return adata


def recluster_excit_inhib(
    adata,
    save_dir=None,
    auto_resolution=True,
    resolution_grid=None,
    n_top_genes=3000,
    n_pcs=30,
    n_neighbors=15,
    random_state=0,
    write_h5ad=False,
):
    """Re-cluster excitatory and inhibitory subsets and save UMAPs.

    Adds `leiden_excit` and `leiden_inhib` to `adata.obs` for cells in each subset.
    """
    if "celltype" not in adata.obs:
        return adata

    ct = adata.obs["celltype"].astype(str)
    exc_mask = ct.str.startswith(("ExN_", "Excit"))
    inh_mask = ct.str.startswith(("InN_", "Inhib"))

    # Run subsets
    adata = _recluster_subset(
        adata,
        mask=exc_mask,
        subset_name="excit",
        save_dir=save_dir,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        resolution=0.6,
        auto_resolution=auto_resolution,
        resolution_grid=resolution_grid,
        random_state=random_state,
        write_h5ad=write_h5ad,
    )

    adata = _recluster_subset(
        adata,
        mask=inh_mask,
        subset_name="inhib",
        save_dir=save_dir,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        resolution=0.8,
        auto_resolution=auto_resolution,
        resolution_grid=resolution_grid,
        random_state=random_state,
        write_h5ad=write_h5ad,
    )

    return adata


def refine_by_subtype_scores(adata, subtype_labels, eligible_celltypes, margin=None):
    """Generic refinement by module scores within an eligible parent class.

    Computes per-cell module scores for each subtype label, selects the best
    label per cell, and assigns it only for cells whose current label is in
    eligible_celltypes. If margin is provided, only assigns when (best -
    second_best) >= margin; otherwise assigns unconditionally within eligible
    cells.

    Args:
        adata: AnnData with 'celltype' column present.
        subtype_labels: List of label keys that exist in get_marker_genes().
        eligible_celltypes: List of parent labels allowed to be refined.
        margin: Optional float; confidence margin threshold.

    Returns:
        AnnData with refined 'celltype' assignments.
    """
    if "celltype" not in adata.obs:
        return adata

    # Determine which gene space to use
    use_raw = getattr(adata, "raw", None) is not None and adata.raw is not None
    var_names = adata.raw.var_names if use_raw else adata.var_names

    # Compute module scores per subtype
    score_cols = []
    for label in subtype_labels:
        genes = [g for g in MARKER_GENES.get(label, []) if g in var_names]
        if not genes:
            continue
        score_name = f"score_{label}"
        sc.tl.score_genes(
            adata, gene_list=genes, score_name=score_name, use_raw=use_raw
        )
        score_cols.append(score_name)

    if not score_cols:
        return adata

    scores = adata.obs[score_cols].to_numpy()
    best_idx = np.argmax(scores, axis=1)
    score_labels = np.array([c.replace("score_", "") for c in score_cols])
    best_labels = score_labels[best_idx]

    # Optional margin gating
    if margin is not None and scores.shape[1] > 1:
        part = np.partition(scores, -2, axis=1)
        second_best = part[:, -2]
        best = scores[np.arange(scores.shape[0]), best_idx]
        confident = (best - second_best) >= float(margin)
    else:
        confident = np.ones(scores.shape[0], dtype=bool)

    eligible_mask = adata.obs["celltype"].isin(eligible_celltypes).to_numpy()
    update_mask = eligible_mask & confident
    if update_mask.any():
        adata.obs.loc[update_mask, "celltype"] = best_labels[update_mask]

    return adata


def assign_major_celltypes_by_scores(adata, margin=0.05):
    """Assign major cell types using module scores with a confidence margin.

    Creates per-panel scores (using raw if available), picks the top-scoring
    label when its score exceeds the next-best by `margin`. Otherwise leaves
    the celltype as NaN for fallback mapping.
    """
    use_raw = getattr(adata, "raw", None) is not None and adata.raw is not None
    var_names = adata.raw.var_names if use_raw else adata.var_names

    score_cols = []
    # Label by clusters rather than by individual cells
    for lbl in MAJOR_LABELS:
        genes = [g for g in MARKER_GENES.get(lbl, []) if g in var_names]
        if not genes:
            continue
        score_name = f"score_{lbl}"
        sc.tl.score_genes(
            adata, gene_list=genes, score_name=score_name, use_raw=use_raw
        )
        score_cols.append(score_name)

    if not score_cols:
        return

    scores = adata.obs[score_cols].to_numpy()
    top_idx = np.argmax(scores, axis=1)
    # Second best via partial sort
    part = np.partition(scores, -2, axis=1)
    second_best = part[:, -2]
    best = scores[np.arange(scores.shape[0]), top_idx]
    labels = np.array([c.replace("score_", "") for c in score_cols])
    winners = labels[top_idx]

    confident = best - second_best >= margin
    # Initialize or update celltype only where confident
    if "celltype" not in adata.obs:
        adata.obs["celltype"] = np.nan
    adata.obs.loc[confident, "celltype"] = winners[confident]


def assign_major_celltypes_by_cluster_scores(adata, margin=0.05, agg="median"):
    """Assign major cell types at the cluster level using module scores.

    This computes per-cell module scores, aggregates them per Leiden cluster,
    selects the top label per cluster, and assigns it to all cells in that
    cluster when the winning score exceeds the runner-up by `margin`.
    """
    use_raw = getattr(adata, "raw", None) is not None and adata.raw is not None
    var_names = adata.raw.var_names if use_raw else adata.var_names

    score_cols = []
    for lbl in MAJOR_LABELS:
        genes = [g for g in MARKER_GENES.get(lbl, []) if g in var_names]
        if not genes:
            continue
        score_name = f"score_{lbl}"
        sc.tl.score_genes(
            adata, gene_list=genes, score_name=score_name, use_raw=use_raw
        )
        score_cols.append(score_name)

    if not score_cols:
        return

    cluster_key = "leiden"
    if cluster_key not in adata.obs:
        return

    if "celltype" not in adata.obs:
        adata.obs["celltype"] = np.nan

    grouped = (
        adata.obs.groupby(cluster_key)[score_cols].median()
        if agg == "median"
        else adata.obs.groupby(cluster_key)[score_cols].mean()
    )

    grouped_vals = grouped.values
    top_idx = np.argmax(grouped_vals, axis=1)
    # second best via partial sort
    part = np.partition(grouped_vals, -2, axis=1)
    second_best = part[:, -2]
    best = grouped_vals[np.arange(grouped_vals.shape[0]), top_idx]
    labels = np.array([c.replace("score_", "") for c in score_cols])
    winners = labels[top_idx]

    confident = best - second_best >= margin

    for cluster_id, is_conf in zip(grouped.index.astype(str), confident):
        if not is_conf:
            continue
        label = winners[grouped.index.astype(str) == cluster_id][0]
        mask = adata.obs[cluster_key].astype(str) == cluster_id
        adata.obs.loc[mask, "celltype"] = label


def plot_cell_type_summary(adata, save_dir=None):
    """Plot summary of cell types across samples

    Args:
        adata: AnnData object with cell type annotations
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.
    """
    # Cell type counts
    celltype_counts = (
        adata.obs.groupby(["orig.ident", "celltype"]).size().unstack(fill_value=0)
    )

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    celltype_counts.plot(kind="bar", stacked=True, ax=ax)
    plt.title("Cell type distribution across samples")
    plt.xlabel("Sample")
    plt.ylabel("Number of cells")
    plt.xticks(rotation=45, ha="right")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if save_dir:
        fig.savefig(
            save_dir / "celltype_distribution.png", dpi=300, bbox_inches="tight"
        )
        print(f"  Saved: {save_dir}/celltype_distribution.png")
        plt.close(fig)
    else:
        plt.show()

    # Print summary table
    print("\nCell type summary:")
    print(adata.obs["celltype"].value_counts().sort_index())


def compute_top_markers_per_cluster(
    adata,
    groupby="leiden",
    method="wilcoxon",
    n_top=30,
    pval_adj_cutoff=None,
    save_dir=None,
    plot=False,
):
    """Compute top marker genes per cluster using differential expression.

    Args:
        adata: AnnData object with clustering results.
        groupby: Column in adata.obs to group by (default: "leiden").
        method: DE method passed to scanpy (e.g., "wilcoxon", "t-test").
        n_top: Number of top genes to rank per group.
        pval_adj_cutoff: Optional adjusted p-value cutoff to filter results.
        save_dir: Optional Path to save a CSV summary and optional plots.
        plot: If True, create a rank_genes_groups plot (saved if save_dir provided).

    Returns:
        Pandas DataFrame with ranked markers across all groups.
    """
    if groupby not in adata.obs:
        raise KeyError(f"Groupby key '{groupby}' not found in adata.obs")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=int(n_top),
        pts=True,
    )

    markers_df = sc.get.rank_genes_groups_df(adata, None)
    # Ensure consistent columns exist across scanpy versions
    if "pvals_adj" in markers_df.columns and pval_adj_cutoff is not None:
        markers_df = markers_df[markers_df["pvals_adj"] <= float(pval_adj_cutoff)]

    if save_dir is not None:
        out_csv = save_dir / "top_markers_by_cluster.csv"
        markers_df.to_csv(out_csv, index=False)
        print(f"  Saved: {out_csv}")

    if plot:
        sc.pl.rank_genes_groups(adata, n_genes=min(n_top, 20), sharey=False, show=False)
        if save_dir is not None:
            out_png = save_dir / "top_markers_ranked.png"
            plt.savefig(out_png, dpi=300, bbox_inches="tight")
            print(f"  Saved: {out_png}")
            plt.close()
        else:
            plt.show()

    return markers_df


def compare_top_markers_to_expected(
    adata,
    markers_df=None,
    groupby="leiden",
    top_n=10,
    panels=None,
    save_dir=None,
    plot=True,
):
    """Compare top DE genes per cluster with expected marker panels.

    Builds overlap metrics between each cluster's top-N DE genes and each
    expected marker gene panel.

    Args:
        adata: AnnData with DE results in .uns["rank_genes_groups"] or provide markers_df.
        markers_df: Optional DataFrame from sc.get.rank_genes_groups_df(adata, None).
        groupby: Cluster column used for DE (default: "leiden").
        top_n: Number of top genes per cluster to evaluate.
        panels: Optional dict mapping panel name -> list of genes. Defaults to MARKER_GENES.
        save_dir: Optional Path to write CSVs and heatmap.
        plot: If True, save a heatmap of precision (overlap/top_n).

    Returns:
        A tuple of (long_df, precision_matrix) where long_df is a list of dicts
        usable as rows for a DataFrame, and precision_matrix is a dict of
        {group: {panel: precision}}.
    """
    if panels is None:
        panels = MARKER_GENES

    # Make sets for faster overlap
    panel_to_genes = {k: set(v) for k, v in panels.items()}

    # Acquire ranked genes per group
    if markers_df is None:
        markers_df = sc.get.rank_genes_groups_df(adata, None)

    # Determine sort key preference
    sort_key = (
        "scores"
        if "scores" in markers_df.columns
        else ("logfoldchanges" if "logfoldchanges" in markers_df.columns else None)
    )
    if sort_key is None:
        raise KeyError("markers_df must contain 'scores' or 'logfoldchanges' column")

    # Build top-N gene sets per group
    var_names = set(adata.var_names)
    groups = sorted(markers_df["group"].unique())
    group_to_top = {}
    for g in groups:
        sub = markers_df[markers_df["group"] == g]
        sub = sub.sort_values(sort_key, ascending=False).head(int(top_n))
        genes = [nm for nm in sub["names"].tolist() if nm in var_names]
        group_to_top[g] = set(genes)

    # Compute overlaps
    long_rows = []
    precision_matrix = {g: {} for g in groups}
    for g in groups:
        top_set = group_to_top[g]
        for panel_name, panel_genes in panel_to_genes.items():
            overlap = len(top_set & panel_genes)
            precision = overlap / max(1, len(top_set))
            recall = overlap / max(1, len(panel_genes))
            denom = len(top_set | panel_genes)
            jaccard = overlap / max(1, denom)
            long_rows.append(
                {
                    "group": g,
                    "panel": panel_name,
                    "overlap": overlap,
                    "top_n": len(top_set),
                    "panel_size": len(panel_genes),
                    "precision": precision,
                    "recall": recall,
                    "jaccard": jaccard,
                }
            )
            precision_matrix[g][panel_name] = precision

    # Optional outputs
    if save_dir is not None:
        try:
            import pandas as _pd  # Local import to avoid hard dependency elsewhere

            long_df = _pd.DataFrame(long_rows)
            out_csv = save_dir / "expected_marker_overlap_long.csv"
            long_df.to_csv(out_csv, index=False)
            print(f"  Saved: {out_csv}")

            # Pivot to matrix for heatmap
            mat = long_df.pivot(index="group", columns="panel", values="precision")
            out_mat = save_dir / "expected_marker_overlap_matrix_precision.csv"
            mat.to_csv(out_mat)
            print(f"  Saved: {out_mat}")

            if plot:
                fig, ax = plt.subplots(
                    figsize=(
                        max(6, len(mat.columns) * 0.6),
                        max(4, len(mat.index) * 0.4),
                    )
                )
                im = ax.imshow(
                    mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=1
                )
                ax.set_xticks(range(len(mat.columns)))
                ax.set_xticklabels(mat.columns, rotation=45, ha="right")
                ax.set_yticks(range(len(mat.index)))
                ax.set_yticklabels(mat.index)
                ax.set_title("Precision: overlap of top-N vs expected markers")
                cbar = fig.colorbar(im, ax=ax)
                cbar.set_label("precision (overlap/top_n)")
                plt.tight_layout()
                out_png = save_dir / "expected_marker_overlap_heatmap.png"
                plt.savefig(out_png, dpi=300, bbox_inches="tight")
                print(f"  Saved: {out_png}")
                plt.close(fig)
        except Exception as e:
            print(f"Warning: could not write overlap CSVs/plot: {e}")

    return long_rows, precision_matrix
