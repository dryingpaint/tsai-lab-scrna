#!/usr/bin/env python3
"""
Cell type annotation utilities for single-cell RNA-seq analysis
Handles marker gene analysis and cell type assignment
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


# Module-level constants: single sources of truth
MARKER_GENES = {
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
    for cell_type, genes in s.items():
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
