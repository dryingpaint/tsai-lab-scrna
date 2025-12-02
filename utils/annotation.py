#!/usr/bin/env python3
"""
Cell type annotation utilities for single-cell RNA-seq analysis
Handles marker gene analysis and cell type assignment
"""

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


def map_subtype_to_major(label):
    """Map subtype labels to their major cell type.
    
    Args:
        label: Cell type label (e.g., "ExN_L2-4", "InN_SST", "Astro")
    
    Returns:
        Major cell type label (e.g., "Excit", "Inhib", or the original label)
    """
    if label.startswith("ExN"):
        return "Excit"
    elif label.startswith("InN"):
        return "Inhib"
    else:
        return label


def plot_marker_genes(adata, marker_genes=MARKER_GENES, save_dir=None):
    """Plot marker genes across clusters

    Args:
        adata: AnnData object with clustering results
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.
    """
    # Plot marker genes
    available_markers = []
    for cell_type, genes in marker_genes.items():
        available = [g for g in genes if g in adata.var_names]
        available_markers.extend(available)

    # Dedupe while preserving order
    if available_markers:
        seen = set()
        available_markers = [
            g for g in available_markers if not (g in seen or seen.add(g))
        ]
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
    
    plt.show()


def create_cluster_aggregated_labels(adata, celltype_col='celltype', cluster_col='leiden', 
                                    purity_threshold=0.60):
    """Create cluster-level aggregated cell type labels with mixed cluster detection.
    
    For each cluster:
    - If dominant cell type is >purity_threshold: assigns that cell type
    - If dominant cell type is <=purity_threshold: labels as "Mixed" and stores top 2-3 cell types
    
    This identifies heterogeneous clusters that may need re-clustering or further analysis.
    
    Args:
        adata: AnnData object with cell type annotations
        celltype_col: Column name containing cell type labels
        cluster_col: Column name containing cluster labels
        purity_threshold: Threshold for cluster purity (default: 0.60 = 60%)
    
    Returns:
        None (modifies adata.obs in place, adds 'celltype_cluster' and 'celltype_cluster_top_types' columns)
    
    Side effects:
        - Adds 'celltype_cluster': cluster-level label ("celltype" or "Mixed")
        - Adds 'celltype_cluster_top_types': top 2-3 cell types for each cluster
        - Adds 'cluster_purity': proportion of dominant cell type in each cluster
        - Prints summary of mixed clusters
    """
    if celltype_col not in adata.obs or cluster_col not in adata.obs:
        return
    
    # Calculate cell type composition per cluster
    composition = pd.crosstab(
        adata.obs[cluster_col],
        adata.obs[celltype_col],
        normalize='index'
    )
    
    # Find dominant cell type and its proportion
    dominant = composition.idxmax(axis=1)
    dominant_prop = composition.max(axis=1)
    
    # For each cluster, get top 2-3 cell types
    top_types_dict = {}
    for cluster_id in composition.index:
        # Sort cell types by proportion (descending)
        sorted_types = composition.loc[cluster_id].sort_values(ascending=False)
        # Take top 3 types with >5% representation
        top_types = sorted_types[sorted_types > 0.05].head(3)
        # Format as string: "Type1 (30%), Type2 (25%), Type3 (20%)"
        top_str = ", ".join([f"{ct} ({prop*100:.1f}%)" for ct, prop in top_types.items()])
        top_types_dict[str(cluster_id)] = top_str
    
    # Determine cluster label: use dominant if pure, else "Mixed"
    cluster_labels = {}
    mixed_clusters = []
    
    for cluster_id in composition.index:
        if dominant_prop[cluster_id] > purity_threshold:
            cluster_labels[str(cluster_id)] = dominant[cluster_id]
        else:
            cluster_labels[str(cluster_id)] = "Mixed"
            mixed_clusters.append(cluster_id)
    
    # Assign to all cells
    adata.obs['celltype_cluster'] = adata.obs[cluster_col].astype(str).map(cluster_labels)
    adata.obs['celltype_cluster_top_types'] = adata.obs[cluster_col].astype(str).map(top_types_dict)
    adata.obs['cluster_purity'] = adata.obs[cluster_col].astype(str).map(
        lambda x: dominant_prop.get(str(x), 0.0)
    )
    
    # Print summary
    print(f"\n{'='*60}")
    print("CLUSTER PURITY ANALYSIS")
    print(f"{'='*60}")
    print(f"Purity threshold: {purity_threshold*100:.0f}%")
    print(f"Pure clusters: {len(cluster_labels) - len(mixed_clusters)}")
    print(f"Mixed clusters: {len(mixed_clusters)}")
    
    if mixed_clusters:
        print(f"\n{'='*60}")
        print("MIXED CLUSTERS (for further investigation)")
        print(f"{'='*60}")
        for cluster_id in mixed_clusters:
            print(f"\nCluster {cluster_id}: {dominant[cluster_id]} ({dominant_prop[cluster_id]*100:.1f}%)")
            print(f"  Top cell types: {top_types_dict[str(cluster_id)]}")
            n_cells = (adata.obs[cluster_col] == cluster_id).sum()
            print(f"  Total cells: {n_cells:,}")
        
        print(f"\n💡 Recommendation:")
        print(f"  - These {len(mixed_clusters)} clusters show heterogeneous cell type composition")
        print(f"  - Consider re-clustering these clusters at higher resolution")
        print(f"  - Or check if they represent transitional/doublet populations")
        print(f"  - Use 'celltype_cluster' == 'Mixed' to filter these cells")
    
    return mixed_clusters


def assign_major_celltypes_by_scores(adata, marker_genes=MARKER_GENES, margin=0.05, major_labels=MAJOR_LABELS):
    """Assign major cell types using module scores with a confidence margin.
    
    Uses a two-stage approach to avoid subtype markers overwhelming major type assignment:
    1. Stage 1: Score only major types (Excit, Inhib, Astro, etc.) - no subtypes
       - ALL cells assigned to best-scoring major type
       - Cells with (best - second_best) >= margin marked as "high" confidence
       - Remaining cells marked as "low" confidence
    2. Stage 2: Within Excit cells, score ExN subtypes; within Inhib cells, score InN subtypes
       - ALL cells assigned to best-scoring subtype
       - Confidence updated based on subtype score margin
    
    This prevents having 12 ExN subtypes compete against 1 Inhib label.

    Args:
        adata: AnnData object
        marker_genes: Dictionary of cell type markers
        margin: Confidence margin between top and second-best scores (default: 0.05)
        major_labels: List of major cell type labels to score in stage 1
        
    Creates columns:
        celltype: Major cell type (Excit, Inhib, Astro, etc.)
        celltype_detail: Detailed subtype (ExN_L5_IT, InN_SST, etc.)
        annotation_confidence: "high" or "low" based on score margin
    """
    use_raw = getattr(adata, "raw", None) is not None and adata.raw is not None
    var_names = adata.raw.var_names if use_raw else adata.var_names

    print("Stage 1: Assigning major cell types...")
    
    # STAGE 1: Score only major cell types
    score_cols = []
    for lbl in major_labels:
        genes = [g for g in marker_genes.get(lbl, []) if g in var_names]
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
    
    # Initialize celltype columns
    adata.obs["celltype"] = np.nan
    adata.obs["celltype_detail"] = np.nan
    adata.obs["annotation_confidence"] = "unassigned"
    
    # Assign major types to confident cells
    adata.obs.loc[confident, "celltype"] = winners[confident]
    adata.obs.loc[confident, "celltype_detail"] = winners[confident]
    adata.obs.loc[confident, "annotation_confidence"] = "high"
    
    # Fallback: assign remaining cells to best-scoring type (even if low confidence)
    low_confidence = ~confident
    if low_confidence.sum() > 0:
        adata.obs.loc[low_confidence, "celltype"] = winners[low_confidence]
        adata.obs.loc[low_confidence, "celltype_detail"] = winners[low_confidence]
        adata.obs.loc[low_confidence, "annotation_confidence"] = "low"
    
    print(f"  ✓ High confidence: {confident.sum():,} cells ({confident.sum()/len(confident)*100:.1f}%)")
    print(f"  ✓ Low confidence: {low_confidence.sum():,} cells ({low_confidence.sum()/len(confident)*100:.1f}%)")
    print(f"  ✓ Total assigned: {len(confident):,} cells (100%)")
    
    # STAGE 2: Refine Excit and Inhib with subtypes
    print("\nStage 2: Refining neuronal subtypes...")
    
    # Refine Excit → ExN subtypes
    excit_mask = adata.obs["celltype"] == "Excit"
    if excit_mask.sum() > 0:
        _refine_subtypes(adata, excit_mask, "ExN", marker_genes, margin, use_raw, var_names)
    
    # Refine Inhib → InN subtypes  
    inhib_mask = adata.obs["celltype"] == "Inhib"
    if inhib_mask.sum() > 0:
        _refine_subtypes(adata, inhib_mask, "InN", marker_genes, margin, use_raw, var_names)
    
    # Print final stats
    print("\n✓ Final cell type distribution:")
    major_counts = adata.obs["celltype"].value_counts()
    for ct, count in major_counts.items():
        print(f"  {ct}: {count:,}")


def _refine_subtypes(adata, mask, prefix, marker_genes, margin, use_raw, var_names):
    """Helper function to refine a major cell type into subtypes.
    
    Args:
        adata: AnnData object
        mask: Boolean mask of cells to refine
        prefix: Subtype prefix (e.g., "ExN" or "InN")
        marker_genes: Dictionary of marker genes
        margin: Confidence margin
        use_raw: Whether to use raw data
        var_names: Variable names to check
    """
    # Get all subtypes with this prefix
    subtypes = [k for k in marker_genes.keys() if k.startswith(prefix)]
    if not subtypes:
        return
    
    # Score subtypes only for cells with this major type
    score_cols = []
    for lbl in subtypes:
        genes = [g for g in marker_genes.get(lbl, []) if g in var_names]
        if not genes:
            continue
        score_name = f"score_{lbl}"
        # Only score if not already scored
        if score_name not in adata.obs.columns:
            sc.tl.score_genes(
                adata, gene_list=genes, score_name=score_name, use_raw=use_raw
            )
        score_cols.append(score_name)
    
    if not score_cols:
        return
    
    # Get scores only for masked cells
    subset_scores = adata.obs.loc[mask, score_cols].to_numpy()
    top_idx = np.argmax(subset_scores, axis=1)
    
    # Check confidence
    if subset_scores.shape[1] > 1:
        part = np.partition(subset_scores, -2, axis=1)
        second_best = part[:, -2]
        best = subset_scores[np.arange(subset_scores.shape[0]), top_idx]
        confident = (best - second_best) >= margin
    else:
        confident = np.ones(subset_scores.shape[0], dtype=bool)
    
    labels = np.array([c.replace("score_", "") for c in score_cols])
    winners = labels[top_idx]
    
    # Update celltype_detail for ALL cells (high and low confidence)
    mask_indices = np.where(mask)[0]
    
    # High confidence cells
    confident_indices = mask_indices[confident]
    adata.obs.loc[adata.obs.index[confident_indices], "celltype_detail"] = winners[confident]
    adata.obs.loc[adata.obs.index[confident_indices], "annotation_confidence"] = "high"
    
    # Low confidence cells - still assign best subtype
    low_conf = ~confident
    if low_conf.sum() > 0:
        low_conf_indices = mask_indices[low_conf]
        adata.obs.loc[adata.obs.index[low_conf_indices], "celltype_detail"] = winners[low_conf]
        # Keep their confidence as "low" from Stage 1, or set to "low" if was "high"
        current_conf = adata.obs.loc[adata.obs.index[low_conf_indices], "annotation_confidence"]
        adata.obs.loc[adata.obs.index[low_conf_indices], "annotation_confidence"] = "low"
    
    major_type = "Excit" if prefix == "ExN" else "Inhib"
    print(f"  {major_type}: {confident.sum():,} high confidence, {low_conf.sum():,} low confidence subtypes ({mask.sum():,} total)")


def assign_major_celltypes_by_cluster_scores(adata, marker_genes=MARKER_GENES, margin=0.05, agg="median", include_subtypes=True):
    """Assign major cell types at the cluster level using module scores.

    Computes scores for major cell types AND subtypes (e.g., ExN_L2-4, InN_SST).
    Aggregates scores per cluster, assigns the best-scoring label to the cluster,
    then maps subtypes back to their parent major type.
    This ensures cells with strong subtype markers are captured (e.g., ExN_L5 → Excit).

    Args:
        adata: AnnData object
        marker_genes: Dictionary of cell type markers to use for annotation.
        margin: Confidence margin between top and second-best scores
        agg: Aggregation method ('median' or 'mean')
        include_subtypes: If True, score subtypes in addition to major types (RECOMMENDED)
    """
    use_raw = getattr(adata, "raw", None) is not None and adata.raw is not None
    var_names = adata.raw.var_names if use_raw else adata.var_names

    # Select which labels to score
    labels_to_score = list(marker_genes.keys())

    score_cols = []
    for lbl in labels_to_score:
        genes = [g for g in marker_genes.get(lbl, []) if g in var_names]
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

    # Initialize celltype columns
    if "celltype" not in adata.obs:
        adata.obs["celltype"] = np.nan
    if include_subtypes and "celltype_detail" not in adata.obs:
        adata.obs["celltype_detail"] = np.nan

    # Aggregate scores per cluster
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

    # Track assignment statistics
    n_via_major = 0
    n_via_subtype = 0
    
    for cluster_id, is_conf in zip(grouped.index.astype(str), confident):
        if not is_conf:
            continue
        label = winners[grouped.index.astype(str) == cluster_id][0]
        mask = adata.obs[cluster_key].astype(str) == cluster_id
        
        if include_subtypes:
            # Store detailed subtype
            adata.obs.loc[mask, "celltype_detail"] = label
            
            # Map to major type
            major_label = map_subtype_to_major(label)
            adata.obs.loc[mask, "celltype"] = major_label
            
            # Track stats
            if major_label != label:
                n_via_subtype += mask.sum()
            else:
                n_via_major += mask.sum()
        else:
            adata.obs.loc[mask, "celltype"] = label
    
    # Report assignment statistics
    n_assigned = adata.obs["celltype"].notna().sum()
    if include_subtypes:
        print(f"✓ Assigned {confident.sum()} / {len(grouped)} clusters")
        print(f"  Total cells assigned: {n_assigned:,} / {adata.n_obs:,} ({n_assigned/adata.n_obs*100:.1f}%)")
        print(f"  Via major type markers: {n_via_major:,}")
        print(f"  Via subtype markers: {n_via_subtype:,}")
    else:
        print(f"✓ Assigned {confident.sum()} / {len(grouped)} clusters")
        print(f"  Total cells assigned: {n_assigned:,} / {adata.n_obs:,} ({n_assigned/adata.n_obs*100:.1f}%)")

