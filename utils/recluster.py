import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from utils.processing import choose_leiden_resolution
from utils.annotation import MARKER_GENES


def _assign_subtypes_to_subset(
    adata, 
    sub, 
    subset_name, 
    parent_type, 
    subtype_labels,
    marker_genes=MARKER_GENES,
    use_cluster_assignment=True,
    confidence_margin=0.05,
):
    """Assign subtypes to a subset based on marker gene scores.
    
    Args:
        adata: Parent AnnData with full gene space in .raw
        sub: Subset AnnData (may be filtered to HVGs)
        subset_name: Name of subset (e.g., "excit", "inhib")
        parent_type: Parent cell type (e.g., "Excit", "Inhib")
        subtype_labels: List of subtype labels to assign
        use_cluster_assignment: If True, assign subtypes at cluster level (respects UMAP structure)
        confidence_margin: Minimum score difference between top and second-best subtype
        
    Returns:
        Updated sub with subtypes assigned
    """
    if not subtype_labels or not parent_type:
        return sub
    
    cluster_key = f"leiden_{subset_name}"
    if use_cluster_assignment and cluster_key not in sub.obs:
        print(f"  Warning: {cluster_key} not found, falling back to per-cell assignment")
        use_cluster_assignment = False
    
    # Get full gene space for marker scoring
    use_base = (
        adata.raw.to_adata() if getattr(adata, "raw", None) is not None else adata
    )
    sub_full = use_base[sub.obs_names].copy()
    
    # Ensure celltype exists in sub_full (copy from parent)
    if "celltype" in adata.obs:
        sub_full.obs["celltype"] = adata.obs.loc[sub_full.obs_names, "celltype"].values
    
    # Copy cluster labels to sub_full
    if cluster_key in sub.obs:
        sub_full.obs[cluster_key] = sub.obs[cluster_key].values
    
    # Compute subtype scores on full gene space
    var_names = sub_full.var_names
    score_cols = []
    for label in subtype_labels:
        genes = [g for g in marker_genes.get(label, []) if g in var_names]
        if not genes:
            continue
        score_name = f"score_{label}"
        sc.tl.score_genes(sub_full, gene_list=genes, score_name=score_name, use_raw=False)
        score_cols.append(score_name)
    
    if not score_cols:
        return sub
    
    score_labels = np.array([c.replace("score_", "") for c in score_cols])
    
    # Initialize assignments
    if "celltype" not in sub_full.obs:
        sub_full.obs["celltype"] = parent_type
    if "celltype" not in sub.obs:
        sub.obs["celltype"] = parent_type
    
    parent_mask = sub_full.obs["celltype"].astype(str) == parent_type
    n_parent = parent_mask.sum()
    
    if use_cluster_assignment and n_parent > 0:
        # CLUSTER-BASED ASSIGNMENT: Assign subtypes at cluster level
        print(f"  Using cluster-based assignment (respects UMAP structure)")
        
        # Get cells that need assignment
        cells_to_assign = sub_full.obs.index[parent_mask]
        
        # Aggregate scores per cluster (median is more robust than mean)
        cluster_scores = sub_full.obs.loc[cells_to_assign].groupby(cluster_key)[score_cols].median()
        
        # Find best subtype per cluster with confidence check
        cluster_assignments = {}
        for cluster_id in cluster_scores.index:
            scores = cluster_scores.loc[cluster_id].values
            top_idx = np.argmax(scores)
            top_score = scores[top_idx]
            
            # Get second best score
            if len(scores) > 1:
                second_best = np.partition(scores, -2)[-2]
                confidence = top_score - second_best
            else:
                confidence = top_score
            
            # Only assign if confident
            if confidence >= confidence_margin:
                cluster_assignments[cluster_id] = score_labels[top_idx]
        
        # Apply cluster assignments
        n_assigned = 0
        n_uncertain = 0
        for cluster_id, subtype in cluster_assignments.items():
            mask = (sub_full.obs[cluster_key] == cluster_id) & parent_mask
            sub_full.obs.loc[mask, "celltype"] = subtype
            sub.obs.loc[mask, "celltype"] = subtype
            n_assigned += mask.sum()
        
        n_uncertain = n_parent - n_assigned
        
        print(f"  Assigned subtypes to {n_assigned:,} cells across {len(cluster_assignments)} confident clusters")
        if n_uncertain > 0:
            print(f"  {n_uncertain:,} cells in uncertain clusters remain as '{parent_type}'")
        
    else:
        # PER-CELL ASSIGNMENT: Original behavior with confidence threshold
        print(f"  Using per-cell assignment")
        
        scores = sub_full.obs[score_cols].to_numpy()
        
        # Calculate confidence for each cell
        top_idx = np.argmax(scores, axis=1)
        top_scores = scores[np.arange(len(scores)), top_idx]
        
        # Get second best scores
        if scores.shape[1] > 1:
            scores_sorted = np.partition(scores, -2, axis=1)
            second_best = scores_sorted[:, -2]
            confidence = top_scores - second_best
        else:
            confidence = top_scores
        
        # Assign only confident cells
        confident_mask = confidence >= confidence_margin
        best_labels = score_labels[top_idx]
        
        # Apply assignments
        assign_mask = parent_mask & confident_mask
        sub_full.obs.loc[assign_mask, "celltype"] = best_labels[parent_mask][confident_mask]
        sub.obs.loc[assign_mask, "celltype"] = best_labels[parent_mask][confident_mask]
        
        n_assigned = assign_mask.sum()
        n_uncertain = parent_mask.sum() - n_assigned
        
        print(f"  Assigned subtypes to {n_assigned:,} confident cells")
        if n_uncertain > 0:
            print(f"  {n_uncertain:,} uncertain cells remain as '{parent_type}'")
    
    # Report final distribution
    if n_assigned > 0:
        subtype_counts = sub.obs.loc[parent_mask, "celltype"].value_counts()
        print(f"  Subtype distribution:")
        for subtype, count in subtype_counts.items():
            print(f"    {subtype}: {count:,}")
    
    return sub


def _map_subset_labels_to_parent(adata, sub, col_name):
    """Map a column from subset back to parent AnnData.
    
    Ensures proper string conversion for h5py compatibility.
    
    Args:
        adata: Parent AnnData
        sub: Subset AnnData
        col_name: Name of column to map
    """
    if col_name not in sub.obs:
        return
    
    # Get labels and convert to Python strings (required for h5py)
    labels = sub.obs[col_name].astype(str)
    labels = labels.replace('nan', '')
    
    # Initialize column in parent if needed
    if col_name not in adata.obs:
        adata.obs[col_name] = ''
    
    # Map labels back to parent
    label_series = pd.Series(
        [str(x) if pd.notna(x) and str(x) != 'nan' else '' for x in labels],
        index=sub.obs_names,
        dtype='object'
    )
    adata.obs.loc[sub.obs_names, col_name] = label_series
    
    # Ensure entire column is object dtype with Python strings (h5py requirement)
    adata.obs[col_name] = adata.obs[col_name].astype(str).astype('object')


def _recluster_subset(
    adata,
    mask,
    subset_name,
    save_dir=None,
    n_top_genes=3000,
    n_pcs=30,
    n_neighbors=15,
    resolution=0.2,
    auto_resolution=False,
    resolution_grid=None,
    random_state=0,
    write_h5ad=False,
    use_cluster_assignment=True,
    confidence_margin=0.05,
    marker_genes=MARKER_GENES,
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
        use_cluster_assignment: If True, assign subtypes at cluster level (better UMAP alignment).
        confidence_margin: Minimum score difference to confidently assign a subtype.
        marker_genes: Dictionary of cell type markers to use for annotation.
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
    
    # Copy celltype from parent (raw data doesn't include all obs columns)
    if "celltype" in adata.obs:
        sub.obs["celltype"] = adata.obs.loc[sub.obs_names, "celltype"].values

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

    # Assign subtypes based on subset type
    # Automatically detect all subtypes from marker_genes
    if subset_name == "excit":
        parent_type = "Excit"
        # Get all ExN_* subtypes from marker_genes
        subtype_labels = [k for k in marker_genes.keys() if k.startswith("ExN")]
        print(f"\n  Detected {len(subtype_labels)} excitatory subtypes: {', '.join(subtype_labels)}")
    elif subset_name == "inhib":
        parent_type = "Inhib"
        # Get all InN_* subtypes from marker_genes
        subtype_labels = [k for k in marker_genes.keys() if k.startswith("InN")]
        print(f"\n  Detected {len(subtype_labels)} inhibitory subtypes: {', '.join(subtype_labels)}")
    else:
        parent_type, subtype_labels = None, []
    
    sub = _assign_subtypes_to_subset(
        adata, 
        sub, 
        subset_name, 
        parent_type, 
        subtype_labels,
        marker_genes=marker_genes,
        use_cluster_assignment=use_cluster_assignment,
        confidence_margin=confidence_margin,
    )

    # Plot and save
    if save_dir is not None:
        # Plot leiden clusters - display in notebook
        print(f"\n  UMAP for {subset_name} (Leiden clusters):")
        sc.pl.umap(
            sub,
            color=f"leiden_{subset_name}",
            legend_loc="right margin",
            title=f"Re-clustered {subset_name} - Leiden clusters",
            show=True,  # Display in notebook
            save=None,  # Don't save via scanpy (we'll save manually)
        )
        # Get current figure and save it
        fig = plt.gcf()
        fig.savefig(save_dir / f"umap_{subset_name}.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/umap_{subset_name}.png")
        plt.close(fig)
        
        # Plot celltype if available - display in notebook
        if "celltype" in sub.obs:
            print(f"\n  UMAP for {subset_name} (Cell types):")
            sc.pl.umap(
                sub,
                color="celltype",
                legend_loc="right margin",
                title=f"Re-clustered {subset_name} - Cell types",
                show=True,  # Display in notebook
                save=None,  # Don't save via scanpy (we'll save manually)
            )
            # Get current figure and save it
            fig = plt.gcf()
            fig.savefig(save_dir / f"umap_{subset_name}_celltype.png", dpi=300, bbox_inches="tight")
            print(f"  Saved: {save_dir}/umap_{subset_name}_celltype.png")
            plt.close(fig)

        if write_h5ad:
            out_path = save_dir / f"subset_{subset_name}.h5ad"
            sub.write(out_path)
            print(f"  Saved: {out_path}")

    # Map labels back to parent AnnData
    _map_subset_labels_to_parent(adata, sub, f"leiden_{subset_name}")
    if "celltype" in sub.obs:
        _map_subset_labels_to_parent(adata, sub, "celltype")

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
    use_cluster_assignment=True,
    confidence_margin=0.05,
    marker_genes=MARKER_GENES,
):
    """Re-cluster excitatory and inhibitory subsets and save UMAPs.

    Adds `leiden_excit` and `leiden_inhib` to `adata.obs` for cells in each subset.
    
    Args:
        adata: AnnData object with 'celltype' column
        save_dir: Directory to save plots and files
        auto_resolution: If True, automatically select optimal Leiden resolution
        resolution_grid: List of resolutions to test if auto_resolution=True
        n_top_genes: Number of HVGs to compute within each subset
        n_pcs: Number of PCs for dimensionality reduction
        n_neighbors: Number of neighbors for kNN graph
        random_state: Random seed for reproducibility
        write_h5ad: If True, save subset h5ad files
        use_cluster_assignment: If True, assign subtypes at cluster level (RECOMMENDED).
            This ensures cell type labels align with UMAP structure.
        confidence_margin: Minimum score difference to assign a subtype (default 0.05).
            Higher values = more conservative, fewer cells assigned.
        marker_genes: Dictionary of cell type markers to use for annotation.
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
        resolution=0.2,
        auto_resolution=auto_resolution,
        resolution_grid=resolution_grid,
        random_state=random_state,
        write_h5ad=write_h5ad,
        use_cluster_assignment=use_cluster_assignment,
        confidence_margin=confidence_margin,
        marker_genes=marker_genes,
    )

    adata = _recluster_subset(
        adata,
        mask=inh_mask,
        subset_name="inhib",
        save_dir=save_dir,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        resolution=0.2,
        auto_resolution=auto_resolution,
        resolution_grid=resolution_grid,
        random_state=random_state,
        write_h5ad=write_h5ad,
        use_cluster_assignment=use_cluster_assignment,
        confidence_margin=confidence_margin,
        marker_genes=marker_genes,
    )

    return adata