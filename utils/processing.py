#!/usr/bin/env python3
"""
Processing utilities for single-cell RNA-seq analysis
Handles normalization, scaling, PCA, UMAP, and clustering
"""

import scanpy as sc
import matplotlib.pyplot as plt
import os


def normalize_and_scale(adata):
    """Normalize and scale data

    Args:
        adata: AnnData object

    Returns:
        Processed AnnData object
    """
    print("Normalizing and scaling data...")

    # Save raw counts
    adata.raw = adata

    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log transform
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # Keep only highly variable genes for downstream analysis
    adata.raw = adata  # Save full data
    adata = adata[:, adata.var.highly_variable]

    # Scale data
    sc.pp.scale(adata, max_value=10)

    return adata


def run_pca_umap_clustering(adata, n_pcs=15, resolution=0.6, save_dir=None):
    """Run PCA, UMAP and clustering

    Args:
        adata: AnnData object
        n_pcs: Number of principal components
        resolution: Clustering resolution
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.

    Returns:
        AnnData object with embeddings and clusters
    """
    print("Running PCA...")
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    # Plot elbow plot
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        # Use scanpy's figdir mechanism to avoid path concatenation issues
        original_figdir = sc.settings.figdir
        sc.settings.figdir = str(save_dir)
        sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False, save="_elbow.png")
        # Rename to a stable filename if present
        try:
            import os as _os

            src = save_dir / "pca_variance_ratio_elbow.png"
            dst = save_dir / "pca_elbow_plot.png"
            if src.exists():
                _os.rename(src, dst)
            print(f"  Saved: {save_dir}/pca_elbow_plot.png")
        finally:
            sc.settings.figdir = original_figdir

    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)

    print("Running UMAP...")
    sc.tl.umap(adata)

    print("Clustering...")
    sc.tl.leiden(adata, resolution=resolution)

    return adata


def plot_embeddings(adata, save_dir=None):
    """Plot UMAP embeddings

    Args:
        adata: AnnData object with UMAP coordinates
        save_dir: Directory to save plots (optional). If provided, plots are saved without display.
    """
    print("Plotting embeddings...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    sc.pl.umap(
        adata,
        color="leiden",
        legend_loc="on data",
        title="Leiden clustering",
        ax=axes[0, 0],
        show=False,
    )
    sc.pl.umap(adata, color="Genotype", title="Genotype", ax=axes[0, 1], show=False)
    sc.pl.umap(adata, color="Sex", title="Sex", ax=axes[1, 0], show=False)
    sc.pl.umap(
        adata, color="Stimulation", title="Stimulation", ax=axes[1, 1], show=False
    )

    plt.tight_layout()

    if save_dir:
        fig.savefig(save_dir / "umap_embeddings.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/umap_embeddings.png")
        plt.close(fig)
    else:
        plt.show()
