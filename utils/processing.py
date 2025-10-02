#!/usr/bin/env python3
"""
Processing utilities for single-cell RNA-seq analysis
Handles normalization, scaling, PCA, UMAP, and clustering
"""

import scanpy as sc
import matplotlib.pyplot as plt


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


def run_pca_umap_clustering(adata, n_pcs=30, resolution=0.3):
    """Run PCA, UMAP and clustering
    
    Args:
        adata: AnnData object
        n_pcs: Number of principal components
        resolution: Clustering resolution
        
    Returns:
        AnnData object with embeddings and clusters
    """
    print("Running PCA...")
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)

    print("Running UMAP...")
    sc.tl.umap(adata)

    print("Clustering...")
    sc.tl.leiden(adata, resolution=resolution)

    return adata


def plot_embeddings(adata):
    """Plot UMAP embeddings
    
    Args:
        adata: AnnData object with UMAP coordinates
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
    plt.show()
