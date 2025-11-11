#!/usr/bin/env python3
"""
Processing utilities for single-cell RNA-seq analysis
Handles normalization, scaling, PCA, UMAP, and clustering
"""

import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score


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


def run_pca_umap_clustering(
    adata,
    n_pcs=15,
    resolution=0.6,
    save_dir=None,
    auto_resolution=False,
    resolution_grid=None,
    min_cluster_size=20,
):
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

    # Optionally perform a resolution sweep and choose an optimal resolution
    # chosen_resolution = resolution
    # if auto_resolution:
    #     print("Performing Leiden resolution sweep...")
    #     chosen_resolution = choose_leiden_resolution(
    #         adata,
    #         resolution_grid=resolution_grid,
    #         min_cluster_size=min_cluster_size,
    #         save_dir=save_dir,
    #     )
    #     adata.uns["leiden_optimal_resolution"] = float(chosen_resolution)
    #     print(f"Chosen Leiden resolution: {chosen_resolution}")

    print("Running UMAP...")
    sc.tl.umap(adata)

    print("Clustering...")
    # Determined 0.8 from the leiden resolution sweep
    sc.tl.leiden(adata, resolution=0.8)

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


def choose_leiden_resolution(
    adata,
    resolution_grid=None,
    min_cluster_size=20,
    save_dir=None,
):
    """Sweep Leiden resolutions and pick a robust choice.

    Strategy:
    - Compute Leiden for a grid of resolutions on the existing kNN graph
    - Evaluate silhouette on PCA space and fraction of cells in small clusters
    - Select the resolution with highest silhouette; among ties within 0.02 of max,
      prefer lower small-cluster fraction, then fewer clusters, then lower resolution

    Side effects:
    - Adds columns `leiden_{res}` to `adata.obs` for each tested resolution
    - Writes sweep metrics CSV and clustree-compatible labels CSV if `save_dir` set
    - Saves diagnostic plot of silhouette and number of clusters versus resolution

    Returns:
    - chosen resolution (float)
    """
    if resolution_grid is None:
        resolution_grid = np.round(np.arange(0.2, 2.05, 0.1), 2)

    # Use PCA embedding for silhouettes if present; otherwise run PCA minimally
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    X = adata.obsm.get("X_pca")

    metrics = []
    label_cols = []
    for res in resolution_grid:
        key = f"leiden_{res:.2f}"
        sc.tl.leiden(adata, resolution=float(res), key_added=key, directed=False)
        labels = adata.obs[key].astype(str)

        # Compute metrics only if >1 cluster and at least one cluster has >1 member
        n_clusters = labels.nunique()
        small_frac = 0.0
        sil = np.nan
        if n_clusters > 1:
            counts = labels.value_counts()
            small_frac = float(
                counts[counts < max(2, int(min_cluster_size))].sum() / len(labels)
            )
            try:
                sil = float(silhouette_score(X, labels))
            except Exception:
                sil = np.nan

        metrics.append(
            {
                "resolution": float(res),
                "n_clusters": int(n_clusters),
                "silhouette": sil,
                "small_cluster_fraction": small_frac,
            }
        )
        label_cols.append(key)

    metrics_df = pd.DataFrame(metrics)

    # Selection rule
    # 1) Take max silhouette; 2) among those within 0.02 of max, minimize small frac,
    # 3) then minimize n_clusters; 4) then choose lowest resolution
    valid = metrics_df.copy()
    max_sil = np.nanmax(valid["silhouette"].values)
    if np.isfinite(max_sil):
        near = valid[np.abs(valid["silhouette"] - max_sil) <= 0.02]
        near = near.sort_values(
            by=["small_cluster_fraction", "n_clusters", "resolution"],
            ascending=[True, True, True],
        )
        chosen = near.iloc[0]
    else:
        # Fallback: choose the lowest resolution with >1 cluster
        candidates = valid[valid["n_clusters"] > 1]
        chosen = (
            candidates.sort_values("resolution").iloc[0]
            if not candidates.empty
            else valid.sort_values("resolution").iloc[0]
        )

    chosen_res = (
        float(chosen["resolution"])
        if "resolution" in chosen
        else float(resolution_grid[0])
    )

    # Persist outputs
    if save_dir:
        try:
            os.makedirs(save_dir, exist_ok=True)
            metrics_path = save_dir / "leiden_resolution_sweep.csv"
            metrics_df.to_csv(metrics_path, index=False)

            # Clustree-compatible wide table
            clustree_df = adata.obs[label_cols].copy()
            clustree_df.insert(0, "cell", adata.obs_names)
            clustree_df.to_csv(save_dir / "clustree_leiden_labels.csv", index=False)

            # Diagnostic plot
            fig, ax1 = plt.subplots(figsize=(7, 4))
            ax2 = ax1.twinx()
            ax1.plot(
                metrics_df["resolution"],
                metrics_df["silhouette"],
                "-o",
                color="#1f77b4",
                label="Silhouette",
            )
            ax2.plot(
                metrics_df["resolution"],
                metrics_df["n_clusters"],
                "-s",
                color="#ff7f0e",
                label="#Clusters",
            )
            ax1.set_xlabel("Leiden resolution")
            ax1.set_ylabel("Silhouette (PCA)", color="#1f77b4")
            ax2.set_ylabel("# clusters", color="#ff7f0e")
            ax1.axvline(chosen_res, color="gray", linestyle="--", linewidth=1)
            fig.tight_layout()
            fig.savefig(
                save_dir / "leiden_sweep_diagnostics.png", dpi=300, bbox_inches="tight"
            )
            plt.close(fig)
            print(f"  Saved: {metrics_path}")
            print(f"  Saved: {save_dir}/clustree_leiden_labels.csv")
            print(f"  Saved: {save_dir}/leiden_sweep_diagnostics.png")
        except Exception:
            pass

    # Ensure `leiden` reflects the chosen resolution labels
    chosen_key = f"leiden_{chosen_res:.2f}"
    if chosen_key in adata.obs:
        adata.obs["leiden"] = adata.obs[chosen_key].astype(str)

    return chosen_res
