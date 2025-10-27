#!/usr/bin/env python3
"""
Alternative QC violin plot implementation using seaborn for true violin plots
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def create_true_violin_plots(adata, save_dir=None):
    """Create proper violin plots for QC metrics

    Args:
        adata: AnnData object with QC metrics
        save_dir: Directory to save plots (optional)
    """
    print("Creating true violin plots...")

    # Extract QC data to DataFrame
    qc_data = pd.DataFrame(
        {
            "n_genes_by_counts": adata.obs["n_genes_by_counts"],
            "total_counts": adata.obs["total_counts"],
            "percent_mt": adata.obs["percent_mt"],
            "percent_ribo": adata.obs["percent_ribo"],
        }
    )

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    # Define metrics and their properties
    metrics = [
        ("n_genes_by_counts", "Genes per cell", axes[0]),
        ("total_counts", "Total counts per cell", axes[1]),
        ("percent_mt", "Mitochondrial %", axes[2]),
        ("percent_ribo", "Ribosomal %", axes[3]),
    ]

    # Create violin plot for each metric
    for metric, title, ax in metrics:
        # Create the violin plot
        sns.violinplot(data=qc_data, y=metric, ax=ax, color="skyblue", inner="box")

        # Customize the plot
        ax.set_ylabel(title)
        ax.set_xlabel("")
        ax.set_title(title)

        # Add horizontal lines for common thresholds
        if metric == "percent_mt":
            ax.axhline(
                y=10, color="red", linestyle="--", alpha=0.5, label="10% threshold"
            )
        elif metric == "n_genes_by_counts":
            ax.axhline(
                y=200, color="red", linestyle="--", alpha=0.5, label="Min threshold"
            )
            ax.axhline(
                y=8000, color="red", linestyle="--", alpha=0.5, label="Max threshold"
            )

    plt.tight_layout()

    if save_dir:
        fig.savefig(save_dir / "true_violin_plots.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/true_violin_plots.png")
        plt.close(fig)
    else:
        plt.show()

    return fig


def create_violin_with_points(adata, save_dir=None):
    """Create violin plots with overlaid strip plots

    Args:
        adata: AnnData object with QC metrics
        save_dir: Directory to save plots (optional)
    """
    print("Creating violin plots with points...")

    # Sample data if too many cells (for visibility)
    if adata.n_obs > 10000:
        import numpy as np

        sample_idx = np.random.choice(adata.n_obs, 5000, replace=False)
        adata_sample = adata[sample_idx]
    else:
        adata_sample = adata

    # Extract QC data
    qc_data = pd.DataFrame(
        {
            "n_genes_by_counts": adata_sample.obs["n_genes_by_counts"],
            "total_counts": adata_sample.obs["total_counts"],
            "percent_mt": adata_sample.obs["percent_mt"],
            "percent_ribo": adata_sample.obs["percent_ribo"],
        }
    )

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    metrics = [
        ("n_genes_by_counts", "Genes per cell", axes[0]),
        ("total_counts", "Total counts per cell", axes[1]),
        ("percent_mt", "Mitochondrial %", axes[2]),
        ("percent_ribo", "Ribosomal %", axes[3]),
    ]

    for metric, title, ax in metrics:
        # Violin plot
        sns.violinplot(data=qc_data, y=metric, ax=ax, color="lightblue", inner=None)

        # Overlay strip plot with small points
        sns.stripplot(
            data=qc_data, y=metric, ax=ax, color="black", alpha=0.3, size=1, jitter=True
        )

        ax.set_ylabel(title)
        ax.set_xlabel("")
        ax.set_title(title)

    plt.tight_layout()

    if save_dir:
        fig.savefig(save_dir / "violin_with_points.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/violin_with_points.png")
        plt.close(fig)
    else:
        plt.show()

    return fig


def compare_distributions(adata, groupby="orig.ident", save_dir=None):
    """Create violin plots comparing distributions across groups

    Args:
        adata: AnnData object
        groupby: Column to group by (e.g., 'orig.ident', 'condition')
        save_dir: Directory to save plots
    """
    # Prepare data
    plot_data = adata.obs[
        [groupby, "n_genes_by_counts", "total_counts", "percent_mt", "percent_ribo"]
    ].copy()

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()

    metrics = [
        ("n_genes_by_counts", "Genes per cell"),
        ("total_counts", "Total counts per cell"),
        ("percent_mt", "Mitochondrial %"),
        ("percent_ribo", "Ribosomal %"),
    ]

    for idx, (metric, title) in enumerate(metrics):
        sns.violinplot(data=plot_data, x=groupby, y=metric, ax=axes[idx])
        axes[idx].set_title(title)
        axes[idx].set_xlabel("")
        axes[idx].tick_params(axis="x", rotation=45)

    plt.tight_layout()

    if save_dir:
        fig.savefig(save_dir / f"violin_by_{groupby}.png", dpi=300, bbox_inches="tight")
        print(f"  Saved: {save_dir}/violin_by_{groupby}.png")
        plt.close(fig)
    else:
        plt.show()

    return fig



