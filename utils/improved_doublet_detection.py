#!/usr/bin/env python3
"""
Improved doublet detection utilities for single-cell RNA-seq analysis
"""

import numpy as np
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt


def detect_doublets_improved(
    adata,
    sample_col="orig.ident",
    expected_doublet_rate=0.06,
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=30,
    manual_threshold=None,
    plot_histograms=True,
    save_dir=None,
):
    """
    Improved doublet detection using Scrublet with better parameter handling

    Args:
        adata: AnnData object (should be after basic QC filtering)
        sample_col: Column name for sample identification
        expected_doublet_rate: Expected doublet rate (default 0.06 for 10x)
        manual_threshold: If set, use this threshold instead of automatic
        plot_histograms: Whether to plot doublet score histograms
        save_dir: Directory to save diagnostic plots

    Returns:
        AnnData object with doublet predictions and scores
    """
    print("Running improved doublet detection...")

    # Store results
    all_scores = np.zeros(adata.n_obs)
    all_predictions = np.zeros(adata.n_obs, dtype=bool)

    # Process each sample separately
    samples = adata.obs[sample_col].unique()

    if plot_histograms and save_dir:
        fig, axes = plt.subplots(4, 4, figsize=(16, 12))
        axes = axes.flatten()

    for idx, sample in enumerate(samples):
        print(f"\nProcessing sample: {sample}")

        # Get sample mask
        mask = adata.obs[sample_col] == sample
        sample_indices = np.where(mask)[0]

        # Extract sample data
        adata_sample = adata[mask].copy()

        # Skip if too few cells
        if adata_sample.n_obs < 100:
            print(f"  Skipping - only {adata_sample.n_obs} cells")
            continue

        # Initialize Scrublet
        scrub = scr.Scrublet(
            adata_sample.X, expected_doublet_rate=expected_doublet_rate
        )

        # Run doublet detection
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=min_counts,
            min_cells=min_cells,
            min_gene_variability_pctl=min_gene_variability_pctl,
            n_prin_comps=n_prin_comps,
            verbose=False,
        )

        # Get the automatic threshold
        auto_threshold = scrub.call_doublets(threshold=None)[1]

        # Use manual threshold if specified
        if manual_threshold is not None:
            threshold = manual_threshold
            predicted_doublets = doublet_scores > threshold
        else:
            threshold = auto_threshold

        # Apply a minimum threshold to avoid missing doublets
        if threshold > 0.4:
            print(f"  Warning: High auto threshold {threshold:.2f}, capping at 0.4")
            threshold = 0.4
            predicted_doublets = doublet_scores > threshold

        # Store results
        all_scores[sample_indices] = doublet_scores
        all_predictions[sample_indices] = predicted_doublets

        # Calculate statistics
        n_doublets = predicted_doublets.sum()
        pct_doublets = n_doublets / len(doublet_scores) * 100

        print(f"  Cells: {len(doublet_scores)}")
        print(f"  Threshold: {threshold:.3f}")
        print(f"  Doublets: {n_doublets} ({pct_doublets:.1f}%)")
        print(f"  Score range: {doublet_scores.min():.3f} - {doublet_scores.max():.3f}")

        # Plot histogram if requested
        if plot_histograms and save_dir and idx < 16:
            ax = axes[idx]
            ax.hist(doublet_scores, bins=50, alpha=0.7, edgecolor="black")
            ax.axvline(
                threshold,
                color="red",
                linestyle="--",
                label=f"Threshold: {threshold:.2f}",
            )
            ax.set_title(f"{sample}\n{n_doublets} doublets ({pct_doublets:.1f}%)")
            ax.set_xlabel("Doublet Score")
            ax.set_ylabel("Frequency")
            ax.legend()

    # Save histogram plot
    if plot_histograms and save_dir:
        plt.tight_layout()
        plt.savefig(
            save_dir / "doublet_score_histograms.png", dpi=300, bbox_inches="tight"
        )
        print(
            f"\nSaved doublet score histograms to {save_dir}/doublet_score_histograms.png"
        )
        plt.close()

    # Add results to adata
    adata.obs["doublet_score"] = all_scores
    adata.obs["predicted_doublet"] = all_predictions

    # Overall summary
    total_doublets = all_predictions.sum()
    total_cells = len(all_predictions)
    overall_rate = total_doublets / total_cells * 100

    print("\nOverall doublet detection summary:")
    print(f"  Total cells: {total_cells}")
    print(f"  Total doublets: {total_doublets}")
    print(f"  Overall rate: {overall_rate:.1f}%")

    # Per-sample summary
    summary = adata.obs.groupby(sample_col).agg(
        {
            "predicted_doublet": [
                "count",
                "sum",
                lambda x: (x.sum() / len(x) * 100).round(1),
            ]
        }
    )
    summary.columns = ["n_cells", "n_doublets", "pct_doublets"]
    print("\nPer-sample summary:")
    print(summary)

    return adata


def plot_doublet_scores_umap(adata, save_dir=None):
    """
    Plot doublet scores on UMAP to visualize their distribution

    Args:
        adata: AnnData object with UMAP and doublet scores
        save_dir: Directory to save plot
    """
    if "X_umap" not in adata.obsm:
        print("No UMAP found, skipping doublet score visualization")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot doublet scores
    sc.pl.umap(
        adata,
        color="doublet_score",
        ax=ax1,
        show=False,
        title="Doublet Scores",
        cmap="Reds",
    )

    # Plot predicted doublets
    sc.pl.umap(
        adata,
        color="predicted_doublet",
        ax=ax2,
        show=False,
        title="Predicted Doublets",
        palette=["lightgray", "red"],
    )

    plt.tight_layout()

    if save_dir:
        plt.savefig(save_dir / "doublet_umap.png", dpi=300, bbox_inches="tight")
        print(f"Saved doublet UMAP to {save_dir}/doublet_umap.png")
        plt.close()
    else:
        plt.show()


def remove_doublets_iterative(adata, max_iterations=3, threshold_increment=0.05):
    """
    Iteratively remove doublets, re-running detection after each round

    Args:
        adata: AnnData object
        max_iterations: Maximum number of iterations
        threshold_increment: Decrease threshold by this amount each iteration

    Returns:
        AnnData object with doublets removed
    """
    print("Running iterative doublet removal...")

    for i in range(max_iterations):
        print(f"\nIteration {i+1}:")

        # Detect doublets
        adata = detect_doublets_improved(
            adata,
            manual_threshold=0.4 - (i * threshold_increment),
            plot_histograms=False,
        )

        # Check if any doublets found
        n_doublets = adata.obs["predicted_doublet"].sum()
        if n_doublets == 0:
            print("No doublets found, stopping iteration")
            break

        # Remove doublets
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        print(f"Removed {n_doublets} doublets, {adata.n_obs} cells remaining")

    return adata
