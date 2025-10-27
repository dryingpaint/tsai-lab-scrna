#!/usr/bin/env python3
"""
Single-cell RNA-seq data preprocessing, QC, and cell type annotation
Refactored modular version

This script performs:
1. CellBender data loading and integration
2. Quality control and doublet detection
3. Normalization and dimensionality reduction
4. Clustering and cell type annotation

uv run python cellbender_qc_annotation.py
"""

import warnings
import argparse
import matplotlib
import scanpy as sc
from pathlib import Path

# Import our custom modules
from utils.data_loader import load_and_merge_cellbender_data, add_metadata
from utils.qc_utils import (
    calculate_qc_metrics,
    plot_qc_metrics,
    filter_cells_and_genes,
)
from utils.processing import (
    normalize_and_scale,
    run_pca_umap_clustering,
    plot_embeddings,
)
from utils.annotation import annotate_cell_types, plot_cell_type_summary
from utils.qc_filters import CELL_FILTERS, DOUBLET_PARAMS, get_filter_summary

# Configure scanpy
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor="white")

# Suppress warnings
warnings.filterwarnings("ignore")


def main(label_mode="cell", plots_dir_path="plots"):
    """Main analysis pipeline

    Args:
        label_mode: "cell" for per-cell labeling or "cluster" for cluster-level labeling.
        plots_dir_path: Directory where plots will be saved.
    """
    print("Starting single-cell analysis pipeline...")

    # Create output directory for plots
    plots_dir = Path(plots_dir_path)
    plots_dir.mkdir(exist_ok=True)
    print(f"Plots will be saved to: {plots_dir.absolute()}")

    # Set matplotlib backend to non-interactive for save-only mode
    matplotlib.use("Agg")
    print("Running in save-only mode - plots will not be displayed")

    # Print filter settings
    print("\n" + get_filter_summary() + "\n")

    # Define parameters
    base_path = "/Users/melissadu/projects/tsai-lab-scrna/data/"
    sample_names = [f"D25-{i}" for i in range(2675, 2691)]
    custom_name = "_processed_feature_bc_matrix_filtered.h5"

    # Step 1: Load and merge data
    adata = load_and_merge_cellbender_data(base_path, sample_names, custom_name)

    # Step 2: Add metadata
    adata = add_metadata(adata, sample_names)

    # Step 3: Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    plot_qc_metrics(adata, save_dir=plots_dir)

    # Step 4: Apply basic QC filters BEFORE doublet detection
    print("\nApplying initial QC filters for doublet detection...")
    # Keep cells with reasonable QC metrics for doublet detection
    adata_for_doublets = adata[
        (adata.obs.n_genes_by_counts >= CELL_FILTERS["min_genes"])
        & (adata.obs.n_genes_by_counts <= CELL_FILTERS["max_genes"])
        & (adata.obs.percent_mt <= CELL_FILTERS["max_mt_pct"])
    ].copy()
    print(
        f"Cells for doublet detection: {adata_for_doublets.n_obs} (from {adata.n_obs})"
    )

    # Step 5: Detect doublets on QC-filtered cells
    from utils.improved_doublet_detection import detect_doublets_improved

    adata_for_doublets = detect_doublets_improved(
        adata_for_doublets,
        expected_doublet_rate=DOUBLET_PARAMS["expected_doublet_rate"],
        manual_threshold=0.35,  # More sensitive threshold
        plot_histograms=True,
        save_dir=plots_dir,
    )

    # Transfer doublet annotations back to original adata
    adata.obs["doublet_score"] = 0.0
    adata.obs["predicted_doublet"] = False
    adata.obs.loc[adata_for_doublets.obs.index, "doublet_score"] = (
        adata_for_doublets.obs["doublet_score"]
    )
    adata.obs.loc[adata_for_doublets.obs.index, "predicted_doublet"] = (
        adata_for_doublets.obs["predicted_doublet"]
    )

    # Step 6: Filter cells and genes (including doublets)
    adata = filter_cells_and_genes(
        adata,
        min_genes=CELL_FILTERS["min_genes"],
        max_genes=CELL_FILTERS["max_genes"],
        max_mt_pct=CELL_FILTERS["max_mt_pct"],
        min_counts=CELL_FILTERS["min_counts"],
        max_counts=CELL_FILTERS["max_counts"],
        max_ribo_pct=CELL_FILTERS["max_ribo_pct"],
    )

    # Step 7: Normalize and scale
    adata = normalize_and_scale(adata)

    # Step 8: PCA, UMAP, clustering
    adata = run_pca_umap_clustering(adata, save_dir=plots_dir)

    # Step 9: Plot embeddings
    plot_embeddings(adata, save_dir=plots_dir)

    # Step 10: Visualize doublets on UMAP
    from utils.improved_doublet_detection import plot_doublet_scores_umap

    plot_doublet_scores_umap(adata, save_dir=plots_dir)

    # Step 11: Annotate cell types
    adata = annotate_cell_types(
        adata,
        save_dir=plots_dir,
        label_mode=label_mode,
        margin=0.05,
        cluster_agg="median",
    )

    # Step 12: Plot cell type summary
    plot_cell_type_summary(adata, save_dir=plots_dir)

    # Save results
    output_path = "annotated_cellbender_data.h5ad"
    adata.write(output_path)
    print(f"Saved annotated data to {output_path}")

    print("Analysis complete!")
    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="scRNA-seq QC, clustering, and annotation"
    )
    parser.add_argument(
        "--label-mode",
        choices=["cell", "cluster"],
        default="cell",
        help="Cell type labeling mode: 'cell' for per-cell or 'cluster' for cluster-level",
    )
    parser.add_argument(
        "--plots-dir",
        default="plots",
        help="Directory to write plots to (default: 'plots')",
    )
    args = parser.parse_args()

    adata = main(label_mode=args.label_mode, plots_dir_path=args.plots_dir)
