#!/usr/bin/env python3
"""
Single-cell RNA-seq data preprocessing, QC, and cell type annotation
Refactored modular version

This script performs:
1. CellBender data loading and integration
2. Quality control and doublet detection
3. Normalization and dimensionality reduction
4. Clustering and cell type annotation
"""

import warnings
import scanpy as sc

# Import our custom modules
from utils.data_loader import load_and_merge_cellbender_data, add_metadata
from utils.qc_utils import (
    calculate_qc_metrics,
    plot_qc_metrics,
    detect_doublets_scrublet,
    filter_cells_and_genes,
)
from utils.processing import (
    normalize_and_scale,
    run_pca_umap_clustering,
    plot_embeddings,
)
from utils.annotation import annotate_cell_types, plot_cell_type_summary

# Configure scanpy
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor="white")

# Suppress warnings
warnings.filterwarnings("ignore")


def main():
    """Main analysis pipeline"""
    print("Starting single-cell analysis pipeline...")

    # Define parameters
    base_path = "/Users/melissadu/Documents/tsai-lab-urop/base/"
    sample_names = [f"D25-{i}" for i in range(2675, 2691)]
    custom_name = "_processed_feature_bc_matrix_filtered.h5"

    # Step 1: Load and merge data
    adata = load_and_merge_cellbender_data(base_path, sample_names, custom_name)

    # Step 2: Add metadata
    adata = add_metadata(adata, sample_names)

    # Step 3: Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    plot_qc_metrics(adata)

    # Step 4: Detect doublets
    adata = detect_doublets_scrublet(adata)

    # Step 5: Filter cells and genes
    adata = filter_cells_and_genes(adata)

    # Step 6: Normalize and scale
    adata = normalize_and_scale(adata)

    # Step 7: PCA, UMAP, clustering
    adata = run_pca_umap_clustering(adata)

    # Step 8: Plot embeddings
    plot_embeddings(adata)

    # Step 9: Annotate cell types
    adata = annotate_cell_types(adata)

    # Step 10: Plot cell type summary
    plot_cell_type_summary(adata)

    # Save results
    output_path = "annotated_cellbender_data.h5ad"
    adata.write(output_path)
    print(f"Saved annotated data to {output_path}")

    print("Analysis complete!")
    return adata


if __name__ == "__main__":
    adata = main()
