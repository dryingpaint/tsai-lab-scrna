#!/usr/bin/env python3
"""
Differential gene expression analysis and pathway enrichment
Refactored modular version

This script performs:
1. Module score analysis for pathway activity
2. Pseudobulk differential expression analysis
3. Gene set enrichment analysis (GSEA)
"""

import warnings
from pathlib import Path
import scanpy as sc

# Import our custom modules
from utils.pathway_anwalysis import (
    load_gene_sets,
    calculate_module_scores,
    plot_module_scores,
    run_gsea_analysis,
    plot_gsea_results,
)
from utils.differential_expression import (
    create_condition_column,
    create_pseudobulk,
    run_differential_expression,
    plot_de_summary,
)

# Configure
sc.settings.verbosity = 1
warnings.filterwarnings("ignore")


def main():
    """Main analysis pipeline"""
    print("Starting differential expression and GSEA analysis...")

    # Load annotated data
    adata_path = "annotated_cellbender_data.h5ad"
    if not Path(adata_path).exists():
        print(
            f"Error: {adata_path} not found. Run 1_cellbender_qc_annotation.py first."
        )
        return None, None, None

    adata = sc.read_h5ad(adata_path)
    print(f"Loaded data: {adata.n_obs} cells, {adata.n_vars} genes")

    # Step 1: Load gene sets
    gene_sets = load_gene_sets()

    # Step 2: Create condition column
    adata = create_condition_column(adata)

    # Step 3: Calculate module scores
    adata = calculate_module_scores(adata, gene_sets)

    # Plot module scores
    score_names = [col for col in adata.obs.columns if col.endswith("_score")]
    if score_names:
        plot_module_scores(adata, score_names)

    # Step 4: Create pseudobulk samples
    pb_df, sample_info_df = create_pseudobulk(adata)

    # Step 5: Run differential expression analysis
    de_results_list = []
    cell_types = ["AST", "ExN", "InN_SST", "InN_VIP", "ODC", "OPC", "Endo"]

    for cell_type in cell_types:
        de_result = run_differential_expression(pb_df, sample_info_df, cell_type)
        if de_result is not None:
            de_results_list.append(de_result)

    if de_results_list:
        all_de_results = pd.concat(de_results_list, ignore_index=True)

        # Save DE results
        all_de_results.to_csv("differential_expression_results.csv", index=False)
        print("Saved DE results to differential_expression_results.csv")

        # Plot DE summary
        plot_de_summary(all_de_results)

        # Step 6: Run GSEA
        gsea_results = run_gsea_analysis(all_de_results, gene_sets)

        if not gsea_results.empty:
            # Save GSEA results
            gsea_results.to_csv("gsea_results.csv", index=False)
            print("Saved GSEA results to gsea_results.csv")

            # Plot example GSEA results
            plot_gsea_results(gsea_results, "AST", "E4_GENUS_vs_Ctrl")

        print("Analysis complete!")
        return adata, all_de_results, gsea_results
    else:
        print("No differential expression results generated")
        return adata, None, None


if __name__ == "__main__":
    import pandas as pd  # Import here to avoid circular imports

    adata, de_results, gsea_results = main()
