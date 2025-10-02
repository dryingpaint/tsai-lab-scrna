#!/usr/bin/env python3
"""
Differential expression analysis utilities for single-cell RNA-seq analysis
Handles pseudobulk creation and statistical testing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


def create_condition_column(adata):
    """Create condition column combining genotype and stimulation

    Args:
        adata: AnnData object

    Returns:
        AnnData object with condition column added
    """
    adata.obs["condition"] = (
        adata.obs["Genotype"].astype(str) + "_" + adata.obs["Stimulation"].astype(str)
    )

    # Set factor levels
    condition_order = ["E3_Ctrl", "E3_GENUS", "E4_Ctrl", "E4_GENUS"]
    adata.obs["condition"] = pd.Categorical(
        adata.obs["condition"], categories=condition_order, ordered=True
    )

    return adata


def create_pseudobulk(adata, groupby=["orig.ident", "celltype"], min_cells=10):
    """Create pseudobulk samples by aggregating cells

    Args:
        adata: AnnData object
        groupby: Columns to group by for pseudobulk creation
        min_cells: Minimum cells required per pseudobulk sample

    Returns:
        Tuple of (pseudobulk_df, sample_info_df)
    """
    print("Creating pseudobulk samples...")

    # Create group identifier
    adata.obs["group_id"] = (
        adata.obs["orig.ident"].astype(str) + "--" + adata.obs["celltype"].astype(str)
    )

    # Get raw counts
    if adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    # Aggregate by group
    pseudobulk_data = []
    sample_info = []

    for group_id in adata.obs["group_id"].unique():
        mask = adata.obs["group_id"] == group_id
        n_cells = mask.sum()

        if n_cells >= min_cells:
            # Sum counts across cells
            if hasattr(X, "toarray"):
                group_counts = X[mask].toarray().sum(axis=0)
            else:
                group_counts = X[mask].sum(axis=0)

            pseudobulk_data.append(group_counts)

            # Get sample metadata
            sample_meta = adata.obs[mask].iloc[0]
            sample_info.append(
                {
                    "group_id": group_id,
                    "sample_id": sample_meta["orig.ident"],
                    "celltype": sample_meta["celltype"],
                    "Genotype": sample_meta["Genotype"],
                    "Stimulation": sample_meta["Stimulation"],
                    "Sex": sample_meta["Sex"],
                    "condition": sample_meta["condition"],
                    "n_cells": n_cells,
                }
            )

    # Create pseudobulk matrix
    pb_matrix = np.array(pseudobulk_data).T  # genes x samples
    pb_df = pd.DataFrame(
        pb_matrix, index=var_names, columns=[info["group_id"] for info in sample_info]
    )

    sample_info_df = pd.DataFrame(sample_info)

    print(f"Created {pb_df.shape[1]} pseudobulk samples from {pb_df.shape[0]} genes")

    return pb_df, sample_info_df


def filter_genes_for_de(pb_df, sample_info_df, min_count=5, min_samples=2):
    """Filter genes for differential expression analysis

    Args:
        pb_df: Pseudobulk expression DataFrame
        sample_info_df: Sample metadata DataFrame
        min_count: Minimum count threshold
        min_samples: Minimum number of samples

    Returns:
        Filtered pseudobulk DataFrame
    """
    print("Filtering genes for DE analysis...")

    # Keep genes expressed above threshold in minimum number of samples
    expressed_mask = (pb_df >= min_count).sum(axis=1) >= min_samples
    pb_filtered = pb_df.loc[expressed_mask]

    print(f"Kept {pb_filtered.shape[0]} genes after filtering")

    return pb_filtered


def run_differential_expression(pb_df, sample_info_df, cell_type):
    """Run differential expression analysis for a specific cell type

    Args:
        pb_df: Pseudobulk expression DataFrame
        sample_info_df: Sample metadata DataFrame
        cell_type: Cell type to analyze

    Returns:
        DataFrame with differential expression results
    """
    print(f"Running DE analysis for {cell_type}...")

    # Subset to cell type
    ct_mask = sample_info_df["celltype"] == cell_type
    if ct_mask.sum() < 4:
        print(f"Not enough samples for {cell_type}")
        return None

    ct_samples = sample_info_df[ct_mask]
    ct_counts = pb_df[ct_samples["group_id"]]

    # Filter genes
    ct_counts = filter_genes_for_de(ct_counts, ct_samples)

    if ct_counts.shape[0] < 100:
        print(f"Not enough genes for {cell_type}")
        return None

    # Normalize (simple CPM normalization)
    lib_sizes = ct_counts.sum(axis=0)
    ct_cpm = ct_counts.div(lib_sizes, axis=1) * 1e6
    ct_log_cpm = np.log2(ct_cpm + 1)

    results = []

    # Define contrasts
    contrasts = [
        ("E3_GENUS_vs_Ctrl", "E3_GENUS", "E3_Ctrl"),
        ("E4_GENUS_vs_Ctrl", "E4_GENUS", "E4_Ctrl"),
        ("E4_vs_E3_Ctrl", "E4_Ctrl", "E3_Ctrl"),
        ("E4_vs_E3_GENUS", "E4_GENUS", "E3_GENUS"),
    ]

    for contrast_name, group1, group2 in contrasts:
        print(f"  Testing {contrast_name}")

        mask1 = ct_samples["condition"] == group1
        mask2 = ct_samples["condition"] == group2

        if mask1.sum() == 0 or mask2.sum() == 0:
            print(f"    Skipping {contrast_name} - missing groups")
            continue

        group1_data = ct_log_cpm.loc[:, ct_samples.loc[mask1, "group_id"]]
        group2_data = ct_log_cpm.loc[:, ct_samples.loc[mask2, "group_id"]]

        # Simple t-test for each gene
        gene_results = []
        for gene in ct_log_cpm.index:
            try:
                stat, pval = stats.ttest_ind(
                    group1_data.loc[gene].values,
                    group2_data.loc[gene].values,
                    equal_var=False,
                )

                mean1 = group1_data.loc[gene].mean()
                mean2 = group2_data.loc[gene].mean()
                logfc = mean1 - mean2

                gene_results.append(
                    {
                        "gene": gene,
                        "logFC": logfc,
                        "P.Value": pval,
                        "cell_type": cell_type,
                        "contrast": contrast_name,
                    }
                )
            except:
                continue

        if gene_results:
            contrast_df = pd.DataFrame(gene_results)
            # Multiple testing correction
            contrast_df["adj.P.Val"] = multipletests(
                contrast_df["P.Value"], method="fdr_bh"
            )[1]
            results.append(contrast_df)

    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return None


def plot_de_summary(de_results):
    """Plot summary of differential expression results

    Args:
        de_results: DataFrame with DE results

    Returns:
        DataFrame with counts summary
    """
    print("Plotting DE summary...")

    # Count significant genes
    sig_genes = de_results[
        (de_results["adj.P.Val"] < 0.05) & (de_results["logFC"].abs() > 0.5)
    ]

    # Count by cell type and contrast
    counts = sig_genes.groupby(["cell_type", "contrast"]).size().reset_index()
    counts.columns = ["cell_type", "contrast", "n_genes"]

    # Pivot for heatmap
    heatmap_data = counts.pivot(
        index="cell_type", columns="contrast", values="n_genes"
    ).fillna(0)

    # Plot heatmap
    plt.figure(figsize=(12, 6))
    sns.heatmap(heatmap_data, annot=True, fmt="g", cmap="Blues")
    plt.title("Number of significant DE genes (adj.P < 0.05, |logFC| > 0.5)")
    plt.xlabel("Contrast")
    plt.ylabel("Cell type")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

    return counts
