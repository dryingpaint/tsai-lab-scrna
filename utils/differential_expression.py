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

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


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
    # E3_Ctrl <-> E4_Ctrl
    # E3_GENUS <-> E4_GENUS
    # E3_Ctrl <-> E3_GENUS
    # E4_Ctrl <-> E4_GENUS
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
            # Convert pandas Series mask to numpy array for sparse matrix indexing
            mask_array = mask.values
            if hasattr(X, "toarray"):
                group_counts = X[mask_array].toarray().sum(axis=0)
            else:
                group_counts = X[mask_array].sum(axis=0)
            
            # Ensure counts are 1D array
            if hasattr(group_counts, "A1"):
                group_counts = group_counts.A1  # For sparse matrix results
            elif len(group_counts.shape) > 1:
                group_counts = group_counts.flatten()

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

def run_de_with_deseq2(counts_df, sample_info_df, contrast_name, group1, group2, 
                        de_params, cell_type):
    """Run DESeq2 differential expression for a single contrast
    
    Args:
        counts_df: Count matrix (genes × samples)
        sample_info_df: Sample metadata DataFrame
        contrast_name: Name of the contrast
        group1: First condition
        group2: Second condition (reference)
        de_params: DE parameters dictionary
        cell_type: Cell type being analyzed
        
    Returns:
        DataFrame with DE results or None
    """
    # Subset to samples in this contrast
    mask = sample_info_df['condition'].isin([group1, group2])
    contrast_samples = sample_info_df[mask].copy()
    contrast_counts = counts_df[contrast_samples['group_id']].copy()
    
    if len(contrast_samples) < 4:
        print(f"  ⚠️  Skipping {contrast_name}: Only {len(contrast_samples)} samples")
        return None
    
    print(f"  Testing {contrast_name} ({(contrast_samples['condition']==group1).sum()} vs {(contrast_samples['condition']==group2).sum()} samples)")
    print(f"    Input: {contrast_counts.shape[0]} genes × {contrast_counts.shape[1]} samples")
    
    # Prepare metadata - ensure condition is categorical with reference level
    contrast_samples['condition'] = pd.Categorical(
        contrast_samples['condition'],
        categories=[group2, group1],  # Reference (group2) first
        ordered=False
    )
    
    # Reorder counts to match metadata
    contrast_counts = contrast_counts[contrast_samples['group_id']]
    
    # Convert counts to integer (DESeq2 requires integer counts)
    counts_int = np.round(contrast_counts.values).astype(int)
    
    # PyDESeq2 expects counts as samples × genes (samples as rows)
    # Our data is genes × samples, so we need to transpose
    counts_transposed = pd.DataFrame(
        counts_int.T,  # Transpose: samples × genes
        index=contrast_counts.columns,  # Sample IDs as row index
        columns=contrast_counts.index   # Gene names as column names
    )
    
    # Prepare metadata with matching index
    metadata_indexed = contrast_samples.set_index('group_id')
    
    # Verify dimensions match
    assert counts_transposed.shape[0] == len(metadata_indexed), \
        f"Count matrix has {counts_transposed.shape[0]} samples but metadata has {len(metadata_indexed)} samples"
    assert all(counts_transposed.index == metadata_indexed.index), \
        "Sample IDs in counts and metadata do not match"
    
    try:
        # Create DESeq2 dataset
        dds = DeseqDataSet(
            counts=counts_transposed,
            metadata=metadata_indexed,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=1  # Single CPU for stability
        )
        
        # Run DESeq2
        dds.deseq2()
        
        # Get results
        stat_res = DeseqStats(dds, contrast=["condition", group1, group2])
        stat_res.summary()
        
        # Extract results
        results_df = stat_res.results_df.copy()
        
        # Rename columns to match expected format
        results_df = results_df.rename(columns={
            'log2FoldChange': 'logFC',
            'pvalue': 'P.Value',
            'padj': 'adj.P.Val',
            'baseMean': 'AveExpr'
        })
        
        # Add metadata
        results_df['gene'] = results_df.index
        results_df['cell_type'] = cell_type
        results_df['contrast'] = contrast_name
        
        # Add significance flags
        results_df['significant'] = (
            (results_df['adj.P.Val'] < de_params['fdr_threshold']) &
            (results_df['logFC'].abs() > de_params['fc_threshold']) &
            (results_df['adj.P.Val'].notna())
        )
        
        results_df['upregulated'] = (
            results_df['significant'] & (results_df['logFC'] > 0)
        )
        results_df['downregulated'] = (
            results_df['significant'] & (results_df['logFC'] < 0)
        )
        
        # Count significant genes
        n_sig = results_df['significant'].sum()
        n_up = results_df['upregulated'].sum()
        n_down = results_df['downregulated'].sum()
        
        print(f"    ✓ {n_sig} significant genes ({n_up} up, {n_down} down)")
        
        return results_df[['gene', 'logFC', 'P.Value', 'adj.P.Val', 'AveExpr',
                          'cell_type', 'contrast', 'significant', 'upregulated', 'downregulated']]
        
    except Exception as e:
        print(f"  ✗ Error running DESeq2: {str(e)}")
        print(f"     Debug: counts shape = {counts_transposed.shape} (samples × genes)")
        print(f"     Debug: metadata shape = {metadata_indexed.shape}")
        print(f"     Try using use_deseq2=False to use t-test method instead")
        return None


def run_de_with_ttest(counts_df, sample_info_df, contrast_name, group1, group2,
                      de_params, cell_type):
    """Fallback t-test method (not recommended for RNA-seq)
    
    Args:
        counts_df: Count matrix (genes × samples)
        sample_info_df: Sample metadata DataFrame
        contrast_name: Name of the contrast
        group1: First condition
        group2: Second condition (reference)
        de_params: DE parameters dictionary
        cell_type: Cell type being analyzed
        
    Returns:
        DataFrame with DE results or None
    """
    mask1 = sample_info_df['condition'] == group1
    mask2 = sample_info_df['condition'] == group2
    
    if mask1.sum() == 0 or mask2.sum() == 0:
        print(f"  ⚠️  Skipping {contrast_name}: Missing groups")
        return None
    
    print(f"  Testing {contrast_name} ({mask1.sum()} vs {mask2.sum()} samples) [t-test]")
    
    # Normalize to CPM and log-transform
    lib_sizes = counts_df.sum(axis=0)
    ct_cpm = counts_df.div(lib_sizes, axis=1) * 1e6
    ct_log_cpm = np.log2(ct_cpm + 1)
    
    group1_data = ct_log_cpm.loc[:, sample_info_df.loc[mask1, 'group_id']]
    group2_data = ct_log_cpm.loc[:, sample_info_df.loc[mask2, 'group_id']]
    
    # Run t-test for each gene
    gene_results = []
    for gene in ct_log_cpm.index:
        try:
            # Welch's t-test (unequal variances)
            stat, pval = stats.ttest_ind(
                group1_data.loc[gene].values,
                group2_data.loc[gene].values,
                equal_var=False
            )
            
            mean1 = group1_data.loc[gene].mean()
            mean2 = group2_data.loc[gene].mean()
            logfc = mean1 - mean2
            
            gene_results.append({
                'gene': gene,
                'logFC': logfc,
                'P.Value': pval,
                'AveExpr': (mean1 + mean2) / 2,
                'cell_type': cell_type,
                'contrast': contrast_name,
            })
        except:
            continue
    
    if gene_results:
        contrast_df = pd.DataFrame(gene_results)
        # Multiple testing correction (FDR)
        contrast_df['adj.P.Val'] = multipletests(
            contrast_df['P.Value'],
            method='fdr_bh'
        )[1]
        
        # Add significance flags
        contrast_df['significant'] = (
            (contrast_df['adj.P.Val'] < de_params['fdr_threshold']) &
            (contrast_df['logFC'].abs() > de_params['fc_threshold'])
        )
        
        contrast_df['upregulated'] = (
            contrast_df['significant'] & (contrast_df['logFC'] > 0)
        )
        contrast_df['downregulated'] = (
            contrast_df['significant'] & (contrast_df['logFC'] < 0)
        )
        
        n_sig = contrast_df['significant'].sum()
        n_up = contrast_df['upregulated'].sum()
        n_down = contrast_df['downregulated'].sum()
        
        print(f"    ✓ {n_sig} significant genes ({n_up} up, {n_down} down)")
        
        return contrast_df
    else:
        return None


def run_de_for_celltype(pb_df, sample_info_df, cell_type, de_params, 
                        min_samples_per_group=2, use_deseq2=True):
    """Run differential expression analysis for a specific cell type
    
    Args:
        pb_df: Pseudobulk expression DataFrame (genes × samples)
        sample_info_df: Sample metadata DataFrame
        cell_type: Cell type to analyze
        de_params: Dictionary of DE parameters (should include 'min_count', 'min_samples_expr', 'fdr_threshold', 'fc_threshold')
        min_samples_per_group: Minimum number of samples per group required (default: 2)
        use_deseq2: Whether to use DESeq2 (True) or t-test fallback (False)
    
    Returns:
        DataFrame with DE results
    """
    print(f"\n{'='*60}")
    print(f"ANALYZING: {cell_type}")
    print(f"{'='*60}")
    
    # Subset to cell type
    ct_mask = sample_info_df['celltype'] == cell_type
    ct_samples = sample_info_df[ct_mask].copy()
    
    if len(ct_samples) < min_samples_per_group * 2:
        print(f"⚠️  Skipping {cell_type}: Only {len(ct_samples)} samples")
        return None
    
    # Get counts for this cell type
    ct_counts = pb_df[ct_samples['group_id']].copy()
    
    # Filter genes
    ct_counts = filter_genes_for_de(
        ct_counts,
        ct_samples,
        min_count=de_params['min_count'],
        min_samples=de_params['min_samples_expr']
    )
    
    if ct_counts.shape[0] < 100:
        print(f"⚠️  Skipping {cell_type}: Only {ct_counts.shape[0]} genes after filtering")
        return None
    
    print(f"  Analyzing {ct_counts.shape[0]:,} genes across {len(ct_samples)} samples")
    if use_deseq2:
        print(f"  Method: DESeq2 (negative binomial model)")
    else:
        print(f"  Method: t-test on log-CPM (fallback)")
    
    # Define contrasts
    contrasts = [
        ('E3_GENUS_vs_Ctrl', 'E3_GENUS', 'E3_Ctrl'),
        ('E4_GENUS_vs_Ctrl', 'E4_GENUS', 'E4_Ctrl'),
        ('E4_vs_E3_Ctrl', 'E4_Ctrl', 'E3_Ctrl'),
        ('E4_vs_E3_GENUS', 'E4_GENUS', 'E3_GENUS'),
    ]
    
    results = []
    
    for contrast_name, group1, group2 in contrasts:
        if use_deseq2:
            result = run_de_with_deseq2(
                ct_counts, ct_samples, contrast_name, 
                group1, group2, de_params, cell_type
            )
        else:
            result = run_de_with_ttest(
                ct_counts, ct_samples, contrast_name,
                group1, group2, de_params, cell_type
            )
        
        if result is not None:
            results.append(result)
    
    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return None


def plot_de_heatmap(pb_df, sample_info_df, de_results, cell_type, contrast, 
                    top_n=50, save_path=None):
    """Plot heatmap of top DE genes
    
    Args:
        pb_df: Pseudobulk expression DataFrame
        sample_info_df: Sample metadata DataFrame
        de_results: DE results DataFrame
        cell_type: Cell type to plot
        contrast: Contrast name
        top_n: Number of top genes to show
        save_path: Path to save figure
    """
    # Get results for this cell type and contrast
    ct_results = de_results[
        (de_results['cell_type'] == cell_type) &
        (de_results['contrast'] == contrast)
    ].copy()
    
    if len(ct_results) == 0:
        print(f"No results for {cell_type} - {contrast}")
        return
    
    # Get top upregulated and downregulated genes
    top_up = ct_results.nlargest(top_n // 2, 'logFC')
    top_down = ct_results.nsmallest(top_n // 2, 'logFC')
    top_genes = pd.concat([top_up, top_down])['gene'].tolist()
    
    # Get samples for this cell type
    ct_samples = sample_info_df[sample_info_df['celltype'] == cell_type]
    
    # Get expression data (log CPM)
    lib_sizes = pb_df[ct_samples['group_id']].sum(axis=0)
    ct_cpm = pb_df[ct_samples['group_id']].div(lib_sizes, axis=1) * 1e6
    ct_log_cpm = np.log2(ct_cpm + 1)
    
    # Subset to top genes
    heatmap_data = ct_log_cpm.loc[top_genes]
    
    # Order samples by condition
    condition_order = ['E3_Ctrl', 'E3_GENUS', 'E4_Ctrl', 'E4_GENUS']
    sample_order = []
    for cond in condition_order:
        cond_samples = ct_samples[ct_samples['condition'] == cond]['group_id'].tolist()
        sample_order.extend([s for s in cond_samples if s in heatmap_data.columns])
    
    heatmap_data = heatmap_data[sample_order]
    
    # Create annotation for samples
    sample_annot = pd.DataFrame({
        'Condition': [ct_samples[ct_samples['group_id'] == s]['condition'].iloc[0] 
                     for s in sample_order]
    }, index=sample_order)
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, max(8, len(top_genes) * 0.3)))
    
    # Color map
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    
    # Plot heatmap
    sns.heatmap(
        heatmap_data,
        cmap=cmap,
        center=heatmap_data.values.mean(),
        xticklabels=False,
        yticklabels=True,
        cbar_kws={'label': 'Log2(CPM + 1)'},
        ax=ax
    )
    
    ax.set_title(f'{cell_type} - {contrast}\nTop {len(top_genes)} DE genes', 
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Samples (grouped by condition)', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    
    plt.show()


def plot_volcano(de_results, cell_type, contrast, fc_threshold=0.5,
                 pval_threshold=0.05, save_path=None):
    """Plot volcano plot for DE results
    
    Args:
        de_results: DE results DataFrame
        cell_type: Cell type to plot
        contrast: Contrast name
        fc_threshold: Log2FC threshold for coloring
        pval_threshold: P-value threshold for coloring
        save_path: Path to save figure
    """
    # Get results for this cell type and contrast
    ct_results = de_results[
        (de_results['cell_type'] == cell_type) &
        (de_results['contrast'] == contrast)
    ].copy()
    
    if len(ct_results) == 0:
        print(f"No results for {cell_type} - {contrast}")
        return
    
    # Calculate -log10(p-value)
    ct_results['neg_log10_pval'] = -np.log10(ct_results['P.Value'] + 1e-300)
    
    # Categorize genes
    ct_results['category'] = 'Not significant'
    ct_results.loc[
        (ct_results['adj.P.Val'] < pval_threshold) & 
        (ct_results['logFC'] > fc_threshold),
        'category'
    ] = 'Upregulated'
    ct_results.loc[
        (ct_results['adj.P.Val'] < pval_threshold) & 
        (ct_results['logFC'] < -fc_threshold),
        'category'
    ] = 'Downregulated'
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot non-significant genes
    ns_data = ct_results[ct_results['category'] == 'Not significant']
    ax.scatter(ns_data['logFC'], ns_data['neg_log10_pval'], 
               c='gray', alpha=0.5, s=20, label='Not significant')
    
    # Plot upregulated genes
    up_data = ct_results[ct_results['category'] == 'Upregulated']
    if len(up_data) > 0:
        ax.scatter(up_data['logFC'], up_data['neg_log10_pval'],
                  c='red', alpha=0.7, s=30, label=f'Upregulated (n={len(up_data)})')
    
    # Plot downregulated genes
    down_data = ct_results[ct_results['category'] == 'Downregulated']
    if len(down_data) > 0:
        ax.scatter(down_data['logFC'], down_data['neg_log10_pval'],
                  c='blue', alpha=0.7, s=30, label=f'Downregulated (n={len(down_data)})')
    
    # Add threshold lines
    ax.axvline(fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(-fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axhline(-np.log10(pval_threshold), color='black', linestyle='--', 
              linewidth=1, alpha=0.5)
    
    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10(P-value)', fontsize=12)
    ax.set_title(f'{cell_type} - {contrast}\nVolcano Plot', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    
    plt.show()
