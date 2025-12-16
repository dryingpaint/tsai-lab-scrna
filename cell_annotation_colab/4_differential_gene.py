# %% [markdown]
# # Notebook 4: Differential Gene Expression Analysis
#
# **Cell Annotation Pipeline - Part 4 of 4**
#
# **📥 Input:** `outputs/annotated_data.h5ad`  
# **📤 Output:** `outputs/differential_expression_results/`
#
# ---
#
# ## Overview
#
# This notebook performs **pseudobulk differential expression analysis** to identify genes that are differentially expressed between experimental conditions (e.g., E3 vs E4, Ctrl vs GENUS) within specific cell types.
#
# **Key Steps:**
# 1. Load annotated data
# 2. Create pseudobulk samples (aggregate cells per donor × cell type)
# 3. Filter cell types (focus on major types, exclude mixed clusters)
# 4. Run differential expression analysis for each cell type
# 5. Visualize results (heatmaps, volcano plots)
# 6. Export results for downstream analysis
#
# **Method:** Pseudobulk aggregation followed by statistical testing (t-test/Wilcoxon). For more robust analysis, results can be exported to R for DESeq2 analysis.
#
# ---

# %% [markdown]
# ## 1. Setup & Load Data
#
# Load the annotated data from Notebook 3 and set up the analysis environment.

# %%
# Verify PyDESeq2 is installed and check module version
import utils.differential_expression as de_module
import scanpy as sc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from utils.differential_expression import (
    create_condition_column,
    create_pseudobulk,
    plot_de_summary,
)


print("✓ Differential expression module loaded from:")
print(f"  {de_module.__file__}")
print()


print()

# Check that run_de_for_celltype exists with PyDESeq2 support
if hasattr(de_module, 'run_de_for_celltype'):
    print("✓ run_de_for_celltype() function available")
    print("✓ Ready for PyDESeq2-based differential expression analysis!")

OUTPUT_DIR = Path('outputs/differential_expression_results/')

# %%
%load_ext autoreload
%autoreload 2

# %%
# Load annotated data
print("Loading annotated data from Notebook 3...")
adata = sc.read_h5ad('outputs/annotated_data.h5ad')

# Validate required columns
required_cols = ['celltype', 'orig.ident', 'Genotype', 'Stimulation', 'leiden']
missing = [c for c in required_cols if c not in adata.obs.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

print(f"\n✓ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  Cell types: {adata.obs['celltype'].nunique()}")
print(f"  Samples: {adata.obs['orig.ident'].nunique()}")
print(f"  Clusters: {adata.obs['leiden'].nunique()}")

# Check if raw data is available
if adata.raw is None:
    print("\n⚠️  Warning: No raw data found. Using normalized data.")
    print("   For best results, use raw counts from Notebook 1.")
else:
    print(f"\n✓ Raw data available: {adata.raw.n_vars:,} genes")

# ============================================================================
# CLUSTER PURITY FILTERING
# ============================================================================
print("\n" + "="*60)
print("CLUSTER PURITY ANALYSIS")
print("="*60)

# Calculate cell type composition per cluster
cluster_composition = pd.crosstab(
    adata.obs['leiden'],
    adata.obs['celltype'],
    normalize='index'  # Normalize by cluster (rows sum to 1)
)

# Find dominant cell type and its proportion for each cluster
dominant_celltype = cluster_composition.idxmax(axis=1)
dominant_proportion = cluster_composition.max(axis=1)

# Create cluster purity dataframe
cluster_purity = pd.DataFrame({
    'cluster': cluster_composition.index,
    'dominant_celltype': dominant_celltype.values,
    'purity': dominant_proportion.values,
    'n_cells': adata.obs['leiden'].value_counts()[cluster_composition.index].values
})

print("\nCluster purity (% of dominant cell type):")
print(cluster_purity.to_string(index=False))

# Define purity threshold
PURITY_THRESHOLD = 0.50  # 🔧 Minimum proportion of dominant cell type

# Filter to pure clusters
pure_clusters = cluster_purity[cluster_purity['purity'] >= PURITY_THRESHOLD]['cluster'].tolist()
mixed_clusters = cluster_purity[cluster_purity['purity'] < PURITY_THRESHOLD]['cluster'].tolist()

print(f"\n✓ Pure clusters (≥{PURITY_THRESHOLD*100:.0f}% one cell type): {len(pure_clusters)} clusters")
print(f"  Clusters: {pure_clusters}")

if mixed_clusters:
    print(f"\n⚠️  Mixed clusters (<{PURITY_THRESHOLD*100:.0f}% any cell type): {len(mixed_clusters)} clusters")
    print(f"  Clusters: {mixed_clusters}")
    print("  These will be EXCLUDED from differential expression analysis")

# Filter adata to only include cells from pure clusters
n_cells_before = adata.n_obs
adata = adata[adata.obs['leiden'].isin(pure_clusters)].copy()
n_cells_after = adata.n_obs
n_cells_removed = n_cells_before - n_cells_after

print(f"\n✓ Filtered to pure clusters:")
print(f"  Cells retained: {n_cells_after:,} ({n_cells_after/n_cells_before*100:.1f}%)")
print(f"  Cells removed: {n_cells_removed:,} ({n_cells_removed/n_cells_before*100:.1f}%)")

# Show cell type distribution after filtering
print(f"\n✓ Cell type distribution (after filtering):")
print(adata.obs['celltype'].value_counts().sort_index())

# Save cluster purity info
cluster_purity.to_csv(OUTPUT_DIR / 'cluster_purity.csv', index=False)
print(f"\n✓ Saved cluster purity info to: {OUTPUT_DIR / 'cluster_purity.csv'}")

# %% [markdown]
# ## 2. Visualize Cluster Purity
#
# Visualize which clusters are pure enough for analysis.

# %%
# Visualize cluster purity
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Plot 1: Barplot of cluster purity
ax = axes[0]
colors = ['green' if p >= PURITY_THRESHOLD else 'red' for p in cluster_purity['purity']]
bars = ax.bar(cluster_purity['cluster'].astype(str), cluster_purity['purity'], color=colors, alpha=0.7)
ax.axhline(PURITY_THRESHOLD, color='black', linestyle='--', linewidth=2, label=f'Threshold ({PURITY_THRESHOLD*100:.0f}%)')
ax.set_xlabel('Cluster', fontsize=12)
ax.set_ylabel('Purity (proportion of dominant cell type)', fontsize=12)
ax.set_title('Cluster Purity Analysis', fontsize=14, fontweight='bold')
ax.set_ylim(0, 1.0)
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Rotate x-axis labels if many clusters
if len(cluster_purity) > 15:
    ax.tick_params(axis='x', rotation=45)

# Plot 2: Show which cell type dominates each cluster
ax = axes[1]
cluster_purity_sorted = cluster_purity.sort_values('purity', ascending=True)
colors = ['green' if p >= PURITY_THRESHOLD else 'red' for p in cluster_purity_sorted['purity']]
bars = ax.barh(cluster_purity_sorted['cluster'].astype(str), 
               cluster_purity_sorted['purity'], 
               color=colors, alpha=0.7)
ax.axvline(PURITY_THRESHOLD, color='black', linestyle='--', linewidth=2, label=f'Threshold ({PURITY_THRESHOLD*100:.0f}%)')
ax.set_ylabel('Cluster', fontsize=12)
ax.set_xlabel('Purity (proportion of dominant cell type)', fontsize=12)
ax.set_title('Clusters by Dominant Cell Type', fontsize=14, fontweight='bold')
ax.set_xlim(0, 1.0)
ax.legend()
ax.grid(axis='x', alpha=0.3)

# Add cell type labels
for i, (cluster, celltype, purity) in enumerate(zip(cluster_purity_sorted['cluster'], 
                                                     cluster_purity_sorted['dominant_celltype'],
                                                     cluster_purity_sorted['purity'])):
    ax.text(purity + 0.02, i, f'{celltype}', va='center', fontsize=9)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'cluster_purity_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"✓ Saved: {OUTPUT_DIR / 'cluster_purity_analysis.png'}")

# Create a summary table
print("\n" + "="*60)
print("CLUSTER FILTERING SUMMARY")
print("="*60)
print(f"\n{'Status':<15} {'Clusters':<10} {'Cells':<10} {'% Cells':<10}")
print("-" * 50)

pure_cluster_cells = cluster_purity[cluster_purity['purity'] >= PURITY_THRESHOLD]['n_cells'].sum()
mixed_cluster_cells = cluster_purity[cluster_purity['purity'] < PURITY_THRESHOLD]['n_cells'].sum()
total_cells = pure_cluster_cells + mixed_cluster_cells

print(f"{'Pure (≥50%)':<15} {len(pure_clusters):<10} {pure_cluster_cells:<10} {pure_cluster_cells/total_cells*100:>8.1f}%")
print(f"{'Mixed (<50%)':<15} {len(mixed_clusters):<10} {mixed_cluster_cells:<10} {mixed_cluster_cells/total_cells*100:>8.1f}%")
print("-" * 50)
print(f"{'Total':<15} {len(cluster_purity):<10} {total_cells:<10} {'100.0%':>8}")

# Save purity visualization data
cluster_purity.to_csv(OUTPUT_DIR / 'cluster_purity_detailed.csv', index=False)
print(f"\n✓ Saved detailed cluster purity to: {OUTPUT_DIR / 'cluster_purity_detailed.csv'}")

# %% [markdown]
# ## 3. Parameter Configuration
#
# Configure which cell types to analyze and set analysis parameters.

# %%
# ============================================================================
# CELL TYPE SELECTION
# ============================================================================

# Automatically use all cell types present in the filtered data
# These come from pure clusters only (>50% one cell type)
celltype_counts = adata.obs['celltype'].value_counts()
CELL_TYPES_TO_ANALYZE = celltype_counts.index.tolist()

print("Cell types available (from pure clusters):")
for ct in CELL_TYPES_TO_ANALYZE:
    n_cells = celltype_counts[ct]
    n_samples = adata.obs[adata.obs['celltype'] == ct]['orig.ident'].nunique()
    print(f"  {ct}: {n_cells:,} cells across {n_samples} samples")

# Optional: Filter to only major cell types (uncomment to use)
# MAJOR_CELL_TYPES = ['Excit', 'Inhib', 'Astro', 'Oligo', 'OPC', 'Micro', 'Endo']
# CELL_TYPES_TO_ANALYZE = [ct for ct in CELL_TYPES_TO_ANALYZE if ct in MAJOR_CELL_TYPES]

# Optional: Filter by minimum cell count (uncomment to use)
MIN_CELLS_PER_TYPE = 100  # 🔧 Minimum cells per cell type
CELL_TYPES_TO_ANALYZE = [ct for ct in CELL_TYPES_TO_ANALYZE if celltype_counts[ct] >= MIN_CELLS_PER_TYPE]

print(f"\n✓ Will analyze {len(CELL_TYPES_TO_ANALYZE)} cell types (≥{MIN_CELLS_PER_TYPE} cells each)")

# ============================================================================
# PSEUDOBULK PARAMETERS
# ============================================================================

PSEUDOBULK_PARAMS = {
    'min_cells': 10,        # 🔧 Minimum cells per pseudobulk sample
    'min_samples': 2,       # 🔧 Minimum samples per condition for DE
}

# ============================================================================
# DIFFERENTIAL EXPRESSION PARAMETERS
# ============================================================================

DE_PARAMS = {
    'min_count': 5,         # 🔧 Minimum count threshold for gene filtering
    'min_samples_expr': 2,  # 🔧 Minimum samples expressing gene
    'fc_threshold': 0.5,    # 🔧 Minimum log2 fold change for significance
    'pval_threshold': 0.05, # 🔧 P-value threshold (before FDR correction)
    'fdr_threshold': 0.05,  # 🔧 FDR threshold (after multiple testing correction)
}

# ============================================================================
# VISUALIZATION PARAMETERS
# ============================================================================

VIZ_PARAMS = {
    'top_n_genes': 50,      # 🔧 Number of top genes for heatmaps
    'volcano_fc_threshold': 0.5,  # 🔧 Log2FC threshold for volcano plot coloring
    'volcano_pval_threshold': 0.05,  # 🔧 P-value threshold for volcano plot
}

print("\nAnalysis configuration:")
print(f"  Cell types to analyze: {len(CELL_TYPES_TO_ANALYZE)}")
print(f"\n  Pseudobulk:")
print(f"    Min cells per sample: {PSEUDOBULK_PARAMS['min_cells']}")
print(f"    Min samples per condition: {PSEUDOBULK_PARAMS['min_samples']}")
print(f"\n  DE filtering:")
print(f"    Min count: {DE_PARAMS['min_count']}")
print(f"    Min samples expressing: {DE_PARAMS['min_samples_expr']}")
print(f"    Significance: |log2FC| > {DE_PARAMS['fc_threshold']}, FDR < {DE_PARAMS['fdr_threshold']}")

# %% [markdown]
# ## 4. Prepare Data for Analysis
#
# Create condition column and prepare for pseudobulk aggregation.

# %%
# Create condition column (Genotype_Stimulation)
adata = create_condition_column(adata)

print("Condition distribution:")
print(adata.obs['condition'].value_counts().sort_index())

# Check sample metadata
print("\nSample metadata:")
sample_info = adata.obs.groupby('orig.ident').agg({
    'Genotype': 'first',
    'Stimulation': 'first',
    'Sex': 'first',
    'condition': 'first'
}).reset_index()
print(sample_info.to_string(index=False))

# %% [markdown]
# ## 5. Create Pseudobulk Samples
#
# Aggregate cells per donor × cell type to create pseudobulk samples. This reduces single-cell noise and enables standard bulk RNA-seq statistical methods.

# %%
# Create pseudobulk samples
print("\n" + "="*60)
print("CREATING PSEUDOBULK SAMPLES")
print("="*60)

pb_df, sample_info_df = create_pseudobulk(
    adata,
    groupby=['orig.ident', 'celltype'],
    min_cells=PSEUDOBULK_PARAMS['min_cells']
)

print(f"\n✓ Created {pb_df.shape[1]} pseudobulk samples")
print(f"  Genes: {pb_df.shape[0]:,}")
print(f"\nSample info summary:")
print(sample_info_df.groupby('celltype').agg({
    'sample_id': 'count',
    'n_cells': ['mean', 'min', 'max']
}).round(1))

# %%
# Check condition distribution per cell type
print("\nCondition distribution per cell type:")
for ct in CELL_TYPES_TO_ANALYZE:
    ct_samples = sample_info_df[sample_info_df['celltype'] == ct]
    if len(ct_samples) > 0:
        print(f"\n{ct}:")
        cond_counts = ct_samples['condition'].value_counts().sort_index()
        print(cond_counts)
        
        # Check if we have enough samples for DE
        n_conditions = ct_samples['condition'].nunique()
        if n_conditions < 2:
            print(f"  ⚠️  Only {n_conditions} condition(s) - cannot run DE analysis")
        elif len(ct_samples) < PSEUDOBULK_PARAMS['min_samples'] * 2:
            print(f"  ⚠️  Only {len(ct_samples)} samples - may have low power")

# %%
## 6. Differential Expression Analysis

# Run differential expression analysis for each cell type. We'll test multiple contrasts:
# - E3_GENUS vs E3_Ctrl (effect of stimulation in E3)
# - E4_GENUS vs E4_Ctrl (effect of stimulation in E4)
# - E4_Ctrl vs E3_Ctrl (genotype effect in control)
# - E4_GENUS vs E3_GENUS (genotype effect in stimulated)
from utils.differential_expression import run_de_for_celltype


# %%
# Run DE analysis for each cell type
all_de_results = []

for cell_type in CELL_TYPES_TO_ANALYZE:
    de_result = run_de_for_celltype(
        pb_df,
        sample_info_df,
        cell_type,
        DE_PARAMS
    )
    
    if de_result is not None:
        all_de_results.append(de_result)

if all_de_results:
    combined_results = pd.concat(all_de_results, ignore_index=True)
    
    # Save results
    output_file = OUTPUT_DIR / 'differential_expression_results.csv'
    combined_results.to_csv(output_file, index=False)
    
    print(f"\n{'='*60}")
    print("DE ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"✓ Total results: {len(combined_results):,} gene-contrast pairs")
    print(f"✓ Significant (FDR < {DE_PARAMS['fdr_threshold']}, |logFC| > {DE_PARAMS['fc_threshold']}): {combined_results['significant'].sum():,}")
    
    # Summary by cell type and contrast
    summary = combined_results.groupby(['cell_type', 'contrast']).agg({
        'significant': 'sum',
        'upregulated': 'sum',
        'downregulated': 'sum',
    }).round(0).astype(int)
    
    print(f"\nSignificant genes by cell type and contrast:")
    print(summary)
    
    print(f"\n✓ Saved results to: {output_file}")
else:
    print("\n⚠️  No DE results generated. Check cell type availability and sample sizes.")
    combined_results = None

# %% [markdown]
# ## 7. Visualization: Heatmaps
#
# Create heatmaps showing upregulated and downregulated genes for each cell type and contrast.

# %%
# Import heatmap plotting function from utils
from utils.differential_expression import plot_de_heatmap

# %%
# Plot heatmaps for each cell type and contrast
if combined_results is not None:
    print("\n" + "="*60)
    print("GENERATING HEATMAPS")
    print("="*60)
    
    for cell_type in CELL_TYPES_TO_ANALYZE:
        ct_results = combined_results[combined_results['cell_type'] == cell_type]
        contrasts = ct_results['contrast'].unique()
        
        for contrast in contrasts:
            # Also plot heatmap of NUMBER of upregulated and downregulated genes, before disaggregating by individual genes
            save_path = OUTPUT_DIR / f'heatmap_{cell_type}_{contrast}.png'
            # plot_de_heatmap(
            #     pb_df,
            #     sample_info_df,
            #     combined_results,
            #     cell_type,
            #     contrast,
            #     top_n=VIZ_PARAMS['top_n_genes'],
            #     save_path=save_path
            # )

# %% [markdown]
# ## 8. Visualization: Volcano Plots
#
# Create volcano plots showing differential expression patterns for each cell type and contrast.

# %% [markdown]
# ## 8. Visualization: Volcano Plots
#
# Create volcano plots showing differential expression patterns for each cell type and contrast.

# %%
# Import volcano plot function from utils
from utils.differential_expression import plot_volcano

# %%
# Plot volcano plots for each cell type and contrast
if combined_results is not None:
    print("\n" + "="*60)
    print("GENERATING VOLCANO PLOTS")
    print("="*60)
    
    for cell_type in CELL_TYPES_TO_ANALYZE:
        ct_results = combined_results[combined_results['cell_type'] == cell_type]
        contrasts = ct_results['contrast'].unique()
        
        for contrast in contrasts:
            save_path = OUTPUT_DIR / f'volcano_{cell_type}_{contrast}.png'
            # plot_volcano(
            #     combined_results,
            #     cell_type,
            #     contrast,
            #     fc_threshold=VIZ_PARAMS['volcano_fc_threshold'],
            #     pval_threshold=VIZ_PARAMS['volcano_pval_threshold'],
            #     save_path=save_path
            # )

# %% [markdown]
# ## 9. Summary Statistics
#
# Compile and export summary statistics for all differential expression results.

# %%
# Create summary statistics
if combined_results is not None:
    print("\n" + "="*60)
    print("DE ANALYSIS SUMMARY STATISTICS")
    print("="*60)
    
    # Summary by cell type and contrast
    summary = combined_results.groupby(['cell_type', 'contrast']).agg({
        'significant': 'sum',
        'upregulated': 'sum',
        'downregulated': 'sum',
        'gene': 'count',  # Total genes tested
    }).round(0).astype(int)
    summary.columns = ['Significant', 'Upregulated', 'Downregulated', 'Total_genes']
    
    print("\nSignificant genes by cell type and contrast:")
    print(summary)
    
    # Save summary
    summary_file = OUTPUT_DIR / 'de_summary_statistics.csv'
    summary.to_csv(summary_file)
    print(f"\n✓ Saved summary to: {summary_file}")
    
    # Top genes overall
    print("\n" + "="*60)
    print("TOP DE GENES (across all cell types)")
    print("="*60)
    
    # Get top upregulated and downregulated genes
    top_up = combined_results.nlargest(20, 'logFC')
    top_down = combined_results.nsmallest(20, 'logFC')
    
    print("\nTop 20 upregulated genes:")
    print(top_up[['gene', 'cell_type', 'contrast', 'logFC', 'adj.P.Val']].to_string(index=False))
    
    print("\nTop 20 downregulated genes:")
    print(top_down[['gene', 'cell_type', 'contrast', 'logFC', 'adj.P.Val']].to_string(index=False))
    
    # Save top genes
    top_genes_file = OUTPUT_DIR / 'top_de_genes.csv'
    top_genes = pd.concat([
        top_up[['gene', 'cell_type', 'contrast', 'logFC', 'adj.P.Val']].assign(direction='up'),
        top_down[['gene', 'cell_type', 'contrast', 'logFC', 'adj.P.Val']].assign(direction='down')
    ])
    top_genes.to_csv(top_genes_file, index=False)
    print(f"\n✓ Saved top genes to: {top_genes_file}")

# %% [markdown]
# ## 9.5. Summary Heatmaps: Upregulated and Downregulated Genes
#
# Create heatmaps showing the number of upregulated and downregulated genes across all cell types and comparisons.

# %%
# Create summary heatmaps for upregulated and downregulated genes
if combined_results is not None:
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Create summary if not already created
    if 'summary' not in locals():
        summary = combined_results.groupby(['cell_type', 'contrast']).agg({
            'significant': 'sum',
            'upregulated': 'sum',
            'downregulated': 'sum',
            'gene': 'count',
        }).round(0).astype(int)
        summary.columns = ['Significant', 'Upregulated', 'Downregulated', 'Total_genes']
    
    # Pivot to create heatmap data: cell types as rows, contrasts as columns
    # Get unique cell types and contrasts in a consistent order
    cell_types = sorted(combined_results['cell_type'].unique())
    contrasts = ['E3_GENUS_vs_Ctrl', 'E4_GENUS_vs_Ctrl', 'E4_vs_E3_Ctrl', 'E4_vs_E3_GENUS']
    
    # Create matrices for upregulated and downregulated
    up_matrix = summary['Upregulated'].unstack(fill_value=0)
    down_matrix = summary['Downregulated'].unstack(fill_value=0)
    
    # Reindex to ensure all cell types and contrasts are present
    up_matrix = up_matrix.reindex(index=cell_types, columns=contrasts, fill_value=0)
    down_matrix = down_matrix.reindex(index=cell_types, columns=contrasts, fill_value=0)
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, len(cell_types) * 0.4)))
    
    # Plot upregulated heatmap
    sns.heatmap(
        up_matrix,
        annot=True,
        fmt='d',
        cmap='Reds',
        cbar_kws={'label': 'Number of upregulated genes'},
        ax=axes[0],
        linewidths=0.5,
        linecolor='gray'
    )
    axes[0].set_title('Upregulated Genes', fontsize=14, fontweight='bold')
    axes[0].set_xlabel('Comparison', fontsize=12)
    axes[0].set_ylabel('Cell Type', fontsize=12)
    axes[0].tick_params(axis='x', rotation=45)
    
    # Plot downregulated heatmap
    sns.heatmap(
        down_matrix,
        annot=True,
        fmt='d',
        cmap='Blues',
        cbar_kws={'label': 'Number of downregulated genes'},
        ax=axes[1],
        linewidths=0.5,
        linecolor='gray'
    )
    axes[1].set_title('Downregulated Genes', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Comparison', fontsize=12)
    axes[1].set_ylabel('Cell Type', fontsize=12)
    axes[1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    # Save figure
    heatmap_summary_path = OUTPUT_DIR / 'de_summary_heatmaps.png'
    plt.savefig(heatmap_summary_path, dpi=300, bbox_inches='tight')
    print(f"✓ Saved summary heatmaps to: {heatmap_summary_path}")
    
    plt.show()
else:
    print("⚠️  No DE results available for heatmap creation.")

# %% [markdown]
# ## 10. Export for DESeq2 (Optional)
#
# For more robust statistical analysis, export pseudobulk data for DESeq2 analysis in R. DESeq2 uses negative binomial models and is generally more powerful than t-tests for RNA-seq data.

# %%
# Export pseudobulk data for DESeq2 analysis
print("\n" + "="*60)
print("EXPORTING FOR DESEQ2")
print("="*60)

# Export counts and sample info for each cell type
for cell_type in CELL_TYPES_TO_ANALYZE:
    ct_samples = sample_info_df[sample_info_df['celltype'] == cell_type]
    
    if len(ct_samples) < 4:
        print(f"⚠️  Skipping {cell_type}: Not enough samples")
        continue
    
    # Get counts for this cell type
    ct_counts = pb_df[ct_samples['group_id']]
    
    # Export counts
    counts_file = OUTPUT_DIR / f'pseudobulk_counts_{cell_type}.csv'
    ct_counts.to_csv(counts_file)
    
    # Export sample info
    sample_file = OUTPUT_DIR / f'sample_info_{cell_type}.csv'
    ct_samples[['group_id', 'sample_id', 'Genotype', 'Stimulation', 
                'Sex', 'condition', 'n_cells']].to_csv(sample_file, index=False)
    
    print(f"✓ Exported {cell_type}: {ct_counts.shape[0]:,} genes × {len(ct_samples)} samples")
    print(f"    Counts: {counts_file}")
    print(f"    Sample info: {sample_file}")

print("\n💡 To run DESeq2 in R:")
print("   1. Load the CSV files in R")
print("   2. Use DESeqDataSetFromMatrix() to create DESeqDataSet")
print("   3. Run DESeq() and results() for each contrast")
print("   4. See r/2_limma_voom_gsea.Rmd for example code")

# %% [markdown]
# ## 11. Analysis Complete! 🎉
#
# ### Output Files
#
# **Main results:**
# - `outputs/differential_expression_results/differential_expression_results.csv` - Full DE results
# - `outputs/differential_expression_results/de_summary_statistics.csv` - Summary statistics
# - `outputs/differential_expression_results/top_de_genes.csv` - Top DE genes
#
# **Pseudobulk data (for DESeq2):**
# - `outputs/differential_expression_results/pseudobulk_counts_*.csv` - Count matrices per cell type
# - `outputs/differential_expression_results/sample_info_*.csv` - Sample metadata per cell type
#
# **Plots:**
# - `plots/differential_expression/heatmap_*.png` - Heatmaps for each cell type and contrast
# - `plots/differential_expression/volcano_*.png` - Volcano plots for each cell type and contrast
#
# ### Next Steps
#
# 1. **Review significant genes** - Check top DE genes for biological relevance
# 2. **Pathway enrichment** - Run GSEA or GO enrichment on significant genes
# 3. **Validate with literature** - Check if DE genes match known pathways
# 4. **DESeq2 analysis** - For more robust statistics, use exported data in R
# 5. **Cell type-specific effects** - Compare DE patterns across cell types
#
# ### Notes
#
# - Current analysis uses t-tests on log-CPM data
# - For publication-quality results, consider using DESeq2 (negative binomial model)
# - Multiple testing correction (FDR) has been applied
# - Focus on genes with both significant p-values AND meaningful fold changes
