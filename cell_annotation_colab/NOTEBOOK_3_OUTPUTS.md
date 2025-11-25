# Notebook 3 Outputs: What's Available for Notebook 4

## Primary Artifact: `annotated_data.h5ad`

The main output from Notebook 3 is `outputs/annotated_data.h5ad`, which is a complete AnnData object containing:

### 📊 **Cell Observations** (`adata.obs`)
Per-cell metadata including:
- **`leiden`** - Cluster assignments (from Notebook 2)
- **`celltype`** - Cell type annotations (e.g., 'Excit', 'Inhib', 'Astro', 'ExN_L6b', 'InN_PVALB', etc.)
- **`celltype_cluster`** - Cluster-level aggregated cell type labels
- **`annotation_confidence`** - 'high' or 'low' confidence in cell type assignment
- **`orig.ident`** - Sample ID (e.g., 'D25-2675')
- **`Genotype`** - E3 or E4
- **`Stimulation`** - Ctrl or GENUS
- **`Sex`** - Male or Female
- **`n_genes_by_counts`** - Number of genes detected per cell
- **`total_counts`** - Total UMI counts per cell
- **`percent_mt`** - Percentage of mitochondrial reads
- **`doublet_score`** - Doublet detection score
- **`leiden_excit`** / **`leiden_inhib`** - Sub-cluster assignments for neuronal subtypes (if re-clustering was performed)

### 🧬 **Gene Expression Data**
- **`adata.X`** - Normalized, log-transformed expression matrix (processed data)
- **`adata.raw`** - Raw count matrix (if preserved from Notebook 1) - **CRITICAL for differential expression**

### 📐 **Dimensionality Reduction** (`adata.obsm`)
- **`X_pca`** - PCA coordinates
- **`X_umap`** - UMAP coordinates for visualization

### 🔬 **Gene Information** (`adata.var`)
- Gene names and filtering information
- Highly variable gene flags

### 📝 **Metadata** (`adata.uns`)
- **`pipeline_params`** - All parameters used in notebooks 1-3
  - QC thresholds
  - Clustering parameters  
  - Annotation parameters
  - Marker gene panels used
- **`rank_genes_groups`** - Marker genes per cluster (from Notebook 2)
- Plotting parameters and color schemes

---

## Secondary Outputs (CSV Files)

### 1. **`outputs/cell_metadata.csv`**
A flat CSV file containing key per-cell metadata:
- leiden (cluster)
- celltype
- orig.ident (sample)
- Genotype
- Sex
- n_genes_by_counts
- total_counts
- percent_mt
- doublet_score

**Use case:** Quick exploration in Excel/R, or integration with non-Python tools

### 2. **`outputs/celltype_counts.csv`**
Cell type distribution across the entire dataset
```
celltype,count
Astro,1504
Endo,247
ExN_L2-3,74
...
```

**Use case:** Quick summary statistics

### 3. **`outputs/analysis_summary.csv`**
High-level dataset statistics:
- Total cells
- Total genes
- Number of clusters
- Number of cell types
- Median genes/cell
- Median UMIs/cell
- Median MT%
- Number of samples

**Use case:** Quality control reporting, manuscript methods

### 4. **`outputs/cluster_purity.csv`** (Generated in Notebook 4)
Shows which clusters are "pure" (dominated by one cell type)
- cluster
- dominant_celltype
- purity (proportion)
- n_cells

**Use case:** Quality assessment of cell type annotations

---

## Visualization Outputs

### Plots from Notebook 3:
- `plots/notebook3/marker_genes_dotplot.png` - Marker gene expression across clusters
- `plots/notebook3/celltype_umap.png` - Cell type annotations on UMAP (3-panel view)
- `plots/notebook3/composition_heatmap.png` - Cell type composition per cluster
- `plots/notebook3/celltype_distribution.png` - Cell types per sample (stacked bar)
- `plots/notebook3/reclustering/umap_excit.png` - Excitatory neuron sub-clustering
- `plots/notebook3/reclustering/umap_inhib.png` - Inhibitory neuron sub-clustering

---

## What Notebook 4 Uses from `annotated_data.h5ad`

### **Essential Data:**
1. **`adata.obs['celltype']`** - Cell type labels for grouping
2. **`adata.obs['leiden']`** - Cluster labels for purity filtering
3. **`adata.obs['orig.ident']`** - Sample IDs for pseudobulk aggregation
4. **`adata.obs['Genotype']`** and **`adata.obs['Stimulation']`** - Experimental conditions for DE contrasts
5. **`adata.raw.X`** or **`adata.X`** - Expression counts for differential expression
   - **Ideally `adata.raw.X`** (raw counts) for proper statistical testing
   - Falls back to `adata.X` (normalized) if raw is not available

### **Why Load the H5AD Instead of CSVs?**
The h5ad file contains:
- ✅ **Gene expression data** (millions of values) - too large for CSV
- ✅ **All cell metadata** in one place
- ✅ **Gene names** aligned with expression matrix
- ✅ **Pipeline provenance** (all parameters used)
- ✅ **Efficient storage** (compressed, indexed)

The CSVs only contain **metadata summaries** - you can't do differential expression without the expression matrix!

---

## Key Design Decision in Notebook 4

**Cluster Purity Filtering (≥50%):**
- Notebook 4 identifies "pure" clusters where a single cell type dominates
- Only cells from pure clusters are used for DE analysis
- This prevents:
  - Mixed cell type signals confounding results
  - Doublets contaminating cell type-specific signatures
  - Low-confidence annotations driving false positives

**Example:**
- Cluster 10: 84.8% Astro → **INCLUDED** (pure)
- Cluster 2: 28.2% OPC, 25% Excit, 20% Oligo → **EXCLUDED** (mixed)

---

## Alternative: What if You Only Had CSVs?

If you only had the CSV files, you would:
❌ **Cannot do differential expression** - no expression data
✅ Can do cell type proportion analysis
✅ Can do quality control checks
✅ Can visualize metadata distributions

**Bottom line:** The h5ad file is essential for any analysis requiring gene expression (DE, pathway analysis, etc.)

---

## Summary Table

| Artifact | Contains | Size | Format | Required for Notebook 4? |
|----------|----------|------|--------|--------------------------|
| `annotated_data.h5ad` | Expression + metadata + embeddings | ~500 MB | Binary (HDF5) | ✅ **YES - Essential** |
| `cell_metadata.csv` | Key cell-level metadata | ~2 MB | Text | ❌ No (redundant with h5ad) |
| `celltype_counts.csv` | Cell type summary | ~1 KB | Text | ❌ No |
| `analysis_summary.csv` | Dataset statistics | <1 KB | Text | ❌ No |
| `cluster_purity.csv` | Cluster composition | ~1 KB | Text | ❌ No (computed in NB4) |

**Note:** Notebook 4 generates its own additional outputs for differential expression results.


