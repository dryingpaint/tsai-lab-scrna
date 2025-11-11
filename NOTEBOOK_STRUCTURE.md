# Multi-Notebook Workflow Structure

## Strategic Breakpoints

### **Notebook 1: Setup, QC & Filtering** (Stages 1-4)
**Input:** Raw CellBender H5 files
**Output:** `qc_filtered_data.h5ad`

**Rationale:**
- Natural checkpoint after QC filtering
- Users often iterate on QC parameters
- Output is significantly smaller (filtered cells/genes)
- Contains: raw counts, QC metrics, doublet scores

**What's saved:**
```python
adata.write('qc_filtered_data.h5ad')
# Includes:
# - .X: raw counts (sparse)
# - .obs: cell metadata (QC metrics, doublet scores, sample info)
# - .var: gene metadata (mt, ribo flags)
```

**Stages:**
1. Data Loading & Integration
2. QC Metrics Calculation & Visualization
3. Doublet Detection
4. Cell & Gene Filtering

---

### **Notebook 2: Normalization, Clustering & Markers** (Stages 5-7)
**Input:** `qc_filtered_data.h5ad`
**Output:** `clustered_data.h5ad`

**Rationale:**
- Natural checkpoint after clustering and marker identification
- Users iterate on clustering parameters (resolution, PCs, neighbors)
- Output includes embeddings and cluster assignments
- Ready for annotation

**What's saved:**
```python
adata.write('clustered_data.h5ad')
# Includes:
# - .X: scaled data
# - .raw: normalized log-counts (for plotting)
# - .obs: clusters, metadata
# - .var: HVG flags, marker gene statistics
# - .obsm: PCA, UMAP embeddings
# - .uns: marker genes, clustering parameters
# - .layers: could store raw, normalized separately
```

**Stages:**
5. Normalization & Scaling
6. PCA, UMAP & Clustering
7. Marker Gene Analysis

---

### **Notebook 3: Annotation, Reclustering & Export** (Stages 8-9)
**Input:** `clustered_data.h5ad`
**Output:** `annotated_data.h5ad` (final)

**Rationale:**
- Final annotation and refinement
- Users may try different annotation strategies
- Output is publication-ready
- Includes all metadata for downstream analysis

**What's saved:**
```python
adata.write('annotated_data.h5ad')
# Includes everything from Notebook 2 plus:
# - .obs['celltype']: cell type annotations
# - .obs['celltype_score']: annotation confidence scores
# - .uns['annotation_params']: parameters used
# - Subtype reclustering results
```

**Stages:**
8. Cell Type Annotation
9. Reclustering & Export
10. Summary & Next Steps

---

## Data Flow Diagram

```
Raw CellBender H5 files
        ↓
   [Notebook 1]
   Setup & QC
        ↓
qc_filtered_data.h5ad  ← Checkpoint 1
        ↓
   [Notebook 2]
   Clustering
        ↓
clustered_data.h5ad    ← Checkpoint 2
        ↓
   [Notebook 3]
   Annotation
        ↓
annotated_data.h5ad    ← Final Output
```

---

## Loading Strategy for Each Notebook

### **Notebook 1:**
```python
# No loading needed - starts fresh
adata = load_and_merge_cellbender_data(BASE_PATH, SAMPLE_NAMES, CUSTOM_NAME)
```

### **Notebook 2:**
```python
# Load filtered data
adata = sc.read_h5ad('qc_filtered_data.h5ad')

# Verify what we have:
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"QC metrics: {[c for c in adata.obs.columns if 'n_genes' in c or 'percent' in c]}")
```

### **Notebook 3:**
```python
# Load clustered data
adata = sc.read_h5ad('clustered_data.h5ad')

# Verify embeddings and clusters:
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Embeddings: {list(adata.obsm.keys())}")
print(f"Clusters: {adata.obs['leiden'].nunique()} clusters")
```

---

## Parameter Sharing Strategy

### **Option A: Shared Config File**
Create `pipeline_config.py`:
```python
# All parameters defined once
SAMPLE_NAMES = [...]
CELL_FILTERS = {...}
DOUBLET_PARAMS = {...}
# etc.
```

Each notebook imports:
```python
from pipeline_config import *
```

### **Option B: Save Parameters to H5AD** (Recommended)
```python
# Notebook 1 saves parameters used:
adata.uns['pipeline_params'] = {
    'cell_filters': CELL_FILTERS,
    'doublet_params': DOUBLET_PARAMS,
    'gene_patterns': GENE_PATTERNS,
}

# Notebook 2 loads and uses:
params = adata.uns['pipeline_params']
print("Previous QC filters used:")
print(params['cell_filters'])
```

### **Option C: Parameters in Each Notebook** (Simplest for Colab)
- Each notebook has its own parameter section
- Copy from previous notebook if re-running
- Most flexible but requires manual tracking

**Recommendation**: Use Option C for simplicity in Colab environment

---

## Iteration Patterns

### **Scenario 1: Adjust QC filters**
```
User runs Notebook 1 → Sees QC results → Adjusts parameters → Re-runs Notebook 1
        ↓ (satisfied with QC)
User runs Notebook 2 (loads fresh qc_filtered_data.h5ad)
```

### **Scenario 2: Adjust clustering**
```
User has qc_filtered_data.h5ad from previous run
        ↓
User runs Notebook 2 → Adjusts resolution → Re-runs clustering section only
        ↓ (satisfied with clusters)
User runs Notebook 3 (loads clustered_data.h5ad)
```

### **Scenario 3: Try different annotations**
```
User has clustered_data.h5ad from previous run
        ↓
User runs Notebook 3 → Tries different marker sets → Re-runs annotation only
```

---

## File Management

### **Recommended Directory Structure:**
```
project/
├── data/
│   ├── D25-2675/
│   │   └── D25-2675_processed_feature_bc_matrix_filtered.h5
│   ├── D25-2676/
│   └── ...
├── outputs/
│   ├── qc_filtered_data.h5ad
│   ├── clustered_data.h5ad
│   └── annotated_data.h5ad
├── plots/
│   ├── notebook1/
│   │   ├── qc_violin_plots.png
│   │   ├── doublet_score_histograms.png
│   │   └── ...
│   ├── notebook2/
│   │   ├── pca_elbow_plot.png
│   │   ├── umap_embeddings.png
│   │   └── ...
│   └── notebook3/
│       ├── cell_type_umap.png
│       └── ...
└── notebooks/
    ├── 1_setup_qc_filtering.ipynb
    ├── 2_clustering_markers.ipynb
    └── 3_annotation_export.ipynb
```

### **In Colab (Google Drive):**
```
MyDrive/
└── scRNA_project/
    ├── data/           (mount from existing location)
    ├── outputs/
    ├── plots/
    └── notebooks/
```

---

## Cross-Notebook Validation

Each notebook includes validation cells to ensure data integrity:

### **Notebook 2 Start:**
```python
# Validate loaded data
assert 'doublet_score' in adata.obs.columns, "Missing doublet scores - did you run Notebook 1?"
assert 'percent_mt' in adata.obs.columns, "Missing QC metrics"
assert adata.obs['predicted_doublet'].sum() < adata.n_obs, "All cells marked as doublets?"

print("✓ Data validation passed")
```

### **Notebook 3 Start:**
```python
# Validate loaded data
assert 'X_umap' in adata.obsm, "Missing UMAP - did you run Notebook 2?"
assert 'leiden' in adata.obs.columns, "Missing clusters"
assert 'rank_genes_groups' in adata.uns, "Missing marker genes"

print("✓ Data validation passed")
```

---

## Notebook Naming Convention

1. **`1_setup_qc_filtering.ipynb`**
   - Clear numeric prefix for ordering
   - Describes what this notebook does
   - Easy to understand in file browser

2. **`2_clustering_markers.ipynb`**
   - Continues numbering
   - Key stages in name

3. **`3_annotation_export.ipynb`**
   - Final stage
   - Indicates this produces final output

---

## Benefits of This Structure

✅ **Faster iteration**: Don't re-run QC when tuning clustering
✅ **Smaller notebooks**: Each <500KB, loads quickly in Colab
✅ **Clear checkpoints**: Natural save/load points
✅ **Parallel workflows**: Can have multiple annotation strategies from same clustered data
✅ **Better git tracking**: Smaller diffs when notebooks change
✅ **Easier debugging**: Isolated stages
✅ **Flexible**: Skip to any stage with appropriate input file

---

## Implementation Notes

- Each notebook will have identical setup sections (imports, helper functions)
- Consider creating `utils.py` for shared functions
- Each notebook saves parameters used in `.uns` for reproducibility
- Plot directories are separate per notebook to avoid overwrites
- Can run notebooks in sequence or individually (with right input files)
