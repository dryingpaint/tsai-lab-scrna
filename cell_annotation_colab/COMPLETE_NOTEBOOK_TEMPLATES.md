# Complete Notebook Templates with Tuning Cells

This document contains the complete structure for all three notebooks with integrated tuning cells.

---

## Notebook 1: `1_setup_qc_filtering.ipynb`

### Cell Structure:

1. **Header** (markdown)
2. **Table of Contents** (markdown)
3. **Setup & Installation** (code)
4. **Mount Google Drive** (code - optional)
5. **Parameter Configuration** (markdown + code cells)
   - Presets
   - Data paths
   - QC filters
   - Doublet params
   - Summary display
6. **Stage 1: Data Loading** (markdown)
   - Custom H5 loader function (code)
   - Load and merge function (code)
   - Metadata function (code)
   - Run loading (code)
7. **Stage 2: QC Metrics** (markdown)
   - Calculate QC metrics (code)
   - Plot QC metrics (code)
   - **🎛️ TUNING CELL** (markdown - collapsible decisions)
   - **📊 Analysis Cell** (code - automated assessment)
8. **Stage 3: Doublet Detection** (markdown)
   - Detect doublets function (code)
   - Prep cells for doublet detection (code)
   - Run doublet detection (code)
   - **🎛️ TUNING CELL** (markdown - per-sample assessment)
   - **📊 Analysis Cell** (code - doublet summary stats)
9. **Stage 4: Filtering** (markdown)
   - Filter function (code)
   - Run filtering (code)
   - Visualization (code)
   - **🎛️ TUNING CELL** (markdown - retention analysis)
   - **📊 Analysis Cell** (code - filtering impact)
10. **Save Output** (code)
11. **Summary & Next Steps** (markdown)

### Key Code Sections:

#### Custom H5 Loader (exact from corrected notebook):
```python
def load_cellbender_h5(file_path):
    \"\"\"Load CellBender processed h5 file\"\"\"
    with h5py.File(file_path, 'r') as f:
        # [Full implementation from corrected notebook]
    return adata
```

#### Doublet Detection (per-sample, with threshold capping):
```python
def detect_doublets_improved(adata, expected_doublet_rate=0.10, manual_threshold=0.35,
                            plot_histograms=True, save_dir=None):
    \"\"\"Detect doublets per-sample\"\"\"
    # [Full implementation with per-sample loop]
    return adata
```

#### Filtering (correct order - doublets last):
```python
def filter_cells_and_genes(adata, ...):
    \"\"\"Filter in correct order: QC filters first, doublets LAST\"\"\"
    # 1. min_genes
    # 2. min_cells
    # 3. max_genes
    # 4. max_mt_pct
    # 5. count filters
    # 6. ribo%
    # 7. doublets (LAST!)
    return adata
```

### Tuning Cells After Stage 2:

```markdown
### 🎛️ Parameter Tuning: QC Results

**Review the plots above and use this guide:**

#### Genes per Cell

<details>
<summary>📊 Two distinct populations (low ~500 and high ~3000)</summary>

**Action:**
\```python
# In Parameter Configuration section:
CELL_FILTERS['min_genes'] = 300
\```
**Then:** Re-run from Stage 2
</details>

<details>
<summary>📊 Long tail to >10,000 genes</summary>

**Action:**
\```python
CELL_FILTERS['max_genes'] = 6000
\```
**Then:** Re-run from Stage 2
</details>

[... more scenarios ...]

#### Mitochondrial %

<details>
<summary>📊 Main population 2-5% with tail to 20%+</summary>

**Action:**
\```python
CELL_FILTERS['max_mt_pct'] = 10  # or 5 for neurons
\```
</details>

[... more scenarios ...]
```

### Analysis Cell After Stage 2:

```python
# QC Metrics Assessment
print("\\n" + "="*60)
print("QC METRICS ASSESSMENT")
print("="*60)

median_genes = adata.obs['n_genes_by_counts'].median()
median_mt = adata.obs['percent_mt'].median()

print(f"Median genes per cell: {median_genes:.0f}")
if median_genes < 2000:
    print("⚠️  LOW: Consider checking sequencing depth")
elif median_genes > 5000:
    print("✅ EXCELLENT: High quality data")
else:
    print("✅ GOOD: Normal range")

print(f"\\nMedian MT%: {median_mt:.2f}%")
if median_mt > 10:
    print("⚠️  HIGH: Consider max_mt_pct = 8")
elif median_mt < 5:
    print("✅ EXCELLENT: Low stress")
else:
    print("✅ GOOD: Acceptable range")

print("\\n💡 NEXT STEPS:")
if median_genes < 2000 or median_mt > 10:
    print("   • Adjust parameters and re-run Stage 2")
else:
    print("   • Proceed to Stage 3: Doublet Detection")
```

### Final Cell - Save Output:

```python
# Save QC-filtered data
output_file = 'outputs/qc_filtered_data.h5ad'

# Store parameters used in uns for reproducibility
adata.uns['pipeline_params'] = {
    'notebook': '1_setup_qc_filtering',
    'cell_filters': CELL_FILTERS,
    'doublet_params': DOUBLET_PARAMS,
    'gene_patterns': GENE_PATTERNS,
    'gene_filters': GENE_FILTERS,
}

adata.write(output_file)

print("\\n" + "="*60)
print("NOTEBOOK 1 COMPLETE")
print("="*60)
print(f"✓ Saved: {output_file}")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Size: {os.path.getsize(output_file) / 1e6:.1f} MB")
print("\\n➡️  NEXT: Open 2_clustering_markers.ipynb")
```

---

## Notebook 2: `2_clustering_markers.ipynb`

### Cell Structure:

1. **Header** (markdown)
2. **Load Previous Data** (code + validation)
3. **Parameter Configuration** (markdown + code)
   - PCA parameters
   - Clustering parameters
   - Marker gene parameters
4. **Stage 5: Normalization** (markdown + code)
5. **Stage 6: PCA & UMAP** (markdown + code + tuning + analysis)
6. **Stage 7: Clustering** (markdown + code + tuning + analysis)
7. **Stage 7: Marker Genes** (markdown + code + tuning + analysis)
8. **Save Output** (code)
9. **Summary** (markdown)

### Key Sections:

#### Load and Validate:

```python
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# Load QC-filtered data from Notebook 1
print("Loading data from Notebook 1...")
adata = sc.read_h5ad('outputs/qc_filtered_data.h5ad')

# Validation
print("\\n" + "="*60)
print("DATA VALIDATION")
print("="*60)

checks = {
    'Has QC metrics': 'percent_mt' in adata.obs.columns,
    'Has doublet scores': 'doublet_score' in adata.obs.columns,
    'Has sample info': 'orig.ident' in adata.obs.columns,
    'Has gene annotations': 'mt' in adata.var.columns,
}

all_passed = True
for check, passed in checks.items():
    status = "✓" if passed else "✗"
    print(f"  {status} {check}")
    if not passed:
        all_passed = False

if not all_passed:
    raise ValueError("Data validation failed! Did you run Notebook 1?")

print("\\n✓ Validation passed")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

# Show parameters used in Notebook 1
if 'pipeline_params' in adata.uns:
    print("\\nParameters from Notebook 1:")
    for key, value in adata.uns['pipeline_params'].items():
        if key != 'notebook':
            print(f"  {key}: {value}")
```

#### Stage 5: Normalization

```python
def normalize_and_scale(adata):
    \"\"\"Standard normalization workflow\"\"\"
    print("Normalizing and scaling...")

    # Save raw
    adata.raw = adata

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVGs
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    n_hvg = adata.var['highly_variable'].sum()
    print(f"  ✓ {n_hvg:,} HVGs identified")

    # Plot
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(PLOTS_DIR / 'highly_variable_genes.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Keep HVGs
    adata = adata[:, adata.var.highly_variable]

    # Scale
    sc.pp.scale(adata, max_value=10)

    return adata
```

#### Tuning Cell After PCA:

```markdown
### 🎛️ Parameter Tuning: PCA Results

**Check the elbow plot above:**

<details>
<summary>📊 Clear elbow at PC 15-20</summary>

**Action:** Current N_PCS = 15 is good, proceed
</details>

<details>
<summary>📊 Elbow at PC 30-40</summary>

**Action:**
\```python
N_PCS = 35
\```
**Then:** Re-run from Stage 6
</details>

<details>
<summary>📊 Elbow at PC 8-10</summary>

**Action:**
\```python
N_PCS = 10
\```
**Check:** May indicate homogeneous population or over-filtering
</details>
```

#### Tuning Cell After Clustering:

```markdown
### 🎛️ Parameter Tuning: Clustering Results

**How many clusters were identified?**

<details>
<summary>📊 <5 clusters</summary>

**Diagnosis:** Under-clustering

**Action:**
\```python
CLUSTERING_PARAMS['resolution'] = 1.0  # Increase
\```
**Then:** Re-run clustering only
</details>

<details>
<summary>📊 5-20 clusters</summary>

**Diagnosis:** ✅ Likely good

**Action:** Validate with marker genes
</details>

<details>
<summary>📊 >30 clusters</summary>

**Diagnosis:** Over-clustering

**Action:**
\```python
CLUSTERING_PARAMS['resolution'] = 0.4  # Decrease
\```
</details>
```

#### Save Output:

```python
# Save clustered data
output_file = 'outputs/clustered_data.h5ad'

# Add parameters to uns
adata.uns['pipeline_params']['notebook'] = '2_clustering_markers'
adata.uns['pipeline_params']['n_pcs'] = N_PCS
adata.uns['pipeline_params']['n_neighbors'] = N_NEIGHBORS
adata.uns['pipeline_params']['clustering'] = CLUSTERING_PARAMS

adata.write(output_file)

print("\\n" + "="*60)
print("NOTEBOOK 2 COMPLETE")
print("="*60)
print(f"✓ Saved: {output_file}")
print(f"  Cells: {adata.n_obs:,}")
print(f"  HVGs: {adata.n_vars:,}")
print(f"  Clusters: {adata.obs['leiden'].nunique()}")
print(f"  Embeddings: {list(adata.obsm.keys())}")
print("\\n➡️  NEXT: Open 3_annotation_export.ipynb")
```

---

## Notebook 3: `3_annotation_export.ipynb`

### Cell Structure:

1. **Header** (markdown)
2. **Load Previous Data** (code + validation)
3. **Parameter Configuration** (markdown + code)
   - Annotation parameters
   - Marker gene dictionary
4. **Stage 8: Cell Type Annotation** (markdown + code + tuning + analysis)
5. **Stage 9: Reclustering** (markdown + code)
6. **Stage 10: Final Export** (code)
7. **Summary & Next Steps** (markdown)

### Key Sections:

#### Load and Validate:

```python
# Load clustered data from Notebook 2
print("Loading data from Notebook 2...")
adata = sc.read_h5ad('outputs/clustered_data.h5ad')

# Validation
print("\\n" + "="*60)
print("DATA VALIDATION")
print("="*60)

checks = {
    'Has UMAP': 'X_umap' in adata.obsm,
    'Has clusters': 'leiden' in adata.obs.columns,
    'Has marker genes': 'rank_genes_groups' in adata.uns,
    'Has raw data': adata.raw is not None,
}

for check, passed in checks.items():
    status = "✓" if passed else "✗"
    print(f"  {status} {check}")
    if not passed:
        raise ValueError(f"Missing: {check}. Did you run Notebook 2?")

print("\\n✓ Validation passed")
print(f"  Clusters: {adata.obs['leiden'].nunique()}")
```

#### Annotation with Scoring:

```python
# Score cells for each cell type
print("Scoring cells for major cell types...")

for cell_type, genes in MARKER_GENES.items():
    available_genes = [g for g in genes if g in adata.raw.var_names]

    if available_genes:
        sc.tl.score_genes(adata, available_genes,
                         score_name=f'{cell_type}_score',
                         use_raw=True)
        print(f"  ✓ {cell_type}: {len(available_genes)}/{len(genes)} markers")

# Assign cell types
score_cols = [col for col in adata.obs.columns if col.endswith('_score')]
scores = adata.obs[score_cols]

# Apply margin threshold
scores_sorted = np.sort(scores.values, axis=1)
max_scores = scores_sorted[:, -1]
second_scores = scores_sorted[:, -2]
confident = (max_scores - second_scores) > ANNOTATION_PARAMS['margin']

# Assign
adata.obs['celltype'] = scores.idxmax(axis=1).str.replace('_score', '')
adata.obs.loc[~confident, 'celltype'] = 'Unlabeled'
adata.obs['celltype_confidence'] = max_scores - second_scores

print(f"\\n✓ {confident.sum():,} / {len(confident):,} confidently labeled ({confident.sum()/len(confident)*100:.1f}%)")
```

#### Tuning Cell After Annotation:

```markdown
### 🎛️ Parameter Tuning: Annotation Results

**What percentage are "Unlabeled"?**

<details>
<summary>📊 <10% unlabeled</summary>

**Diagnosis:** ✅ Excellent coverage

**Action:** Proceed to export
</details>

<details>
<summary>📊 >20% unlabeled</summary>

**Diagnosis:** Poor coverage

**Actions:**
\```python
# Lower confidence threshold:
ANNOTATION_PARAMS['margin'] = 0.03

# OR try cluster-level:
ANNOTATION_PARAMS['label_mode'] = 'cluster'
\```
**Then:** Re-run annotation
</details>

**Do cell types cluster spatially on UMAP?**

<details>
<summary>📊 Yes - spatially coherent</summary>

**Diagnosis:** ✅ Biologically meaningful

**Action:** Proceed with confidence
</details>

<details>
<summary>📊 No - scattered</summary>

**Diagnosis:** Check marker specificity

**Action:** Review marker genes or try different annotation approach
</details>
```

#### Final Export:

```python
# Final export with all metadata
output_file = 'outputs/annotated_data.h5ad'

# Add final parameters
adata.uns['pipeline_params']['notebook'] = '3_annotation_export'
adata.uns['pipeline_params']['annotation'] = ANNOTATION_PARAMS
adata.uns['pipeline_params']['marker_genes_used'] = list(MARKER_GENES.keys())

adata.write(output_file)

# Also export metadata CSV
metadata_cols = ['leiden', 'celltype', 'celltype_confidence',
                'orig.ident', 'Genotype', 'Sex', 'Stimulation',
                'n_genes_by_counts', 'total_counts', 'percent_mt', 'doublet_score']
adata.obs[metadata_cols].to_csv('outputs/cell_metadata.csv')

# Summary stats
summary = pd.DataFrame({
    'Metric': ['Total cells', 'Total genes', 'Clusters', 'Cell types',
               'Median genes/cell', 'Median UMIs/cell'],
    'Value': [
        f"{adata.n_obs:,}",
        f"{adata.n_vars:,}",
        adata.obs['leiden'].nunique(),
        adata.obs['celltype'].nunique(),
        f"{adata.obs['n_genes_by_counts'].median():.0f}",
        f"{adata.obs['total_counts'].median():.0f}",
    ]
})
summary.to_csv('outputs/analysis_summary.csv', index=False)

print("\\n" + "="*60)
print("PIPELINE COMPLETE!")
print("="*60)
print(f"✓ Annotated data: {output_file}")
print(f"✓ Metadata CSV: outputs/cell_metadata.csv")
print(f"✓ Summary: outputs/analysis_summary.csv")
print(f"\\n  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  {adata.obs['leiden'].nunique()} clusters")
print(f"  {adata.obs['celltype'].nunique()} cell types")
print("\\n📊 Cell type distribution:")
print(adata.obs['celltype'].value_counts())
print("\\n🎉 Ready for downstream analysis!")
```

---

## Summary

Each notebook:
- ✅ Has clear input/output files
- ✅ Includes data validation
- ✅ Has parameter configuration sections
- ✅ Integrates tuning cells after each stage
- ✅ Includes automated assessment code cells
- ✅ Saves parameters used for reproducibility
- ✅ Has clear next steps

The three notebooks work as a pipeline but can also be run independently with the right input files.
