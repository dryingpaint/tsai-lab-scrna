#!/usr/bin/env python3
"""
Complete the three notebooks with all code and tuning cells.

This script:
1. Splits cell_annotation_pipeline_corrected.ipynb at checkpoints
2. Adds tuning cells after each stage
3. Adds assessment cells
4. Adds remaining stages (6-9) from original pipeline
5. Creates complete, working notebooks
"""

import json
from pathlib import Path

def load_notebook(path):
    """Load a notebook file"""
    with open(path) as f:
        return json.load(f)

def save_notebook(nb, path):
    """Save a notebook file"""
    with open(path, 'w') as f:
        json.dump(nb, f, indent=1)

def create_cell(cell_type, source, metadata=None):
    """Create a notebook cell"""
    cell = {
        "cell_type": cell_type,
        "metadata": metadata or {},
        "source": source if isinstance(source, list) else [source]
    }
    if cell_type == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    return cell

print("="*70)
print("BUILDING COMPLETE NOTEBOOKS")
print("="*70)

# Load source
source_nb = load_notebook("cell_annotation_pipeline_corrected.ipynb")
print(f"\\n✓ Loaded source: {len(source_nb['cells'])} cells")

# Define where to split
# Looking at the analysis above:
# Cell 0-13: Setup & parameters (goes to all notebooks with modifications)
# Cell 14-16: Stage 1 (Notebook 1)
# Cell 17-20: Stage 2 (Notebook 1)
# Cell 21-24: Stage 3 (Notebook 1)
# Cell 25-27: Stage 4 (Notebook 1)
# Cell 28-29: Stage 5 (Notebook 2)
# Cell 30: Continuation note (skip)

# Map out the structure
structure = {
    "setup": (0, 14),      # Cells 0-13: Shared setup
    "stage1": (14, 17),    # Data loading
    "stage2": (17, 21),    # QC metrics
    "stage3": (21, 25),    # Doublet detection
    "stage4": (25, 28),    # Filtering
    "stage5": (28, 30),    # Normalization
}

print("\\nSource structure:")
for key, (start, end) in structure.items():
    print(f"  {key}: cells {start}-{end-1}")

OUTPUT_DIR = Path("cell_annotation_colab")
OUTPUT_DIR.mkdir(exist_ok=True)

# ===================================================================
# NOTEBOOK 1: Setup, QC & Filtering
# ===================================================================
print("\\n" + "="*70)
print("CREATING NOTEBOOK 1")
print("="*70)

nb1_cells = []

# Take setup cells (0-13)
nb1_cells.extend(source_nb['cells'][0:14])
print(f"  Added setup cells: 0-13")

# Take Stage 1-4 cells
nb1_cells.extend(source_nb['cells'][14:28])
print(f"  Added Stage 1-4: cells 14-27")

# Add tuning cell after Stage 2 (QC)
nb1_cells.append(create_cell("markdown", """### 🎛️ Parameter Tuning: QC Results

Review the plots above. What did you observe?

<details>
<summary>📊 Two distinct populations in genes per cell plot</summary>

**Diagnosis:** Low population likely empty droplets/dead cells

**Action:**
```python
# Scroll back to Parameter Configuration and update:
CELL_FILTERS['min_genes'] = 300  # Increase threshold
```
**Then:** Re-run from Stage 2
</details>

<details>
<summary>📊 MT% tail extends to 20%+</summary>

**Diagnosis:** Stressed/dying cells in tail

**Action:**
```python
CELL_FILTERS['max_mt_pct'] = 10  # or 5 for neurons
```
**Then:** Re-run from Stage 2
</details>

<details>
<summary>📊 Metrics look good (tight distributions)</summary>

**Diagnosis:** ✅ High quality data

**Action:** Proceed to Stage 3
</details>

**See `../PARAMETER_TUNING_GUIDE.md` for 35+ scenarios**
"""))

# Add assessment cell after Stage 2
nb1_cells.append(create_cell("code", """# QC Metrics Assessment
print("\\n" + "="*60)
print("QC METRICS ASSESSMENT")
print("="*60)

median_genes = adata.obs['n_genes_by_counts'].median()
median_mt = adata.obs['percent_mt'].median()

print(f"Median genes/cell: {median_genes:.0f}")
if median_genes < 2000:
    print("⚠️  LOW: Check sequencing depth")
elif median_genes > 5000:
    print("✅ EXCELLENT: High quality")
else:
    print("✅ GOOD: Normal range")

print(f"\\nMedian MT%: {median_mt:.2f}%")
if median_mt > 10:
    print("⚠️  HIGH: Consider max_mt_pct = 8")
elif median_mt < 5:
    print("✅ EXCELLENT: Low stress")
else:
    print("✅ GOOD: Acceptable")

print("\\n💡 NEXT STEPS:")
if median_genes < 2000 or median_mt > 10:
    print("   • Adjust parameters and re-run Stage 2")
else:
    print("   • Proceed to Stage 3")
"""))

# Add tuning cell after Stage 3 (Doublets)
nb1_cells.append(create_cell("markdown", """### 🎛️ Parameter Tuning: Doublet Detection

Review the per-sample histograms above.

<details>
<summary>📊 Doublet rates 6-10% per sample</summary>

**Diagnosis:** ✅ Normal for 10x

**Action:** Proceed to Stage 4
</details>

<details>
<summary>📊 One or more samples >15% doublets</summary>

**Diagnosis:** Sample overload or quality issue

**Action:**
```python
# If unexpected:
DOUBLET_PARAMS['manual_threshold'] = 0.30  # More stringent
```
**Then:** Re-run from Stage 3
</details>

<details>
<summary>📊 All samples <3% doublets</summary>

**Diagnosis:** ⚠️ Threshold may be too permissive

**Action:**
```python
DOUBLET_PARAMS['manual_threshold'] = 0.25  # Lower
```
**Then:** Re-run from Stage 3
</details>
"""))

# Add doublet summary cell
nb1_cells.append(create_cell("code", """# Doublet Detection Summary
print("\\n" + "="*60)
print("DOUBLET SUMMARY")
print("="*60)

summary = adata.obs.groupby('orig.ident').agg({
    'predicted_doublet': ['count', 'sum']
})
summary.columns = ['n_cells', 'n_doublets']
summary['pct'] = (summary['n_doublets'] / summary['n_cells'] * 100).round(1)

print("\\nPer-sample:")
print(summary)

overall = (adata.obs['predicted_doublet'].sum() / len(adata.obs)) * 100
print(f"\\nOverall: {overall:.1f}%")

if overall < 3:
    print("⚠️  Very low - consider threshold = 0.25")
elif overall > 20:
    print("⚠️  Very high - investigate quality")
elif 6 <= overall <= 10:
    print("✅ Expected range")
else:
    print("✅ Acceptable")
"""))

# Add tuning after Stage 4 (Filtering)
nb1_cells.append(create_cell("markdown", """### 🎛️ Parameter Tuning: Filtering Results

Check cell retention above.

<details>
<summary>📊 <5,000 cells retained</summary>

**Diagnosis:** May lack power for clustering

**Action:** Review if filters too stringent
```python
CELL_FILTERS['min_genes'] = 150  # Relax
CELL_FILTERS['max_mt_pct'] = 12
```
**Then:** Re-run from Stage 2
</details>

<details>
<summary>📊 5,000-50,000 cells retained</summary>

**Diagnosis:** ✅ Sufficient

**Action:** Proceed to save and Notebook 2
</details>

<details>
<summary>📊 One sample has <500 cells</summary>

**Diagnosis:** Sample-specific issue

**Action:** Investigate that sample's QC metrics
</details>
"""))

# Add filtering assessment
nb1_cells.append(create_cell("code", """# Filtering Assessment
print("\\n" + "="*60)
print("FILTERING ASSESSMENT")
print("="*60)

print(f"Total cells: {adata.n_obs:,}")
print(f"Total genes: {adata.n_vars:,}")

sample_counts = adata.obs['orig.ident'].value_counts()
print(f"\\nCells per sample:")
print(f"  Min: {sample_counts.min():,}")
print(f"  Max: {sample_counts.max():,}")
print(f"  Mean: {sample_counts.mean():.0f}")

if adata.n_obs < 5000:
    print("\\n⚠️  <5k cells - may need to relax filters")
elif adata.n_obs > 50000:
    print("\\n✅ >50k cells - excellent power!")
else:
    print("\\n✅ Sufficient cells for analysis")

if sample_counts.min() < 500:
    low_samples = sample_counts[sample_counts < 500].index.tolist()
    print(f"\\n⚠️  Low-count samples: {low_samples}")
    print("   → Investigate these samples")
"""))

# Add save cell
nb1_cells.append(create_cell("markdown", """## Save Output

Save QC-filtered data for Notebook 2."""))

nb1_cells.append(create_cell("code", """# Save QC-filtered data
output_file = 'outputs/qc_filtered_data.h5ad'
Path('outputs').mkdir(exist_ok=True)

# Store parameters used
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
print(f"  Size: {Path(output_file).stat().st_size / 1e6:.1f} MB")
print("\\n➡️  NEXT: Open 2_clustering_markers.ipynb")
"""))

# Create Notebook 1
nb1 = {
    "cells": nb1_cells,
    "metadata": source_nb['metadata'],
    "nbformat": 4,
    "nbformat_minor": 0
}

output1 = OUTPUT_DIR / "1_setup_qc_filtering.ipynb"
save_notebook(nb1, output1)
print(f"\\n✓ Created {output1.name} with {len(nb1_cells)} cells")

# ===================================================================
# NOTEBOOK 2: Clustering & Markers
# ===================================================================
print("\\n" + "="*70)
print("CREATING NOTEBOOK 2")
print("="*70)

nb2_cells = []

# Header
nb2_cells.append(create_cell("markdown", """# Notebook 2: Clustering & Markers

**Cell Annotation Pipeline - Part 2 of 3**

**Stages:** 5-7
**📥 Input:** `outputs/qc_filtered_data.h5ad`
**📤 Output:** `outputs/clustered_data.h5ad`
**➡️ Next:** `3_annotation_export.ipynb`

---
"""))

# Load and validate
nb2_cells.append(create_cell("markdown", """## Load Data from Notebook 1

Load the QC-filtered data and validate."""))

nb2_cells.append(create_cell("code", """import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load
print("Loading data from Notebook 1...")
adata = sc.read_h5ad('outputs/qc_filtered_data.h5ad')

# Validate
print("\\n" + "="*60)
print("DATA VALIDATION")
print("="*60)

checks = {
    'QC metrics': 'percent_mt' in adata.obs.columns,
    'Doublet scores': 'doublet_score' in adata.obs.columns,
    'Sample info': 'orig.ident' in adata.obs.columns,
}

for check, passed in checks.items():
    print(f"  {'✓' if passed else '✗'} {check}")
    if not passed:
        raise ValueError(f"Missing {check} - run Notebook 1!")

print(f"\\n✓ Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
"""))

# Parameters
nb2_cells.append(create_cell("markdown", """## Parameter Configuration

Set parameters for clustering."""))

nb2_cells.append(create_cell("code", """# Clustering parameters
N_PCS = 15
N_NEIGHBORS = 10
CLUSTERING_PARAMS = {'resolution': 0.8}

# Create plots directory
PLOTS_DIR = Path('plots/notebook2')
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

print("Parameters:")
print(f"  N_PCS: {N_PCS}")
print(f"  N_NEIGHBORS: {N_NEIGHBORS}")
print(f"  Resolution: {CLUSTERING_PARAMS['resolution']}")
"""))

# Stage 5 (from source cells 28-29)
nb2_cells.extend(source_nb['cells'][28:30])

# Stage 6: PCA/UMAP/Clustering
nb2_cells.append(create_cell("markdown", """## Stage 6: PCA, UMAP & Clustering

Reduce dimensionality and identify cell clusters."""))

nb2_cells.append(create_cell("code", """# Run PCA
print("Running PCA...")
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Elbow plot
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.axvline(N_PCS, color='r', linestyle='--', label=f'Selected: {N_PCS}')
plt.legend()
plt.savefig(PLOTS_DIR / 'pca_elbow_plot.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"✓ Using {N_PCS} PCs")
"""))

# PCA tuning cell
nb2_cells.append(create_cell("markdown", """### 🎛️ Parameter Tuning: PCA

<details>
<summary>📊 Clear elbow at PC 15-20</summary>

**Action:** Current N_PCS = 15 is good
</details>

<details>
<summary>📊 Elbow at PC 30-40</summary>

**Action:**
```python
N_PCS = 35
```
**Then:** Re-run from Stage 6
</details>
"""))

nb2_cells.append(create_cell("code", """# Compute neighbors and UMAP
print("Computing neighborhood graph...")
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)

print("Running UMAP...")
sc.tl.umap(adata)

print("Clustering...")
sc.tl.leiden(adata, resolution=CLUSTERING_PARAMS['resolution'])

print(f"✓ Identified {adata.obs['leiden'].nunique()} clusters")
"""))

# Plot embeddings
nb2_cells.append(create_cell("code", """# Plot UMAP
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

sc.pl.umap(adata, color='leiden', legend_loc='on data',
          title='Clusters', ax=axes[0,0], show=False)
sc.pl.umap(adata, color='orig.ident', title='Sample',
          ax=axes[0,1], show=False)
sc.pl.umap(adata, color='Genotype', title='Genotype',
          ax=axes[1,0], show=False)
sc.pl.umap(adata, color='Sex', title='Sex',
          ax=axes[1,1], show=False)

plt.tight_layout()
plt.savefig(PLOTS_DIR / 'umap_embeddings.png', dpi=300, bbox_inches='tight')
plt.show()
"""))

# Clustering tuning
nb2_cells.append(create_cell("markdown", """### 🎛️ Parameter Tuning: Clustering

<details>
<summary>📊 <5 clusters</summary>

**Action:**
```python
CLUSTERING_PARAMS['resolution'] = 1.0  # Increase
```
</details>

<details>
<summary>📊 >30 clusters</summary>

**Action:**
```python
CLUSTERING_PARAMS['resolution'] = 0.4  # Decrease
```
</details>
"""))

# Stage 7: Markers
nb2_cells.append(create_cell("markdown", """## Stage 7: Marker Gene Analysis

Identify differentially expressed genes per cluster."""))

nb2_cells.append(create_cell("code", """# Compute markers
print("Computing marker genes...")
sc.tl.rank_genes_groups(adata, groupby='leiden',
                        method='wilcoxon', use_raw=True)

# Plot
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
plt.savefig(PLOTS_DIR / 'top_marker_genes.png', dpi=300, bbox_inches='tight')
plt.show()

# Export to CSV
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

markers_dict = {}
for group in groups:
    markers_dict[f'Cluster_{group}_genes'] = result['names'][group][:30]
    markers_dict[f'Cluster_{group}_scores'] = result['scores'][group][:30]

markers_df = pd.DataFrame(markers_dict)
markers_df.to_csv(PLOTS_DIR.parent / 'top_markers_by_cluster.csv', index=False)

print(f"✓ Saved markers to top_markers_by_cluster.csv")
"""))

# Save
nb2_cells.append(create_cell("code", """# Save clustered data
output_file = 'outputs/clustered_data.h5ad'

adata.uns['pipeline_params']['notebook'] = '2_clustering_markers'
adata.uns['pipeline_params']['n_pcs'] = N_PCS
adata.uns['pipeline_params']['n_neighbors'] = N_NEIGHBORS
adata.uns['pipeline_params']['clustering'] = CLUSTERING_PARAMS

adata.write(output_file)

print("\\n" + "="*60)
print("NOTEBOOK 2 COMPLETE")
print("="*60)
print(f"✓ Saved: {output_file}")
print(f"  Clusters: {adata.obs['leiden'].nunique()}")
print("\\n➡️  NEXT: Open 3_annotation_export.ipynb")
"""))

# Create Notebook 2
nb2 = {
    "cells": nb2_cells,
    "metadata": source_nb['metadata'],
    "nbformat": 4,
    "nbformat_minor": 0
}

output2 = OUTPUT_DIR / "2_clustering_markers.ipynb"
save_notebook(nb2, output2)
print(f"\\n✓ Created {output2.name} with {len(nb2_cells)} cells")

# ===================================================================
# NOTEBOOK 3: Annotation & Export
# ===================================================================
print("\\n" + "="*70)
print("CREATING NOTEBOOK 3")
print("="*70)

nb3_cells = []

# Header
nb3_cells.append(create_cell("markdown", """# Notebook 3: Annotation & Export

**Cell Annotation Pipeline - Part 3 of 3**

**Stages:** 8-9
**📥 Input:** `outputs/clustered_data.h5ad`
**📤 Output:** `outputs/annotated_data.h5ad` (FINAL)

---
"""))

# Load
nb3_cells.append(create_cell("code", """import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load
print("Loading data from Notebook 2...")
adata = sc.read_h5ad('outputs/clustered_data.h5ad')

# Validate
checks = {
    'UMAP': 'X_umap' in adata.obsm,
    'Clusters': 'leiden' in adata.obs.columns,
    'Markers': 'rank_genes_groups' in adata.uns,
}

for check, passed in checks.items():
    print(f"  {'✓' if passed else '✗'} {check}")
    if not passed:
        raise ValueError(f"Missing {check} - run Notebook 2!")

print(f"\\n✓ Loaded: {adata.n_obs:,} cells, {adata.obs['leiden'].nunique()} clusters")
"""))

# Parameters
nb3_cells.append(create_cell("code", """# Annotation parameters
ANNOTATION_PARAMS = {
    'margin': 0.05,
    'label_mode': 'cell',
}

# Marker genes (mouse brain)
MARKER_GENES = {
    "Excit": ["Slc17a7", "Camk2a", "Satb2"],
    "Inhib": ["Gad1", "Gad2", "Slc6a1"],
    "Astro": ["Slc1a2", "Aqp4", "Gfap"],
    "Oligo": ["Plp1", "Mog", "Mbp"],
    "OPC": ["Pdgfra", "Cspg4"],
    "Micro": ["P2ry12", "Cx3cr1", "Csf1r"],
    "Endo": ["Pecam1", "Kdr", "Flt1"],
}

PLOTS_DIR = Path('plots/notebook3')
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

print("Annotation parameters set")
"""))

# Stage 8: Annotation
nb3_cells.append(create_cell("markdown", """## Stage 8: Cell Type Annotation

Score cells for each cell type based on marker expression."""))

nb3_cells.append(create_cell("code", """# Score cells
print("Scoring cells...")

for cell_type, genes in MARKER_GENES.items():
    available = [g for g in genes if g in adata.raw.var_names]
    if available:
        sc.tl.score_genes(adata, available,
                         score_name=f'{cell_type}_score',
                         use_raw=True)
        print(f"  ✓ {cell_type}: {len(available)}/{len(genes)} markers")

# Assign cell types
score_cols = [col for col in adata.obs.columns if col.endswith('_score')]
scores = adata.obs[score_cols]

# Apply margin
scores_sorted = np.sort(scores.values, axis=1)
max_scores = scores_sorted[:, -1]
second_scores = scores_sorted[:, -2] if scores.shape[1] > 1 else scores_sorted[:, -1]
confident = (max_scores - second_scores) > ANNOTATION_PARAMS['margin']

adata.obs['celltype'] = scores.idxmax(axis=1).str.replace('_score', '')
adata.obs.loc[~confident, 'celltype'] = 'Unlabeled'

print(f"\\n✓ {confident.sum():,} / {len(confident):,} confidently labeled")
print(f"\\nCell type distribution:")
print(adata.obs['celltype'].value_counts())
"""))

# Plot
nb3_cells.append(create_cell("code", """# Plot cell types
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sc.pl.umap(adata, color='celltype', ax=axes[0], show=False,
          title='Cell type annotations')
sc.pl.umap(adata, color='leiden', legend_loc='on data',
          ax=axes[1], show=False, title='Clusters')

plt.tight_layout()
plt.savefig(PLOTS_DIR / 'cell_type_umap.png', dpi=300, bbox_inches='tight')
plt.show()

# Composition heatmap
composition = pd.crosstab(adata.obs['leiden'], adata.obs['celltype'],
                         normalize='index')

plt.figure(figsize=(10, 6))
sns.heatmap(composition, annot=True, fmt='.2f', cmap='YlOrRd')
plt.title('Cell type composition per cluster')
plt.tight_layout()
plt.savefig(PLOTS_DIR / 'composition_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()
"""))

# Export
nb3_cells.append(create_cell("markdown", """## Final Export

Save annotated data and metadata."""))

nb3_cells.append(create_cell("code", """# Save annotated data
output_file = 'outputs/annotated_data.h5ad'

adata.uns['pipeline_params']['notebook'] = '3_annotation_export'
adata.uns['pipeline_params']['annotation'] = ANNOTATION_PARAMS

adata.write(output_file)

# Export metadata CSV
metadata_cols = ['leiden', 'celltype', 'orig.ident', 'Genotype', 'Sex',
                'n_genes_by_counts', 'total_counts', 'percent_mt']
adata.obs[metadata_cols].to_csv('outputs/cell_metadata.csv')

# Summary
summary = pd.DataFrame({
    'Metric': ['Total cells', 'Clusters', 'Cell types',
               'Median genes/cell', 'Median UMIs/cell'],
    'Value': [
        f"{adata.n_obs:,}",
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
print(f"✓ Metadata: outputs/cell_metadata.csv")
print(f"✓ Summary: outputs/analysis_summary.csv")
print(f"\\n  {adata.n_obs:,} cells")
print(f"  {adata.obs['leiden'].nunique()} clusters")
print(f"  {adata.obs['celltype'].nunique()} cell types")
print("\\n🎉 Ready for downstream analysis!")

display(summary)
"""))

# Create Notebook 3
nb3 = {
    "cells": nb3_cells,
    "metadata": source_nb['metadata'],
    "nbformat": 4,
    "nbformat_minor": 0
}

output3 = OUTPUT_DIR / "3_annotation_export.ipynb"
save_notebook(nb3, output3)
print(f"\\n✓ Created {output3.name} with {len(nb3_cells)} cells")

print("\\n" + "="*70)
print("ALL NOTEBOOKS CREATED")
print("="*70)
print(f"\\nLocation: {OUTPUT_DIR}/")
print("\\nNotebooks:")
print(f"  1. {output1.name} ({len(nb1_cells)} cells)")
print(f"  2. {output2.name} ({len(nb2_cells)} cells)")
print(f"  3. {output3.name} ({len(nb3_cells)} cells)")
print("\\n✓ All notebooks are complete and ready for Colab!")
