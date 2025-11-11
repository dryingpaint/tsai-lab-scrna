# Comparison: Notebook vs Original Pipeline

## Critical Differences Found

### 1. **Data Loading Method**

**Original Pipeline:**
- Uses custom `load_cellbender_h5()` function that reads H5 files with h5py
- Handles matrix transposition explicitly (genes × cells → cells × genes)
- Uses `anndata.concat()` with `join="outer"` and `fill_value=0`
- Adds sample prefix to cell barcodes: `f"{sample}_{barcode}"`

**Notebook:**
- Uses `sc.read_10x_h5()` directly (simpler but assumes standard 10x format)
- Uses `adata.concatenate()` with `batch_key` and `index_unique='_'`
- May not handle CellBender-specific H5 structure correctly

**Issue:** CellBender outputs may have different H5 structure than standard 10x files
**Fix Required:** Use the custom loading function from `utils/data_loader.py`

---

### 2. **Metadata Addition**

**Original Pipeline (`utils/data_loader.py`):**
```python
metadata_df = pd.DataFrame({
    "orig.ident": sample_ids,
    "Genotype": ["E3", "E4", "E3", "E4"] * 4,
    "Stimulation": ["Ctrl"] * 8 + ["GENUS"] * 8,
    "Sex": ["M", "M", "F", "F"] * 4,
})
```
- Uses `orig.ident` column (not `sample`)
- Has specific metadata pattern for 16 samples

**Notebook:**
```python
metadata = pd.DataFrame({
    'sample': sample_names,
    'Genotype': ['WT'] * len(sample_names),  # Placeholder
    'Sex': ['M'] * len(sample_names),
    'Stimulation': ['Control'] * len(sample_names),
})
```
- Uses placeholder values
- Maps from `sample` column

**Issue:** Different column name and placeholder metadata
**Fix Required:** Use `orig.ident` consistently, document actual metadata pattern

---

### 3. **QC Metrics Calculation**

**Original Pipeline (`utils/qc_utils.py:35-41`):**
```python
adata.obs["percent_mt"] = (
    adata[:, adata.var["mt"]].X.sum(axis=1).A1 / adata.obs["total_counts"]
) * 100

adata.obs["percent_ribo"] = (
    adata[:, adata.var["ribo"]].X.sum(axis=1).A1 / adata.obs["total_counts"]
) * 100
```
- Manual calculation using matrix operations
- `.A1` converts sparse matrix to 1D array

**Notebook:**
```python
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)
adata.obs['percent_mt'] = adata.obs['pct_counts_mt']
adata.obs['percent_ribo'] = adata.obs['pct_counts_ribo']
```
- Uses scanpy's built-in function with `qc_vars` parameter
- Renames columns from `pct_counts_*` to `percent_*`

**Status:** Both should produce identical results
**Note:** Original doesn't use `var_type="genes"` parameter

---

### 4. **Doublet Detection - Per-Sample Processing**

**Original Pipeline (`utils/improved_doublet_detection.py`):**
- Processes each sample separately using `sample_col="orig.ident"`
- Uses `.copy()` when subsetting: `adata_sample = adata[mask].copy()`
- Implements manual threshold capping at 0.4 if auto-threshold too high
- Stores results in pre-allocated arrays for all cells
- Plots histograms for all samples (up to 16)

**Notebook:**
- Only shows simplified doublet detection
- Doesn't process per-sample
- Missing the threshold capping logic

**Issue:** Critical difference - doublets should be detected per-sample
**Fix Required:** Use the full `detect_doublets_improved()` function

---

### 5. **Doublet Detection Workflow**

**Original Pipeline (`cellbender_qc_annotation.py:86-117`):**
```python
# Step 1: Apply basic QC filters BEFORE doublet detection
adata_for_doublets = adata[
    (adata.obs.n_genes_by_counts >= CELL_FILTERS["min_genes"])
    & (adata.obs.n_genes_by_counts <= CELL_FILTERS["max_genes"])
    & (adata.obs.percent_mt <= CELL_FILTERS["max_mt_pct"])
].copy()

# Step 2: Detect doublets on QC-filtered cells
adata_for_doublets = detect_doublets_improved(...)

# Step 3: Transfer doublet annotations back to original adata
adata.obs["doublet_score"] = 0.0
adata.obs["predicted_doublet"] = False
adata.obs.loc[adata_for_doublets.obs.index, "doublet_score"] = ...
adata.obs.loc[adata_for_doublets.obs.index, "predicted_doublet"] = ...

# Step 4: Then apply full filtering including doublet removal
adata = filter_cells_and_genes(...)
```

**Notebook:**
- Implements the same workflow structure ✓
- Properly subsets data for doublet detection ✓
- Transfers results back to full dataset ✓

**Status:** Workflow is correct

---

### 6. **Cell and Gene Filtering Order**

**Original Pipeline (`utils/qc_utils.py:200-226`):**
```python
# 1. Filter cells by min_genes
sc.pp.filter_cells(adata, min_genes=min_genes)
# 2. Filter genes by min_cells
sc.pp.filter_genes(adata, min_cells=GENE_FILTERS["min_cells"])
# 3. Filter cells by max_genes
adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
# 4. Filter cells by max_mt_pct
adata = adata[adata.obs.percent_mt < max_mt_pct, :]
# 5. Optional count filters
# 6. Remove doublets LAST
adata = adata[~adata.obs.predicted_doublet, :]
```

**Notebook:**
```python
# 1. Remove doublets FIRST
adata = adata[~adata.obs['predicted_doublet']].copy()
# 2. Filter cells by min_genes
sc.pp.filter_cells(adata, min_genes=min_genes)
# 3. Filter by max_genes, counts, MT%, etc.
# 4. Filter genes
sc.pp.filter_genes(adata, min_cells=GENE_FILTERS['min_cells'])
```

**Issue:** Different filtering order
- Original: Removes doublets LAST (after all QC filters)
- Notebook: Removes doublets FIRST

**Impact:** Minor - results should be similar but not identical
**Recommendation:** Match original order for consistency

---

### 7. **Normalization and Scaling**

**Both versions:**
- Save raw counts
- Normalize to 10,000 reads per cell
- Log1p transform
- Find highly variable genes
- Keep only HVGs
- Scale with max_value=10

**Status:** Identical ✓

---

### 8. **PCA and Clustering**

**Original Pipeline:**
- Fixed resolution of 0.8 (determined from previous sweep)
- `auto_resolution` code is commented out
- Uses `choose_leiden_resolution()` function when enabled

**Notebook:**
- Shows auto_resolution workflow but simplified
- Doesn't implement full `choose_leiden_resolution()` function
- Uses fixed resolution 0.8 as fallback

**Status:** Simplified but uses same final resolution ✓

---

### 9. **Marker Gene Scoring**

**Original Pipeline (`utils/annotation.py`):**
- Uses `sc.tl.score_genes()` with `use_raw=True`
- Implements complex assignment logic with margin thresholds
- Has separate functions for cluster-level vs cell-level scoring
- Includes refinement for subtypes (excitatory layers, inhibitory subtypes)

**Notebook:**
- Shows simplified scoring approach
- Missing the full assignment logic
- Doesn't implement subtype refinement

**Issue:** Simplified version missing key logic
**Fix Required:** Note that full implementation is in utils/annotation.py

---

## Summary of Required Fixes

### High Priority (Affects Results)
1. ✅ **Data Loading**: Use custom `load_cellbender_h5()` function
2. ✅ **Doublet Detection**: Use full `detect_doublets_improved()` with per-sample processing
3. ✅ **Filtering Order**: Match original order (doublets removed last)
4. ⚠️ **Cell Type Annotation**: Document that simplified version is shown, full logic in utils

### Medium Priority (Best Practices)
5. ⚠️ **Metadata**: Update to use `orig.ident` and document actual metadata structure
6. ⚠️ **QC Calculation**: Both methods work, but document difference

### Low Priority (Documentation)
7. ✓ **PCA Elbow**: Resolution selection simplified but uses correct final value
8. ✓ **Normalization**: Identical implementation

---

## Recommendation

The notebook should either:
1. **Option A**: Import and use the actual utility functions from `utils/` modules
2. **Option B**: Include disclaimer that code is simplified for readability, refer to pipeline for production

I recommend Option A for critical sections (data loading, doublet detection) and Option B with clear notes for complex sections (annotation).
