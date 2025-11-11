# Parameter Tuning Guide: Data-Driven Recommendations

This guide provides specific, actionable recommendations for tuning parameters based on what you observe in the results.

---

## Stage 2: QC Metrics - When to Adjust

### After viewing QC violin plots:

#### **Genes per cell (n_genes_by_counts)**

**If you see:**
- **Two distinct populations**: One low (<500) and one higher (2000-5000)
  - ✅ **Action**: Set `min_genes = 200-300` to exclude the low population
  - ✅ **Reason**: Low population likely represents empty droplets or dead cells

- **Long tail extending to >10,000 genes**
  - ✅ **Action**: Lower `max_genes = 6000-7000`
  - ✅ **Reason**: High gene count cells are likely doublets; be more stringent

- **Main population centered around 3000-4000 with tight distribution**
  - ✅ **Action**: Current parameters OK, proceed

- **Very wide distribution (500-8000) without clear peaks**
  - ⚠️ **Action**: Check if this is biological (e.g., different cell types) or technical
  - ⚠️ **Consider**: Batch effects, sample quality issues

**Specific thresholds:**
```python
# Conservative (high quality only):
CELL_FILTERS['min_genes'] = 500
CELL_FILTERS['max_genes'] = 6000

# Standard (balanced):
CELL_FILTERS['min_genes'] = 200
CELL_FILTERS['max_genes'] = 8000

# Permissive (include more cells):
CELL_FILTERS['min_genes'] = 100
CELL_FILTERS['max_genes'] = 10000
```

---

#### **Mitochondrial % (percent_mt)**

**If you see:**
- **Main population at 2-5% with tail extending to 20%+**
  - ✅ **Action**: Set `max_mt_pct = 8-10%` to exclude tail
  - ✅ **Reason**: Tail represents stressed/dying cells

- **Bimodal distribution**: One peak at ~5%, another at ~15%
  - ⚠️ **Action**: Investigate second peak - could be:
    - Biological (certain cell types naturally higher MT)
    - Technical (sample degradation)
  - ⚠️ **Decision**: If biological, keep `max_mt_pct = 15-20%`; if technical, use `max_mt_pct = 8%`

- **Main population already <5%**
  - ✅ **Action**: Use `max_mt_pct = 5-8%` (very high quality)
  - ✅ **Reason**: You have excellent quality cells

- **Broad distribution 5-20% without clear cutoff**
  - ⚠️ **Problem**: Sample quality issues or tissue-specific characteristics
  - ⚠️ **Action for neurons**: Use `max_mt_pct = 5%` (neurons should be low)
  - ⚠️ **Action for other tissues**: Research typical MT% for your tissue type

**Tissue-specific recommendations:**
```python
# Neurons (cortex, hippocampus):
CELL_FILTERS['max_mt_pct'] = 5

# Mixed brain (includes glia):
CELL_FILTERS['max_mt_pct'] = 10

# Peripheral tissues, high-metabolism:
CELL_FILTERS['max_mt_pct'] = 15

# Stressed/diseased tissue:
CELL_FILTERS['max_mt_pct'] = 20  # Be cautious
```

---

#### **Total counts**

**If you see:**
- **Long left tail (<1000 counts)**
  - ✅ **Action**: Increase `min_counts = 1500-2000`
  - ✅ **Reason**: Low count cells have insufficient information

- **Outliers with >100,000 counts**
  - ✅ **Action**: Lower `max_counts = 60000-80000`
  - ✅ **Reason**: Extremely high counts may be doublets or technical artifacts

- **Tight distribution around 10,000-30,000**
  - ✅ **Action**: Current parameters OK
  - ✅ **Note**: This is ideal - consistent sequencing depth

---

### After viewing scatter plots:

#### **Counts vs Genes scatter**

**If you see:**
- **Linear relationship but some points way above the line**
  - ✅ **Action**: These are doublets; ensure doublet detection catches them
  - ✅ **Check**: After doublet detection, verify these points are flagged

- **Two distinct linear relationships**
  - ⚠️ **Means**: Two populations with different RNA content
  - ⚠️ **Could be**:
    - Biological (e.g., neurons vs glia)
    - Technical (batch effect)
  - ⚠️ **Action**: Don't over-filter; investigate after clustering

- **Curvature at high counts** (points flatten out)
  - ✅ **Normal**: Saturation effect at high sequencing depth
  - ✅ **Action**: No adjustment needed

#### **Counts vs MT% scatter**

**If you see:**
- **High MT% cells clustered at LOW counts**
  - ✅ **Action**: Set `max_mt_pct = 10%` to exclude dying cells
  - ✅ **Reason**: Classic signature of dying cells (lost cytoplasmic RNA, retained MT RNA)

- **High MT% cells have NORMAL counts**
  - ⚠️ **Action**: Investigate - may be biologically relevant
  - ⚠️ **Example**: Highly metabolic cell types (e.g., cardiomyocytes, muscle)
  - ⚠️ **Decision**: Consider keeping with higher `max_mt_pct`

- **No correlation between counts and MT%**
  - ✅ **Good**: No obvious cell stress pattern
  - ✅ **Action**: Proceed with standard `max_mt_pct = 10%`

---

### After filtering: Expected cell retention

**Typical retention rates:**
- **Excellent quality**: 70-85% of cells pass filters
- **Good quality**: 50-70% of cells pass filters
- **Problematic**: <50% of cells pass filters → investigate data quality

**If you retain <40% of cells:**
1. Check QC plots for anomalies
2. Verify CellBender ran correctly
3. Consider more permissive thresholds
4. Check for sample-specific issues

---

## Stage 3: Doublet Detection - When to Adjust

### After viewing per-sample histograms:

#### **If doublet scores show clear separation**
```
Example: Peak at 0.1-0.2, gap, then doublets at 0.4-0.6
```
- ✅ **Action**: Use `manual_threshold = 0.35` (default)
- ✅ **Reason**: Clear separation allows confident calling

#### **If doublet scores are bimodal with overlap**
```
Example: Singlet peak at 0.15, doublet peak at 0.35, overlapping at 0.25-0.3
```
- ⚠️ **Action**:
  - Conservative: `manual_threshold = 0.25` (removes more, some false positives)
  - Permissive: `manual_threshold = 0.40` (keeps more, some false negatives)
- ⚠️ **Decision**: Depends on downstream analysis
  - If doing DE: Be more stringent (0.25-0.30)
  - If exploratory: Be more permissive (0.35-0.40)

#### **If most cells have low scores (<0.3) with few outliers**
```
Example: Tight peak at 0.05-0.15, then scattered cells at 0.5+
```
- ✅ **Action**: Use `manual_threshold = 0.30-0.35`
- ✅ **Reason**: Outliers are clear doublets; threshold doesn't matter much

#### **If scores are continuously distributed (no clear peaks)**
```
Example: Smooth distribution from 0.0 to 0.6
```
- ⚠️ **Problem**: Difficult to separate singlets from doublets
- ⚠️ **Causes**:
  - Very homogeneous cell population (hard to detect doublets)
  - Poor sample quality
  - Too few cells in sample (<500)
- ⚠️ **Action**:
  1. Use expected doublet rate: `threshold = expected_doublet_rate + 0.05`
  2. Visually inspect UMAP after clustering for doublet clusters
  3. Consider stricter threshold (`0.25-0.30`) and validate with marker genes

---

### After checking per-sample doublet rates:

#### **If rates vary widely across samples**
```
Example: Sample A: 5%, Sample B: 15%, Sample C: 8%
```
- ⚠️ **Investigate**:
  - Check loading concentration (Sample B may have been overloaded)
  - Check cell counts per sample
  - Verify sample quality
- ⚠️ **Action**:
  - If Sample B had higher cell loading: Expected, keep current settings
  - If Sample B has quality issues: Consider excluding or being more stringent

**Decision tree:**
```
Sample doublet rate > 15%?
├─ Yes: Was this sample loaded at higher concentration?
│  ├─ Yes: OK, expected
│  └─ No: Consider stricter threshold (0.25-0.30) or exclude sample
└─ No: Proceed with current settings
```

#### **If all samples have very low doublet rates (<3%)**
- ⚠️ **Warning**: Detection may be too permissive
- ⚠️ **Action**: Lower `manual_threshold = 0.25-0.30`
- ⚠️ **Reason**: Typical 10x rates are 6-10%; <3% is suspicious

#### **If all samples have very high doublet rates (>20%)**
- ⚠️ **Warning**: Detection may be too stringent OR sample quality poor
- ⚠️ **Actions**:
  1. Check histograms - is there clear separation?
  2. Try raising `manual_threshold = 0.40-0.45`
  3. Check if `expected_doublet_rate` is too high
  4. Consider sample quality issues

---

### Platform-specific expected rates:

```python
# 10x Chromium v2 (older):
DOUBLET_PARAMS['expected_doublet_rate'] = 0.06
DOUBLET_PARAMS['manual_threshold'] = 0.35

# 10x Chromium v3 (standard):
DOUBLET_PARAMS['expected_doublet_rate'] = 0.08
DOUBLET_PARAMS['manual_threshold'] = 0.35

# 10x Chromium Next GEM (high throughput):
DOUBLET_PARAMS['expected_doublet_rate'] = 0.10
DOUBLET_PARAMS['manual_threshold'] = 0.35

# 10x Chromium with >10,000 cells/sample:
DOUBLET_PARAMS['expected_doublet_rate'] = 0.10
DOUBLET_PARAMS['manual_threshold'] = 0.30  # More stringent
```

---

## Stage 4: Filtering Summary - Post-Filter Checks

### After filtering, check cell counts:

#### **Total cells retained**

**If you have <5,000 cells total:**
- ⚠️ **Problem**: May not have enough power for clustering
- ⚠️ **Actions**:
  1. Review filters - were they too stringent?
  2. Consider relaxing `min_genes`, `max_mt_pct`
  3. Check if doublet removal was too aggressive

**If you have 5,000-20,000 cells:**
- ✅ **Good**: Sufficient for most analyses
- ✅ **Action**: Proceed

**If you have >50,000 cells:**
- ✅ **Excellent**: Great statistical power
- ℹ️ **Note**: Clustering may take longer; consider subsampling for initial exploration

---

### Cells per sample distribution:

#### **If one sample has <500 cells after filtering:**
```
Example: All samples: 2000-3000 cells, Sample X: 200 cells
```
- ⚠️ **Problem**: Sample X may have quality issues
- ⚠️ **Actions**:
  1. Check QC plots specifically for Sample X
  2. Check doublet rate for Sample X
  3. Consider excluding Sample X if it's an outlier
  4. If excluding: Update metadata and re-run from Stage 1

#### **If samples have balanced cell counts (CV < 30%):**
- ✅ **Good**: No major sample-specific issues
- ✅ **Action**: Proceed

#### **If cell counts vary 5-fold or more:**
```
Example: Sample A: 500 cells, Sample B: 3000 cells
```
- ⚠️ **Investigate**:
  - Was loading concentration different?
  - Did one sample have lower quality?
  - Biological difference (e.g., one sample has fewer viable cells)?
- ⚠️ **For downstream analysis**: Consider this imbalance when interpreting results

---

### Gene count distribution after filtering:

#### **Median genes per cell:**

**If median < 2000:**
- ⚠️ **Low**: May have low sequencing depth or high background
- ⚠️ **Check**: Were sequencing reads sufficient?
- ⚠️ **Action**:
  - If depth was adequate: Likely biological (low-complexity cell types)
  - If depth was low: Consider requesting deeper sequencing

**If median 2000-5000:**
- ✅ **Good**: Standard for most scRNA-seq
- ✅ **Action**: Proceed

**If median >6000:**
- ✅ **Excellent**: Very high quality
- ✅ **Note**: Rich dataset for detecting subtle differences

---

## Stage 6: PCA/UMAP - When to Adjust

### After viewing PCA elbow plot:

#### **If elbow is clear at PC 20:**
```
Variance explained drops sharply after PC 20, then plateaus
```
- ✅ **Action**: Set `N_PCS = 20-25`
- ✅ **Reason**: PCs after 20 mainly capture noise

#### **If variance decreases gradually (no clear elbow):**
```
Smooth decrease from PC 1 to PC 50
```
- ⚠️ **Means**: Complex dataset with many sources of variation
- ⚠️ **Action**: Use `N_PCS = 30-40`
- ⚠️ **Reason**: Need more PCs to capture biological variation

#### **If elbow is at PC 10:**
```
Sharp drop after PC 10
```
- ⚠️ **Unusual**: May indicate:
  - Very homogeneous cell population
  - Over-filtering removed diversity
  - Technical issues (batch effects dominating)
- ⚠️ **Action**:
  - Use `N_PCS = 10-15`
  - Check PC loadings to understand what's captured
  - Review if filtering was too stringent

---

### After initial UMAP:

#### **If cells cluster by sample (not by cell type):**
```
UMAP colored by 'sample' shows distinct groups per sample
```
- ⚠️ **Problem**: Batch effects dominate biological signal
- ⚠️ **Actions**:
  1. **First try**: Increase `N_PCS = 30-40` (capture more variation)
  2. **If still present**: Use batch correction (Harmony, scVI, or Seurat integration)
  3. **Check**: Were samples processed in different batches?

#### **If UMAP shows no clear structure (cloud of points):**
```
All cells mixed together, no visible clusters
```
- ⚠️ **Possible causes**:
  - `N_NEIGHBORS` too high (over-smoothing)
  - `N_PCS` too low (insufficient variation captured)
  - Truly homogeneous population
- ⚠️ **Actions**:
  1. Lower `N_NEIGHBORS = 5-10`
  2. Increase `N_PCS = 30-40`
  3. Try different resolution for clustering
  4. Check if population is expected to be homogeneous

#### **If UMAP shows many tiny scattered clusters:**
```
Dozens of small, isolated groups
```
- ⚠️ **Means**:
  - `N_NEIGHBORS` too low (under-smoothing)
  - Low-quality cells creating noise
- ⚠️ **Actions**:
  1. Increase `N_NEIGHBORS = 15-20`
  2. Review if QC was too permissive
  3. Lower clustering resolution

#### **If UMAP shows clear, separated clusters:**
- ✅ **Good**: Proceed to clustering and annotation
- ✅ **Action**: Current parameters working well

---

## Stage 7: Clustering - Resolution Selection

### After clustering at resolution 0.8:

#### **If you get <5 clusters:**
```
Example: 3 major clusters
```
- ⚠️ **Means**: Under-clustering
- ⚠️ **Action**:
  - Increase resolution: Try `0.6 → 1.0 → 1.2` in steps
  - Stop when biologically meaningful splits appear

#### **If you get 5-20 clusters:**
- ✅ **Good**: Likely biologically meaningful
- ✅ **Action**: Validate with marker genes

#### **If you get >30 clusters:**
```
Example: 45 small clusters
```
- ⚠️ **Means**: Over-clustering
- ⚠️ **Action**:
  - Decrease resolution: Try `0.8 → 0.6 → 0.4` in steps
  - Look for clusters with similar markers (can merge)

---

### Cluster size distribution:

#### **If most clusters have >100 cells:**
- ✅ **Good**: Sufficient cells for DE analysis
- ✅ **Action**: Proceed

#### **If many clusters have <20 cells:**
```
Example: 10 clusters with <20 cells each
```
- ⚠️ **Problem**:
  - Over-clustering
  - Low-quality cells forming singleton clusters
- ⚠️ **Actions**:
  1. Decrease resolution to `0.4-0.6`
  2. Check if small clusters are doublets (high counts, mixed markers)
  3. Consider merging similar small clusters

#### **If one large cluster (>50% of cells) and many tiny ones:**
```
Example: Cluster 0: 10,000 cells, Clusters 1-10: <100 cells each
```
- ⚠️ **Means**: Uneven clustering
- ⚠️ **Actions**:
  1. Check if large cluster is homogeneous (marker genes)
  2. If heterogeneous: Increase resolution
  3. If tiny clusters are outliers: May be low-quality cells or rare populations

---

## Stage 8: Marker Genes - Validation Checklist

### Top markers per cluster:

#### **If top markers are ribosomal/mitochondrial genes:**
```
Example: Cluster 5 top markers: Rpl3, Rps5, mt-Co1, mt-Nd4
```
- ⚠️ **Problem**: Low-quality cluster or stressed cells
- ⚠️ **Actions**:
  1. Check cluster size (likely small)
  2. Check QC metrics for this cluster (`percent_mt`, `percent_ribo`)
  3. If high MT%: These are dying cells → exclude cluster
  4. Consider stricter QC filters and re-run

#### **If top markers are cell cycle genes:**
```
Example: Cluster 8 top markers: Mki67, Top2a, Cdk1
```
- ✅ **Means**: Proliferating cells (normal in some tissues)
- ✅ **Actions**:
  - If expected (e.g., neurogenesis, development): Keep cluster
  - If not expected: Check if these are doublets or contamination
  - Consider regressing out cell cycle effects if it obscures biology

#### **If clusters have overlapping top markers:**
```
Example: Clusters 3 and 5 share 8/10 top markers
```
- ⚠️ **Means**: Over-clustering; these should probably be merged
- ⚠️ **Actions**:
  1. Decrease resolution
  2. Manually merge clusters post-hoc
  3. Use hierarchical clustering to identify mergeable clusters

#### **If clusters have distinct, specific markers:**
```
Example: Cluster 1: Slc17a7, Camk2a; Cluster 2: Gad1, Gad2
```
- ✅ **Excellent**: Clear biological identities
- ✅ **Action**: Proceed to annotation

---

## Stage 9: Cell Type Annotation - Quality Checks

### After annotation:

#### **If >20% cells are "Unlabeled":**
```
Example: 25% of cells have no confident cell type assignment
```
- ⚠️ **Actions**:
  1. Lower confidence `margin = 0.03` (less stringent)
  2. Check if marker genes are appropriate for your tissue
  3. Add tissue-specific markers to `MARKER_GENES`
  4. Consider cluster-level annotation (`label_mode='cluster'`)

#### **If most clusters are a single cell type:**
```
Example: All 15 clusters are annotated as "Excit"
```
- ⚠️ **Problem**:
  - Markers not specific enough
  - Over-clustering of homogeneous population
- ⚠️ **Actions**:
  1. Add more specific subtype markers
  2. Check if subtypes exist (examine top cluster markers)
  3. Decrease clustering resolution

#### **If one cluster has mixed cell types:**
```
Example: Cluster 5: 40% Excit, 30% Inhib, 30% Astro
```
- ⚠️ **Problem**:
  - Possible doublet cluster
  - Under-clustering
  - Transitional/intermediate state
- ⚠️ **Actions**:
  1. Check doublet scores for this cluster
  2. Increase clustering resolution
  3. Examine marker gene expression on UMAP
  4. Consider if biological intermediate state

#### **If annotation matches UMAP spatial organization:**
```
Similar cell types cluster together on UMAP
```
- ✅ **Excellent**: Biologically coherent
- ✅ **Action**: Proceed with confidence

---

## Quick Decision Matrix

| Observation | Most Likely Cause | First Action | If That Fails |
|-------------|-------------------|--------------|---------------|
| <50% cells pass filters | QC too stringent | Relax thresholds | Check sample quality |
| Sample clustering on UMAP | Batch effects | Increase N_PCS | Use batch correction |
| >30 tiny clusters | Over-clustering | Lower resolution | Increase N_NEIGHBORS |
| <5 large clusters | Under-clustering | Increase resolution | Check if homogeneous |
| High doublet rate (>20%) | Threshold too low | Raise threshold to 0.40 | Check sample loading |
| Low doublet rate (<3%) | Threshold too high | Lower threshold to 0.25 | Verify with UMAP |
| Many unlabeled cells | Markers not specific | Add tissue markers | Use cluster-level mode |
| MT genes as top markers | Low-quality cells | Exclude cluster | Stricter QC |

---

## Iterative Tuning Workflow

1. **First pass**: Use default parameters
2. **Check QC plots** → Adjust QC thresholds if needed → Re-run from Stage 2
3. **Check doublet histograms** → Adjust threshold → Re-run from Stage 3
4. **Check PCA elbow** → Adjust N_PCS → Re-run from Stage 6
5. **Check UMAP** → Adjust N_NEIGHBORS if needed → Re-run from Stage 6
6. **Check clustering** → Adjust resolution → Re-run clustering only
7. **Check markers** → Validate biology → If issues, go back to step 2

**Time-saving tip**: Only go back to the earliest stage that needs adjustment. Most parameters can be tuned from their respective stage forward.
