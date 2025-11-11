# Tuning Recommendation Cells to Add After Each Stage

Add these markdown + code cells after each stage in the notebook for specific parameter tuning guidance.

---

## After Stage 2: QC Metrics

### Markdown Cell:

```markdown
### 🎛️ Parameter Tuning: What to Adjust Based on Your Results

Review the plots above and use this decision tree to adjust parameters:

#### **Genes per Cell**

**What did you observe?**

<details>
<summary>📊 Two distinct populations (low ~500 and high ~3000)</summary>

**Diagnosis**: Empty droplets or dead cells in low population

**Action**:
```python
# Go back to Parameter Configuration section and update:
CELL_FILTERS['min_genes'] = 300  # Increase to exclude low population
```

**Then**: Re-run from Stage 2
</details>

<details>
<summary>📊 Long tail extending to >10,000 genes</summary>

**Diagnosis**: Likely contains doublets

**Action**:
```python
# Go back to Parameter Configuration section and update:
CELL_FILTERS['max_genes'] = 6000  # Lower threshold
```

**Then**: Re-run from Stage 2
</details>

<details>
<summary>📊 Tight distribution around 2000-5000 genes</summary>

**Diagnosis**: ✅ Good quality!

**Action**: Current parameters work well, proceed to next stage
</details>

---

#### **Mitochondrial %**

**What did you observe?**

<details>
<summary>📊 Main population at 2-5% with tail to 20%+</summary>

**Diagnosis**: Tail represents stressed/dying cells

**Action**:
```python
# Standard filtering:
CELL_FILTERS['max_mt_pct'] = 10

# Or more stringent for neurons:
CELL_FILTERS['max_mt_pct'] = 5
```

**Then**: Re-run from Stage 2
</details>

<details>
<summary>📊 Bimodal distribution (peaks at ~5% and ~15%)</summary>

**Diagnosis**: Could be biological OR technical

**Investigate**: Check if second peak correlates with specific samples or cell types

**Actions**:
```python
# If biological (certain cell types):
CELL_FILTERS['max_mt_pct'] = 15

# If technical (sample degradation):
CELL_FILTERS['max_mt_pct'] = 8
```

**Then**: Re-run from Stage 2
</details>

<details>
<summary>📊 Most cells already <5%</summary>

**Diagnosis**: ✅ Excellent quality!

**Action**:
```python
CELL_FILTERS['max_mt_pct'] = 5  # Can be stringent
```

**Then**: Re-run from Stage 2
</details>

---

#### **Expected Cell Retention**

Check how many cells passed the filters above.
```

### Code Cell:

```python
# Cell retention analysis
print("\\n" + "="*60)
print("CELL RETENTION ANALYSIS")
print("="*60)

# This assumes adata_original was saved before filtering
# If you want to track this, add after Stage 2 QC calculation:
# adata_original = adata.copy()

# For now, we can estimate from the current state
retention_checks = {
    'Total cells in filtered dataset': f"{adata.n_obs:,}",
    'Total genes in filtered dataset': f"{adata.n_vars:,}",
    'Median genes per cell': f"{adata.obs['n_genes_by_counts'].median():.0f}",
    'Median UMI counts per cell': f"{adata.obs['total_counts'].median():.0f}",
    'Median MT%': f"{adata.obs['percent_mt'].median():.2f}%",
}

for metric, value in retention_checks.items():
    print(f"{metric:.<50} {value}")

print("="*60)

# Interpretation guide
print("\\n📊 INTERPRETATION:")

if adata.obs['n_genes_by_counts'].median() < 2000:
    print("⚠️  WARNING: Median genes per cell is low (<2000)")
    print("   → Check if sequencing depth was sufficient")
    print("   → May be low-complexity cell types (expected for some tissues)")
elif adata.obs['n_genes_by_counts'].median() > 6000:
    print("✅ EXCELLENT: Very high genes per cell (>6000)")
    print("   → High-quality, information-rich dataset")
else:
    print("✅ GOOD: Median genes per cell in normal range (2000-6000)")

if adata.obs['percent_mt'].median() > 10:
    print("\\n⚠️  WARNING: Median MT% is high (>10%)")
    print("   → May indicate cell stress or sample quality issues")
    print("   → Consider more stringent max_mt_pct threshold")
elif adata.obs['percent_mt'].median() < 5:
    print("\\n✅ EXCELLENT: Low median MT% (<5%)")
    print("   → High-quality viable cells")
else:
    print("\\n✅ GOOD: Median MT% in acceptable range (5-10%)")

print("\\n💡 NEXT STEPS:")
print("   • If metrics look good: Proceed to Stage 3 (Doublet Detection)")
print("   • If issues detected: Adjust parameters above and re-run Stage 2")
```
```

---

## After Stage 3: Doublet Detection

### Markdown Cell:

```markdown
### 🎛️ Parameter Tuning: Doublet Detection Results

Review the per-sample histograms and summary statistics above.

#### **Per-Sample Doublet Rates**

**What did you observe?**

<details>
<summary>📊 Rates between 6-10% across all samples</summary>

**Diagnosis**: ✅ Normal for 10x Chromium

**Action**: Current threshold working well, proceed to Stage 4
</details>

<details>
<summary>📊 One or more samples >15% doublets</summary>

**Diagnosis**: Sample may have been overloaded or has quality issues

**Investigate**:
1. Check if high-doublet samples had more cells loaded
2. Check QC metrics specifically for those samples

**Actions**:
```python
# If overloading was intentional:
# → Keep current settings, expected behavior

# If unexpected:
# → Consider excluding problematic samples OR
# → Use stricter threshold for those samples
DOUBLET_PARAMS['manual_threshold'] = 0.30  # More stringent
```

**Then**: Re-run from Stage 3
</details>

<details>
<summary>📊 All samples <3% doublets</summary>

**Diagnosis**: ⚠️ Threshold may be too permissive

**Action**:
```python
# Lower threshold to detect more doublets:
DOUBLET_PARAMS['manual_threshold'] = 0.25
```

**Then**: Re-run from Stage 3

**Verify**: After clustering, check if intermediate clusters appear (likely doublets)
</details>

<details>
<summary>📊 All samples >20% doublets</summary>

**Diagnosis**: ⚠️ Threshold too stringent OR poor sample quality

**Actions**:
```python
# Try more permissive threshold:
DOUBLET_PARAMS['manual_threshold'] = 0.40

# OR adjust expected rate:
DOUBLET_PARAMS['expected_doublet_rate'] = 0.08  # Lower expectation
```

**Then**: Re-run from Stage 3
</details>

---

#### **Score Distributions (from histograms)**

**What did you observe in the histograms?**

<details>
<summary>📊 Clear separation: singlets peak ~0.15, gap, doublets peak ~0.45</summary>

**Diagnosis**: ✅ Excellent - clear singlet/doublet populations

**Action**: Current threshold (0.35) is optimal, proceed
</details>

<details>
<summary>📊 Overlapping peaks: singlets ~0.15, doublets ~0.35, overlap at 0.25-0.30</summary>

**Diagnosis**: Difficult to separate - will have some false positives/negatives

**Decision**:
```python
# For differential expression (be stringent):
DOUBLET_PARAMS['manual_threshold'] = 0.25

# For exploratory analysis (be permissive):
DOUBLET_PARAMS['manual_threshold'] = 0.40
```

**Then**: Re-run from Stage 3
</details>

<details>
<summary>📊 No clear peaks: smooth distribution from 0.0 to 0.6</summary>

**Diagnosis**: ⚠️ Difficult dataset
- May have homogeneous cell population
- Poor sample quality
- Too few cells

**Actions**:
1. Use expected rate-based threshold:
```python
DOUBLET_PARAMS['manual_threshold'] = DOUBLET_PARAMS['expected_doublet_rate'] + 0.05
```
2. Plan to validate with UMAP after clustering
3. Check for doublet clusters using marker genes

**Then**: Re-run from Stage 3, verify in Stage 6-7
</details>
```

### Code Cell:

```python
# Doublet detection summary analysis
print("\\n" + "="*60)
print("DOUBLET DETECTION SUMMARY")
print("="*60)

# Per-sample doublet statistics
doublet_summary = adata.obs.groupby('orig.ident').agg({
    'predicted_doublet': ['count', 'sum'],
    'doublet_score': ['mean', 'median', 'max']
}).round(3)

doublet_summary.columns = ['n_cells', 'n_doublets', 'mean_score', 'median_score', 'max_score']
doublet_summary['pct_doublets'] = (doublet_summary['n_doublets'] / doublet_summary['n_cells'] * 100).round(1)

print("\\nPer-sample statistics:")
print(doublet_summary)

# Overall statistics
overall_rate = (adata.obs['predicted_doublet'].sum() / len(adata.obs)) * 100
mean_sample_rate = doublet_summary['pct_doublets'].mean()
cv_sample_rate = (doublet_summary['pct_doublets'].std() / doublet_summary['pct_doublets'].mean()) * 100

print("\\n" + "="*60)
print(f"Overall doublet rate: {overall_rate:.1f}%")
print(f"Mean per-sample rate: {mean_sample_rate:.1f}%")
print(f"CV across samples: {cv_sample_rate:.1f}%")
print("="*60)

# Interpretation
print("\\n📊 INTERPRETATION:")

if overall_rate < 3:
    print("⚠️  WARNING: Very low doublet rate (<3%)")
    print("   → Consider lowering threshold to 0.25-0.30")
    print("   → Verify after clustering that no doublet clusters remain")
elif overall_rate > 20:
    print("⚠️  WARNING: Very high doublet rate (>20%)")
    print("   → Check if threshold (0.35) is too stringent")
    print("   → Investigate sample quality")
    print("   → Consider raising threshold to 0.40")
elif 6 <= overall_rate <= 10:
    print("✅ EXCELLENT: Doublet rate in expected range (6-10%)")
    print("   → Typical for 10x Chromium platform")
else:
    print("✅ GOOD: Doublet rate acceptable (3-15%)")

if cv_sample_rate > 50:
    print("\\n⚠️  HIGH VARIABILITY: Doublet rates vary widely across samples")
    print(f"   → CV = {cv_sample_rate:.1f}%")
    print("   → Investigate high-doublet samples:")
    high_doublet_samples = doublet_summary[doublet_summary['pct_doublets'] > 15]
    if len(high_doublet_samples) > 0:
        print(f"   → Samples with >15% doublets: {list(high_doublet_samples.index)}")
else:
    print(f"\\n✅ CONSISTENT: Doublet rates similar across samples (CV = {cv_sample_rate:.1f}%)")

print("\\n💡 NEXT STEPS:")
print("   • If rates look normal: Proceed to Stage 4 (Filtering)")
print("   • If rates are concerning: Adjust threshold and re-run Stage 3")
print("   • After clustering: Validate by checking for doublet clusters on UMAP")
```
```

---

## After Stage 4: Filtering

### Markdown Cell:

```markdown
### 🎛️ Parameter Tuning: Filtering Results

Review the filtering summary and cell retention statistics above.

#### **Total Cell Retention**

**How many cells do you have after filtering?**

<details>
<summary>📊 <5,000 cells total</summary>

**Diagnosis**: ⚠️ May not have enough power for fine-grained clustering

**Actions**:
1. Review if filters were too stringent
2. Consider relaxing:
```python
CELL_FILTERS['min_genes'] = 150  # Was 200
CELL_FILTERS['max_mt_pct'] = 12  # Was 10
DOUBLET_PARAMS['manual_threshold'] = 0.40  # Was 0.35
```

3. Check if doublet removal was too aggressive

**Then**: Re-run from appropriate stage
</details>

<details>
<summary>📊 5,000-20,000 cells</summary>

**Diagnosis**: ✅ Sufficient for most analyses

**Action**: Proceed to Stage 5 (Normalization)
</details>

<details>
<summary>📊 >50,000 cells</summary>

**Diagnosis**: ✅ Excellent statistical power!

**Note**:
- Clustering may take longer (5-10 min)
- Consider subsampling for initial exploration if >100k cells
- Full dataset recommended for final analysis

**Action**: Proceed to Stage 5 (Normalization)
</details>

---

#### **Per-Sample Cell Counts**

**Check the bar plot above: Are sample sizes balanced?**

<details>
<summary>📊 One sample has <500 cells (others have >2000)</summary>

**Diagnosis**: ⚠️ Outlier sample with potential quality issues

**Investigate**:
```python
# Check QC metrics for the low-count sample
low_sample = 'D25-2675'  # Replace with your sample ID
sample_qc = adata.obs[adata.obs['orig.ident'] == low_sample][
    ['n_genes_by_counts', 'total_counts', 'percent_mt', 'doublet_score']
].describe()
print(sample_qc)
```

**Actions**:
- If QC metrics are poor: Consider excluding this sample
- If QC metrics are normal: May be biological (fewer viable cells in this sample)

**To exclude a sample**:
```python
# Go back to Parameter Configuration and update:
SAMPLE_NAMES = [s for s in SAMPLE_NAMES if s != 'D25-2675']  # Remove problematic sample
```

**Then**: Re-run from Stage 1
</details>

<details>
<summary>📊 All samples have 1000-5000 cells (balanced)</summary>

**Diagnosis**: ✅ No major sample-specific issues

**Action**: Proceed to Stage 5
</details>

<details>
<summary>📊 Cell counts vary >5-fold across samples</summary>

**Diagnosis**: Large imbalance

**Consider**:
- Was loading concentration different? (Expected)
- Sample quality differences? (Investigate)
- Biological differences? (Document)

**Impact**: Downstream analyses may be biased toward high-count samples

**Action**: Proceed but document this for interpretation
</details>
```

### Code Cell:

```python
# Filtering impact analysis
print("\\n" + "="*60)
print("FILTERING IMPACT ANALYSIS")
print("="*60)

# Create detailed summary
filtering_summary = {
    'Metric': [],
    'Value': [],
    'Assessment': []
}

# Total cells
total_cells = adata.n_obs
filtering_summary['Metric'].append('Total cells after filtering')
filtering_summary['Value'].append(f"{total_cells:,}")
if total_cells < 5000:
    filtering_summary['Assessment'].append('⚠️ Low - consider relaxing filters')
elif total_cells > 50000:
    filtering_summary['Assessment'].append('✅ Excellent statistical power')
else:
    filtering_summary['Assessment'].append('✅ Sufficient for analysis')

# Median QC metrics
median_genes = adata.obs['n_genes_by_counts'].median()
filtering_summary['Metric'].append('Median genes per cell')
filtering_summary['Value'].append(f"{median_genes:.0f}")
if median_genes < 2000:
    filtering_summary['Assessment'].append('⚠️ Low - check sequencing depth')
elif median_genes > 5000:
    filtering_summary['Assessment'].append('✅ Excellent quality')
else:
    filtering_summary['Assessment'].append('✅ Normal range')

median_mt = adata.obs['percent_mt'].median()
filtering_summary['Metric'].append('Median MT%')
filtering_summary['Value'].append(f"{median_mt:.2f}%")
if median_mt > 10:
    filtering_summary['Assessment'].append('⚠️ High - consider stricter filter')
elif median_mt < 5:
    filtering_summary['Assessment'].append('✅ Excellent quality')
else:
    filtering_summary['Assessment'].append('✅ Acceptable range')

# Sample balance
sample_counts = adata.obs['orig.ident'].value_counts()
min_sample = sample_counts.min()
max_sample = sample_counts.max()
cv_samples = (sample_counts.std() / sample_counts.mean()) * 100

filtering_summary['Metric'].append('Sample size range')
filtering_summary['Value'].append(f"{min_sample:,} - {max_sample:,}")
if min_sample < 500:
    filtering_summary['Assessment'].append('⚠️ Some samples very small')
elif cv_samples > 50:
    filtering_summary['Assessment'].append('⚠️ High variability across samples')
else:
    filtering_summary['Assessment'].append('✅ Reasonably balanced')

# Display summary
summary_df = pd.DataFrame(filtering_summary)
print("\\n")
display(summary_df)

print("\\n" + "="*60)
print("💡 RECOMMENDATIONS:")
print("="*60)

# Generate recommendations
recommendations = []

if total_cells < 5000:
    recommendations.append("• Consider relaxing QC filters to retain more cells")
    recommendations.append("  → Increase max_mt_pct by 2-5%")
    recommendations.append("  → Decrease min_genes by 50-100")
    recommendations.append("  → Increase doublet threshold to 0.40")

if median_genes < 2000:
    recommendations.append("• Low median genes per cell detected")
    recommendations.append("  → Check if sequencing depth was sufficient")
    recommendations.append("  → May be expected for certain cell types")

if median_mt > 10:
    recommendations.append("• High median MT% detected")
    recommendations.append("  → Consider lowering max_mt_pct to 8%")
    recommendations.append("  → May indicate sample quality issues")

if min_sample < 500:
    recommendations.append("• Some samples have very few cells")
    recommendations.append(f"  → Samples with <500 cells: {list(sample_counts[sample_counts < 500].index)}")
    recommendations.append("  → Consider excluding low-count samples")

if len(recommendations) == 0:
    print("✅ All metrics look good! Proceed to Stage 5 (Normalization)")
else:
    print("⚠️ Issues detected:")
    for rec in recommendations:
        print(rec)
    print("\\n   Decide whether to:")
    print("   → Adjust parameters and re-run")
    print("   → Proceed with current data (document limitations)")
```
```

---

## After Stage 6: PCA/UMAP

### Markdown Cell:

```markdown
### 🎛️ Parameter Tuning: PCA and UMAP Results

Review the PCA elbow plot and UMAP embeddings above.

#### **PCA Elbow Plot**

**Where is the elbow (where variance plateaus)?**

<details>
<summary>📊 Clear elbow at PC 15-20</summary>

**Diagnosis**: ✅ Standard complexity

**Action**: Current `N_PCS = 15` is appropriate, proceed
</details>

<details>
<summary>📊 Elbow at PC 30-40</summary>

**Diagnosis**: Complex dataset with many sources of variation

**Action**:
```python
# Update in Parameter Configuration:
N_PCS = 30  # Or 35, 40 based on elbow
```

**Then**: Re-run from Stage 6
</details>

<details>
<summary>📊 Elbow at PC 8-10</summary>

**Diagnosis**: ⚠️ Unusual - very homogeneous OR over-filtered

**Investigate**:
- Check PC loadings to see what's captured
- May be single cell type or batch effect dominant

**Action**:
```python
N_PCS = 10  # Use fewer PCs
```

**Then**: Re-run from Stage 6
</details>

<details>
<summary>📊 No clear elbow (gradual decrease)</summary>

**Diagnosis**: Many dimensions capture meaningful variation

**Action**:
```python
N_PCS = 30  # Use more PCs to capture variation
```

**Then**: Re-run from Stage 6
</details>

---

#### **UMAP Structure**

**What pattern do you see when colored by sample?**

<details>
<summary>📊 Cells cluster by sample (batch effects)</summary>

**Diagnosis**: ⚠️ Batch effects dominate biological signal

**Actions**:
1. First try increasing PCs:
```python
N_PCS = 40  # Capture more variation
```

2. If still present, use batch correction:
```python
# Option A: Harmony (simple)
import harmonypy
adata_corrected = harmonypy.run_harmony(adata, 'orig.ident')

# Option B: Scanpy's combat
import scanpy.external as sce
sce.pp.harmony_integrate(adata, 'orig.ident')
```

**Then**: Re-run from Stage 6
</details>

<details>
<summary>📊 Cells mix well across samples (no batch effect)</summary>

**Diagnosis**: ✅ Excellent - biological signal dominates

**Action**: Proceed to clustering
</details>

---

**What structure do you see overall?**

<details>
<summary>📊 No clear structure (cloud of points)</summary>

**Diagnosis**: ⚠️ Over-smoothing OR homogeneous population

**Actions**:
```python
# Try lower neighbors:
N_NEIGHBORS = 5

# Try more PCs:
N_PCS = 30
```

**Then**: Re-run from Stage 6
</details>

<details>
<summary>📊 Many tiny scattered clusters</summary>

**Diagnosis**: ⚠️ Under-smoothing OR noisy data

**Actions**:
```python
# Try higher neighbors:
N_NEIGHBORS = 20

# Check if QC was too permissive
```

**Then**: Re-run from Stage 6
</details>

<details>
<summary>📊 Clear, separated clusters</summary>

**Diagnosis**: ✅ Good structure for clustering

**Action**: Proceed to clustering (Stage 7)
</details>
```

---

## After Stage 7: Clustering and Markers

### Markdown Cell:

```markdown
### 🎛️ Parameter Tuning: Clustering and Marker Results

Review the cluster UMAP and top marker genes above.

#### **Number of Clusters**

**How many clusters were identified?**

<details>
<summary>📊 <5 clusters</summary>

**Diagnosis**: ⚠️ Under-clustering

**Action**:
```python
# Increase resolution:
CLUSTERING_PARAMS['resolution'] = 1.0  # Or 1.2
```

**Then**: Re-run clustering (can skip PCA/UMAP)
</details>

<details>
<summary>📊 5-20 clusters</summary>

**Diagnosis**: ✅ Likely biologically meaningful

**Action**: Validate with marker genes, proceed if they make sense
</details>

<details>
<summary>📊 >30 clusters</summary>

**Diagnosis**: ⚠️ Over-clustering

**Action**:
```python
# Decrease resolution:
CLUSTERING_PARAMS['resolution'] = 0.4  # Or 0.6
```

**Then**: Re-run clustering
</details>

---

#### **Cluster Sizes**

**Check cluster size distribution:**

<details>
<summary>📊 Most clusters have >100 cells</summary>

**Diagnosis**: ✅ Good - sufficient for DE analysis

**Action**: Proceed to annotation
</details>

<details>
<summary>📊 Many clusters have <20 cells</summary>

**Diagnosis**: ⚠️ Over-clustering or outlier cells

**Actions**:
```python
# Lower resolution:
CLUSTERING_PARAMS['resolution'] = 0.4

# Check if tiny clusters are doublets or low-quality
```

**Then**: Re-run clustering
</details>

---

####  **Top Marker Genes**

**What are the top markers for your clusters?**

<details>
<summary>📊 Ribosomal/mitochondrial genes (Rpl*, Rps*, mt-*)</summary>

**Diagnosis**: ⚠️ Low-quality cluster

**Actions**:
1. Check QC metrics for this cluster
2. If high MT%: Exclude cluster or use stricter QC
3. Consider this cluster as artifact

**To exclude**:
```python
# After clustering:
adata = adata[adata.obs['leiden'] != '5']  # Replace '5' with cluster ID
```
</details>

<details>
<summary>📊 Cell cycle genes (Mki67, Top2a, Cdk1)</summary>

**Diagnosis**: Proliferating cells

**Decision**:
- If expected (development, neurogenesis): Keep
- If obscuring biology: Regress out cell cycle

**To regress**:
```python
# Before normalization (Stage 5):
cell_cycle_genes = ['Mki67', 'Top2a', 'Cdk1', 'Mcm5', 'Pcna']
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
```
</details>

<details>
<summary>📊 Distinct, specific markers per cluster</summary>

**Diagnosis**: ✅ Excellent - biologically meaningful

**Action**: Proceed to annotation (Stage 8)
</details>

<details>
<summary>📊 Overlapping markers between clusters</summary>

**Diagnosis**: ⚠️ Over-clustering

**Action**:
```python
# Decrease resolution:
CLUSTERING_PARAMS['resolution'] = 0.5
```

**Then**: Re-run clustering
</details>
```

---

## After Stage 8: Cell Type Annotation

### Markdown Cell:

```markdown
### 🎛️ Parameter Tuning: Annotation Results

Review the cell type UMAP and composition heatmap above.

#### **Unlabeled Cells**

**What percentage of cells are "Unlabeled"?**

<details>
<summary>📊 <10% unlabeled</summary>

**Diagnosis**: ✅ Excellent annotation coverage

**Action**: Proceed to analysis
</details>

<details>
<summary>📊 10-20% unlabeled</summary>

**Diagnosis**: ✅ Acceptable - may be rare types or intermediates

**Action**: Can proceed or refine annotations
</details>

<details>
<summary>📊 >20% unlabeled</summary>

**Diagnosis**: ⚠️ Poor annotation coverage

**Actions**:
1. Lower confidence margin:
```python
ANNOTATION_PARAMS['margin'] = 0.03  # Was 0.05
```

2. Add tissue-specific markers to `MARKER_GENES`

3. Try cluster-level annotation:
```python
ANNOTATION_PARAMS['label_mode'] = 'cluster'  # Was 'cell'
```

**Then**: Re-run annotation
</details>

---

#### **Annotation-UMAP Agreement**

**Do similar cell types cluster together spatially?**

<details>
<summary>📊 Yes - similar types are near each other on UMAP</summary>

**Diagnosis**: ✅ Biologically coherent

**Action**: Proceed with confidence
</details>

<details>
<summary>📊 No - same type scattered across UMAP</summary>

**Diagnosis**: ⚠️ Possible issues:
- Markers not specific enough
- Batch effects
- Heterogeneous annotation

**Actions**:
1. Check marker gene specificity
2. Review clustering resolution
3. Check for batch effects
</details>

---

#### **Cluster Purity**

**Check the composition heatmap: Are clusters mostly one cell type?**

<details>
<summary>📊 Most clusters are >80% single cell type</summary>

**Diagnosis**: ✅ Pure clusters - good clustering

**Action**: Proceed to downstream analysis
</details>

<details>
<summary>📊 Several clusters are mixed (30-40% of 2-3 types)</summary>

**Diagnosis**: ⚠️ Possible issues:
- Doublet cluster
- Under-clustering
- Transitional state

**Actions**:
1. Check doublet scores for mixed clusters
2. Increase clustering resolution
3. Examine marker expression carefully
</details>

---

#### **Expected Cell Types**

**Do you see the cell types you expect for your tissue?**

<details>
<summary>📊 Yes - major expected types are present</summary>

**Diagnosis**: ✅ Successful annotation

**Action**: Proceed to biological interpretation
</details>

<details>
<summary>📊 No - missing expected cell types</summary>

**Diagnosis**: ⚠️ Possible issues:
- Types were filtered out
- Types are present but unlabeled
- Markers not appropriate

**Actions**:
1. Check if types might be in "Unlabeled"
2. Add specific markers for missing types
3. Review if QC filters were too stringent
</details>

<details>
<summary>📊 Unexpected cell types appear</summary>

**Diagnosis**: Could be:
- Contamination
- Mis-annotation
- Interesting biology

**Actions**:
1. Check marker expression for unexpected types
2. Review literature for tissue composition
3. Consider if markers need refinement
</details>
```

This provides specific, actionable guidance at each stage based on what the user actually observes in their data!
