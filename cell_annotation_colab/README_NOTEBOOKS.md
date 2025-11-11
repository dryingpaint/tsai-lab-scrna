# Cell Annotation Pipeline - Multi-Notebook Workflow

## 📚 Complete Documentation Package

This directory contains everything needed to create and use the three-notebook pipeline for single-cell RNA-seq cell annotation.

---

## 🎯 Quick Start

### **For Using the Pipeline:**

1. **Start with Notebook 1:** `1_setup_qc_filtering.ipynb`
   - Load CellBender data
   - Calculate QC metrics
   - Detect doublets (per-sample)
   - Filter cells and genes
   - **Output:** `qc_filtered_data.h5ad`

2. **Continue with Notebook 2:** `2_clustering_markers.ipynb`
   - Load filtered data
   - Normalize and scale
   - PCA, UMAP, clustering
   - Identify marker genes
   - **Output:** `clustered_data.h5ad`

3. **Finish with Notebook 3:** `3_annotation_export.ipynb`
   - Load clustered data
   - Annotate cell types
   - Recluster subtypes
   - Export final results
   - **Output:** `annotated_data.h5ad`

---

## 📖 Reference Documents

### **Core References:**

1. **`NOTEBOOK_STRUCTURE.md`**
   - Strategic breakpoints and rationale
   - Data flow between notebooks
   - Loading strategy for each notebook
   - File management and directory structure

2. **`COMPLETE_NOTEBOOK_TEMPLATES.md`**
   - Complete cell-by-cell structure for all 3 notebooks
   - Every code cell with full implementation
   - Every tuning cell with decision trees
   - Every assessment cell with automated analysis

3. **`NOTEBOOK_TUNING_CELLS.md`**
   - Ready-to-paste tuning cells for each stage
   - Collapsible decision trees (<details> tags)
   - Code snippets for each scenario
   - Automated assessment cells

4. **`PARAMETER_TUNING_GUIDE.md`**
   - Comprehensive 35+ scenarios
   - Platform-specific recommendations
   - Tissue-specific thresholds
   - Decision matrices and workflows

### **Code Validation:**

5. **`NOTEBOOK_COMPARISON.md`**
   - Detailed comparison with original pipeline
   - What's different and why
   - Critical fixes applied

6. **`FIXES_FOR_NOTEBOOK.md`**
   - Exact code corrections
   - Side-by-side comparison
   - Critical vs optional fixes

7. **`cell_annotation_pipeline_corrected.ipynb`**
   - Source code with all corrections applied
   - Use this as reference for functions

---

## 🔧 Key Features

### ✅ **Accurate Calculations**
- Custom CellBender H5 loader (handles matrix transposition)
- Per-sample doublet detection (critical!)
- Correct filtering order (doublets removed LAST)
- Manual MT%/ribo% calculation (exact match to original)
- Uses `orig.ident` column throughout

### ✅ **Interactive Tuning**
After each stage:
- **🎛️ Tuning Cell** (collapsible decision trees)
- **📊 Assessment Cell** (automated analysis)
- Specific actions based on observed results
- Clear "re-run from X" instructions

### ✅ **Checkpointed Workflow**
- Natural save points between notebooks
- Can iterate on any stage independently
- Smaller files load faster
- Parameters saved in `.uns` for reproducibility

---

## 📂 File Organization

```
project/
├── notebooks/
│   ├── 1_setup_qc_filtering.ipynb
│   ├── 2_clustering_markers.ipynb
│   └── 3_annotation_export.ipynb
├── outputs/
│   ├── qc_filtered_data.h5ad      ← After Notebook 1
│   ├── clustered_data.h5ad        ← After Notebook 2
│   ├── annotated_data.h5ad        ← Final output
│   ├── cell_metadata.csv
│   └── analysis_summary.csv
├── plots/
│   ├── notebook1/
│   ├── notebook2/
│   └── notebook3/
└── data/
    └── [CellBender H5 files]
```

---

## 🎨 Notebook Structure

Each notebook follows the same pattern:

1. **Header** (what this notebook does)
2. **Setup** (imports, parameters)
3. **Data Loading** (with validation in notebooks 2-3)
4. **Processing Stages** (with inline documentation)
5. **After each stage:**
   - Results visualization
   - **🎛️ Tuning Cell** (decision tree)
   - **📊 Assessment Cell** (automated analysis)
6. **Save Output** (with parameters in `.uns`)
7. **Summary** (what was done, what's next)

---

## 💡 Using the Tuning Cells

### Example: After QC Metrics in Notebook 1

**Step 1:** View your QC violin plots

**Step 2:** Expand relevant scenario in tuning cell:
```
📊 Mitochondrial % shows tail extending to 20%+
   → Diagnosis: Stressed/dying cells
   → Action: Set max_mt_pct = 10
```

**Step 3:** Scroll back to Parameter Configuration section

**Step 4:** Update parameter:
```python
CELL_FILTERS['max_mt_pct'] = 10
```

**Step 5:** Re-run from Stage 2 forward

**Step 6:** Check assessment cell output:
```
✅ GOOD: Median MT% in acceptable range (5-10%)
💡 NEXT STEPS: Proceed to Stage 3
```

---

## 🔄 Iteration Workflow

### **Scenario 1: Adjusting QC**
```
Run Notebook 1 → Check results → Adjust params → Re-run Notebook 1
                                                        ↓
                                                (satisfied with QC)
                                                        ↓
                                        Run Notebook 2 with new output
```

### **Scenario 2: Adjusting Clustering**
```
Have qc_filtered_data.h5ad from previous run
                ↓
        Run Notebook 2 → Try resolution 0.6 → Re-run clustering
                ↓
        (satisfied) → Run Notebook 3
```

### **Scenario 3: Multiple Annotation Strategies**
```
Have clustered_data.h5ad
        ↓
    Try annotation with margin=0.05 → Save as annotated_v1.h5ad
        ↓
    Try annotation with margin=0.03 → Save as annotated_v2.h5ad
        ↓
    Compare results, choose best
```

---

## 🚨 Critical Differences from Standard Workflows

| Aspect | Standard | This Pipeline |
|--------|----------|---------------|
| H5 loading | `sc.read_10x_h5()` | Custom `load_cellbender_h5()` |
| Doublet detection | All cells together | **Per-sample** with threshold cap |
| Filtering order | Varies | **Doublets removed LAST** |
| MT% calculation | Automatic | Manual for exact reproducibility |
| Sample ID | Variable | **`orig.ident`** consistently |

---

## 📊 Expected Outputs

### **After Notebook 1:**
- `qc_filtered_data.h5ad` (~50-200 MB typical)
  - Contains: raw counts, QC metrics, doublet scores
  - Filtered: doublets removed, low-quality cells excluded

### **After Notebook 2:**
- `clustered_data.h5ad` (~50-200 MB typical)
  - Contains: scaled data, PCA, UMAP, clusters, markers
  - Ready for annotation

### **After Notebook 3:**
- `annotated_data.h5ad` (~50-200 MB typical)
  - Contains: everything + cell type labels, confidence scores
  - Publication-ready
- `cell_metadata.csv` (easy viewing in Excel)
- `analysis_summary.csv` (key statistics)

---

## 🎓 Learning Resources

### **Built-in Documentation:**
- Each function has detailed docstrings
- Tuning cells explain parameter effects
- Assessment cells teach what to look for
- Comments explain the "why" not just "what"

### **External References:**
- [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/)
- [Single-cell best practices](https://www.sc-best-practices.org/)
- [Original pipeline](../../../) (for comparison)

---

## 🐛 Troubleshooting

### **"Data validation failed!"**
- Did you run the previous notebook?
- Check that output file exists: `ls outputs/`
- Verify file isn't corrupted: try re-running previous notebook

### **"Missing doublet scores"**
- Notebook 1 wasn't completed fully
- Re-run Notebook 1 from start

### **"No clear elbow in PCA plot"**
- Normal for complex datasets
- Use N_PCS = 30-40
- See PARAMETER_TUNING_GUIDE.md for details

### **"Very few cells pass filters"**
- QC thresholds too stringent
- See tuning cells after Stage 2 for specific adjustments
- Check sample quality

### **"Doublet rates vary widely"**
- Normal if samples loaded at different densities
- See tuning cells after Stage 3
- May need sample-specific investigation

---

## 🎯 Next Steps After Pipeline

Once you have `annotated_data.h5ad`:

1. **Differential Expression**
   - Compare conditions/genotypes
   - Use `sc.tl.rank_genes_groups()`
   - Or export to R for limma/DESeq2

2. **Pathway Analysis**
   - GSEA with marker genes
   - GO term enrichment
   - See `limma_voom_gsea.py` in original pipeline

3. **Cell-Cell Communication**
   - CellPhoneDB
   - NicheNet
   - CellChat

4. **Trajectory Analysis**
   - Pseudotime (Monocle, PAGA)
   - RNA velocity
   - Developmental trajectories

5. **Integration**
   - Combine with other datasets
   - Compare to published atlases
   - Meta-analysis

---

## 📝 Citation

If you use this pipeline, please cite:
- Scanpy: Wolf et al., Genome Biology 2018
- Scrublet: Wolock et al., Cell Systems 2019
- Original methods from your publications

---

## 🤝 Support

For questions or issues:
1. Check PARAMETER_TUNING_GUIDE.md first
2. Review COMPLETE_NOTEBOOK_TEMPLATES.md for examples
3. Consult original pipeline code
4. Contact your bioinformatics core

---

## ✨ Summary

This multi-notebook workflow provides:
- ✅ Accurate calculations matching original pipeline
- ✅ Interactive parameter tuning at every stage
- ✅ Automated assessment and recommendations
- ✅ Clear checkpoints for efficient iteration
- ✅ Complete documentation and examples
- ✅ Educational content for learning
- ✅ Production-ready for publication

**Start with Notebook 1 and follow the tuning cells at each stage!**
