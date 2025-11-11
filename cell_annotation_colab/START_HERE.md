# Cell Annotation Pipeline - Colab Notebooks

## 🚀 Quick Start

Upload these three notebooks to Google Colab and run them in order:

1. **`1_setup_qc_filtering.ipynb`** → Outputs `qc_filtered_data.h5ad`
2. **`2_clustering_markers.ipynb`** → Outputs `clustered_data.h5ad`
3. **`3_annotation_export.ipynb`** → Outputs `annotated_data.h5ad` (FINAL)

Each notebook includes:
- ✅ Exact calculations matching the original pipeline
- ✅ Interactive parameter tuning cells (🎛️)
- ✅ Automated assessment cells (📊)
- ✅ Clear next steps

---

## 📁 What's in This Folder

### **Notebooks (Ready for Colab):**
- `1_setup_qc_filtering.ipynb` (36 cells) - Data loading, QC, doublet detection, filtering
- `2_clustering_markers.ipynb` (16 cells) - Normalization, PCA, UMAP, clustering, markers
- `3_annotation_export.ipynb` (8 cells) - Cell type annotation and final export

### **Documentation:**
- **`README_NOTEBOOKS.md`** ⭐ - Complete guide (START HERE for details)
- `NOTEBOOK_STRUCTURE.md` - Why notebooks are split this way
- `PARAMETER_TUNING_GUIDE.md` - 35+ tuning scenarios with exact actions
- `NOTEBOOK_TUNING_CELLS.md` - Reference for tuning cells
- `COMPLETE_NOTEBOOK_TEMPLATES.md` - Full structure documentation

### **Validation:**
- `NOTEBOOK_COMPARISON.md` - How notebooks compare to original pipeline
- `FIXES_FOR_NOTEBOOK.md` - Critical corrections applied
- `NOTEBOOK_FILES_SUMMARY.md` - Complete file inventory

---

## 🎯 Using the Notebooks

### **Step 1: Upload to Colab**
1. Go to [Google Colab](https://colab.research.google.com/)
2. Upload `1_setup_qc_filtering.ipynb`
3. Mount your Google Drive (if data is there)

### **Step 2: Configure Parameters**
In the "Parameter Configuration" section, update:
```python
# Data paths
BASE_PATH = "/path/to/your/data/"
SAMPLE_NAMES = [f"D25-{i}" for i in range(2675, 2691)]  # Your samples

# QC thresholds (adjust based on your data)
CELL_FILTERS = {
    'min_genes': 200,
    'max_genes': 8000,
    'max_mt_pct': 10,
    ...
}
```

### **Step 3: Run and Tune**
1. Run all cells sequentially
2. After each stage, check the **🎛️ Tuning Cell**
3. Check the **📊 Assessment Cell** for automated recommendations
4. If needed, adjust parameters and re-run from that stage

### **Step 4: Continue to Next Notebook**
Once satisfied with results, the notebook saves a checkpoint file. Load this in the next notebook.

---

## 🎨 Key Features

### **1. Interactive Tuning Cells**
After each stage, expandable decision trees guide you:

```markdown
### 🎛️ Parameter Tuning: QC Results

<details>
<summary>📊 MT% tail extends to 20%+</summary>

**Diagnosis:** Stressed/dying cells
**Action:** Set max_mt_pct = 10
**Then:** Re-run from Stage 2
</details>
```

### **2. Automated Assessment**
Code cells analyze your results and provide verdicts:

```python
if median_mt > 10:
    print("⚠️  HIGH: Consider max_mt_pct = 8")
else:
    print("✅ GOOD: Acceptable range")
```

### **3. Checkpointed Workflow**
- Save after QC/filtering
- Save after clustering
- Save final annotated data
- Iterate on any stage independently

### **4. Exact Calculations**
- Custom CellBender H5 loader
- Per-sample doublet detection
- Correct filtering order
- Matches original pipeline precisely

---

## 📊 Expected Outputs

After running all three notebooks:

```
outputs/
├── qc_filtered_data.h5ad      (~50-200 MB)
├── clustered_data.h5ad         (~50-200 MB)
├── annotated_data.h5ad         (~50-200 MB) ← FINAL
├── cell_metadata.csv
└── analysis_summary.csv

plots/
├── notebook1/
│   ├── qc_violin_plots.png
│   ├── doublet_score_histograms.png
│   └── filtered_data_summary.png
├── notebook2/
│   ├── pca_elbow_plot.png
│   ├── umap_embeddings.png
│   └── top_marker_genes.png
└── notebook3/
    ├── cell_type_umap.png
    └── composition_heatmap.png
```

---

## 💡 Tuning Workflow Example

### Scenario: High MT% in QC plots

**Step 1:** Run Notebook 1, view QC plots
- Observe: MT% distribution has tail extending to 25%

**Step 2:** Expand relevant tuning cell
```
📊 MT% tail extends to 20%+
→ Diagnosis: Stressed/dying cells
→ Action: Set max_mt_pct = 10
```

**Step 3:** Scroll back to Parameter Configuration
```python
CELL_FILTERS['max_mt_pct'] = 10  # Changed from 15
```

**Step 4:** Re-run from "Stage 2: QC Metrics" onward

**Step 5:** Check assessment cell
```
✅ GOOD: Median MT% in acceptable range (5-10%)
💡 Proceed to Stage 3
```

---

## 🔄 Iteration Patterns

### **Pattern 1: Adjust QC Filters**
```
Run Notebook 1 → See results → Adjust CELL_FILTERS → Re-run Notebook 1
                                                              ↓
                                                    (satisfied with QC)
                                                              ↓
                                             Run Notebook 2 with new file
```

### **Pattern 2: Try Different Clustering**
```
Have qc_filtered_data.h5ad
        ↓
Run Notebook 2 → Try resolution 0.6 → Re-run clustering cells
        ↓
(satisfied) → Run Notebook 3
```

### **Pattern 3: Multiple Annotation Strategies**
```
Have clustered_data.h5ad
        ↓
Try different marker sets → Save as annotated_v1.h5ad
Try different confidence → Save as annotated_v2.h5ad
        ↓
Compare and choose best
```

---

## 🚨 Critical Differences from Standard Workflows

| Aspect | Standard | This Pipeline |
|--------|----------|---------------|
| H5 loading | `sc.read_10x_h5()` | Custom loader for CellBender |
| Doublet detection | All cells | Per-sample with threshold cap |
| Filtering order | Variable | Doublets removed LAST |
| Sample ID | Variable | `orig.ident` consistently |
| MT% calculation | Automatic | Manual for reproducibility |

---

## 📚 Documentation Guide

- **Quick reference:** This file (START_HERE.md)
- **Complete guide:** README_NOTEBOOKS.md
- **Parameter help:** PARAMETER_TUNING_GUIDE.md (35+ scenarios)
- **Structural info:** NOTEBOOK_STRUCTURE.md
- **Validation:** NOTEBOOK_COMPARISON.md

---

## 🐛 Troubleshooting

### **"Data validation failed!"**
→ Previous notebook didn't complete. Check output files exist.

### **"Missing doublet scores"**
→ Notebook 1 wasn't fully run. Re-run from start.

### **"Very few cells pass filters"**
→ QC too stringent. See tuning cells after Stage 2 for adjustments.

### **"Doublet rates vary widely"**
→ Normal if samples loaded differently. See tuning cells after Stage 3.

### **"No clear elbow in PCA"**
→ Complex dataset. Use N_PCS = 30-40. See PARAMETER_TUNING_GUIDE.md.

---

## ✨ What Makes These Notebooks Special

1. **Exact Reproducibility** - Matches original pipeline calculations precisely
2. **Interactive Guidance** - Data-driven recommendations at every stage
3. **Efficient Iteration** - Checkpointed workflow, don't re-run everything
4. **Self-Documenting** - Parameters saved in files for reproducibility
5. **Educational** - Learn while analyzing
6. **Production-Ready** - Suitable for publication

---

## 🎓 Next Steps After Pipeline

Once you have `annotated_data.h5ad`:

1. **Differential Expression** - Compare conditions/genotypes
2. **Pathway Analysis** - GSEA, GO enrichment
3. **Cell-Cell Communication** - CellPhoneDB, NicheNet
4. **Trajectory Analysis** - Pseudotime, RNA velocity
5. **Integration** - Combine with other datasets

See README_NOTEBOOKS.md for detailed next steps.

---

## 📞 Need Help?

1. Check **PARAMETER_TUNING_GUIDE.md** for specific scenarios
2. Review **COMPLETE_NOTEBOOK_TEMPLATES.md** for examples
3. Consult **NOTEBOOK_COMPARISON.md** for validation questions
4. Contact your bioinformatics core

---

## 📄 Citation

If you use these notebooks, please cite:
- Scanpy: Wolf et al., Genome Biology 2018
- Scrublet: Wolock et al., Cell Systems 2019
- Your original methods papers

---

**🎉 You're ready to start! Upload `1_setup_qc_filtering.ipynb` to Colab and begin!**

For detailed instructions, see **README_NOTEBOOKS.md**
