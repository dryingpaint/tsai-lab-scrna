# Critical Fixes Required for Notebook

## 1. Data Loading (Stage 1)

### Current Code (INCORRECT):
```python
def load_and_merge_cellbender_data(base_path, sample_names, custom_name):
    """Load and merge CellBender H5 files from multiple samples"""
    adatas = []
    for sample in sample_names:
        filepath = f"{base_path}{sample}/{custom_name}"
        adata_sample = sc.read_10x_h5(filepath)  # ❌ WRONG
        adata_sample.var_names_make_unique()
        adata_sample.obs['sample'] = sample
        adatas.append(adata_sample)

    adata = adatas[0].concatenate(...)  # ❌ WRONG METHOD
```

### Corrected Code:
```python
# Add this helper function first
def load_cellbender_h5(file_path):
    """Load CellBender processed h5 file

    CellBender outputs may have different structure than standard 10x files.
    This function handles the specific H5 format properly.
    """
    import h5py
    from scipy import sparse
    import anndata

    with h5py.File(file_path, 'r') as f:
        # Get the matrix data
        matrix = f['matrix']
        features = f['matrix']['features']
        barcodes = f['matrix']['barcodes']
        data = f['matrix']['data']
        indices = f['matrix']['indices']
        indptr = f['matrix']['indptr']
        shape = f['matrix']['shape']

        # Read the actual values
        data_vals = data[:]
        indices_vals = indices[:]
        indptr_vals = indptr[:]
        shape_vals = tuple(shape[:])

        # Create sparse matrix
        X = sparse.csc_matrix((data_vals, indices_vals, indptr_vals), shape=shape_vals)

        # Get feature names and barcodes
        gene_names = [x.decode('utf-8') for x in features['name'][:]]
        gene_ids = [x.decode('utf-8') for x in features['id'][:]]
        cell_barcodes = [x.decode('utf-8') for x in barcodes[:]]

        # Create AnnData object (transpose if needed to get cells x genes)
        if X.shape[0] == len(gene_names) and X.shape[1] == len(cell_barcodes):
            # Matrix is genes x cells, transpose to cells x genes
            adata = anndata.AnnData(X.T.tocsr())
        else:
            # Matrix is already cells x genes
            adata = anndata.AnnData(X.tocsr())

        adata.var_names = gene_names
        adata.var['gene_ids'] = gene_ids
        adata.obs_names = cell_barcodes
        adata.var_names_make_unique()

    return adata


def load_and_merge_cellbender_data(base_path, sample_names, custom_name):
    """Load and merge CellBender H5 files from multiple samples"""
    print(f"Loading {len(sample_names)} samples...")

    adatas = []
    for sample in sample_names:
        filepath = f"{base_path}{sample}/{custom_name}"
        try:
            # Use custom loader for CellBender format
            adata_sample = load_cellbender_h5(filepath)

            # Add sample metadata
            adata_sample.obs['sample'] = sample
            adata_sample.obs['orig.ident'] = sample  # ✓ Use orig.ident

            # Add sample prefix to cell barcodes for uniqueness
            adata_sample.obs_names = [f"{sample}_{barcode}" for barcode in adata_sample.obs_names]

            adatas.append(adata_sample)
            print(f"  ✓ {sample}: {adata_sample.n_obs} cells, {adata_sample.n_vars} genes")
        except Exception as e:
            print(f"  ✗ Failed to load {sample}: {e}")

    if not adatas:
        raise ValueError("No data loaded! Check your paths.")

    # Merge using anndata.concat (not concatenate!)
    import anndata
    adata = anndata.concat(adatas, join='outer', fill_value=0)
    adata.var_names_make_unique()

    print(f"\n✓ Merged dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata
```

---

## 2. Metadata Addition (Stage 1)

### Current Code (INCORRECT):
```python
def add_metadata(adata, sample_names):
    """Add metadata columns based on sample names"""
    metadata = pd.DataFrame({
        'sample': sample_names,  # ❌ Should be 'orig.ident'
        'Genotype': ['WT'] * len(sample_names),  # ❌ Placeholder
        'Sex': ['M'] * len(sample_names),
        'Stimulation': ['Control'] * len(sample_names),
    })
```

### Corrected Code:
```python
def add_metadata(adata, sample_names):
    """Add experimental metadata to AnnData object

    ⚠️ IMPORTANT: Customize this function for your experiment!
    The pattern below is specific to the example dataset.
    """
    print("\\nAdding metadata...")

    # Example metadata pattern (for 16 samples with specific design)
    # 🔧 CUSTOMIZE THIS FOR YOUR EXPERIMENT
    if len(sample_names) == 16:
        # Original pipeline pattern
        metadata = pd.DataFrame({
            'orig.ident': sample_names,
            'Genotype': ['E3', 'E4', 'E3', 'E4'] * 4,
            'Stimulation': ['Ctrl'] * 8 + ['GENUS'] * 8,
            'Sex': ['M', 'M', 'F', 'F'] * 4,
        })
    else:
        # Generic placeholder - YOU MUST CUSTOMIZE THIS
        print("  ⚠️ WARNING: Using placeholder metadata!")
        print("  ⚠️ Edit this function to match your experimental design!")
        metadata = pd.DataFrame({
            'orig.ident': sample_names,
            'Genotype': ['Unknown'] * len(sample_names),
            'Sex': ['Unknown'] * len(sample_names),
            'Stimulation': ['Unknown'] * len(sample_names),
        })

    # Map metadata to cells using orig.ident
    for col in ['Genotype', 'Sex', 'Stimulation']:
        adata.obs[col] = adata.obs['orig.ident'].map(
            dict(zip(metadata['orig.ident'], metadata[col]))
        )

    print("  ✓ Metadata added")
    return adata
```

---

## 3. QC Metrics Calculation (Stage 2)

### Current Code (ACCEPTABLE but document):
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

### Add This Note:
```python
# NOTE: Original pipeline calculates MT% manually:
# adata.obs["percent_mt"] = (adata[:, adata.var["mt"]].X.sum(axis=1).A1 / adata.obs["total_counts"]) * 100
# Both methods produce identical results.
```

---

## 4. Doublet Detection (Stage 3) - CRITICAL FIX

### Current Code (WRONG - doesn't process per-sample):
```python
scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(...)
```

### Corrected Code:
```python
def detect_doublets_improved(adata, expected_doublet_rate=0.10, manual_threshold=0.35,
                            plot_histograms=True, save_dir=None):
    """Detect doublets using Scrublet with per-sample processing

    IMPORTANT: Doublets must be detected per-sample to account for
    sample-specific doublet rates and characteristics.
    """
    import numpy as np

    print("\\nDetecting doublets...")
    print(f"  Expected doublet rate: {expected_doublet_rate*100}%")
    print(f"  Manual threshold: {manual_threshold}")

    # Store results for all cells
    all_scores = np.zeros(adata.n_obs)
    all_predictions = np.zeros(adata.n_obs, dtype=bool)

    # Process each sample separately - CRITICAL!
    samples = adata.obs['orig.ident'].unique()

    if plot_histograms and save_dir:
        fig, axes = plt.subplots(4, 4, figsize=(16, 12))
        axes = axes.flatten()

    for idx, sample in enumerate(samples):
        print(f"\\nProcessing sample: {sample}")

        # Get sample mask
        mask = adata.obs['orig.ident'] == sample
        sample_indices = np.where(mask)[0]

        # Extract sample data (must use .copy()!)
        adata_sample = adata[mask].copy()

        # Skip if too few cells
        if adata_sample.n_obs < 100:
            print(f"  Skipping - only {adata_sample.n_obs} cells")
            continue

        # Initialize Scrublet for this sample
        scrub = scr.Scrublet(adata_sample.X, expected_doublet_rate=expected_doublet_rate)

        # Run doublet detection
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=DOUBLET_PARAMS['min_counts'],
            min_cells=DOUBLET_PARAMS['min_cells'],
            min_gene_variability_pctl=DOUBLET_PARAMS['min_gene_variability_pctl'],
            n_prin_comps=DOUBLET_PARAMS['n_prin_comps'],
            verbose=False,
        )

        # Get automatic threshold
        auto_threshold = scrub.call_doublets(threshold=None)[1]

        # Use manual threshold if specified
        if manual_threshold is not None:
            threshold = manual_threshold
            predicted_doublets = doublet_scores > threshold
        else:
            threshold = auto_threshold

        # Cap threshold to avoid missing obvious doublets
        if threshold > 0.4:
            print(f"  Warning: High auto threshold {threshold:.2f}, capping at 0.4")
            threshold = 0.4
            predicted_doublets = doublet_scores > threshold

        # Store results for this sample's cells
        all_scores[sample_indices] = doublet_scores
        all_predictions[sample_indices] = predicted_doublets

        # Statistics
        n_doublets = predicted_doublets.sum()
        pct_doublets = n_doublets / len(doublet_scores) * 100

        print(f"  Cells: {len(doublet_scores)}")
        print(f"  Threshold: {threshold:.3f}")
        print(f"  Doublets: {n_doublets} ({pct_doublets:.1f}%)")

        # Plot histogram
        if plot_histograms and save_dir and idx < 16:
            ax = axes[idx]
            ax.hist(doublet_scores, bins=50, alpha=0.7, edgecolor='black')
            ax.axvline(threshold, color='red', linestyle='--',
                      label=f'Threshold: {threshold:.2f}')
            ax.set_title(f"{sample}\\n{n_doublets} doublets ({pct_doublets:.1f}%)")
            ax.set_xlabel('Doublet Score')
            ax.set_ylabel('Frequency')
            ax.legend(fontsize=8)

    # Save plot
    if plot_histograms and save_dir:
        plt.tight_layout()
        plt.savefig(save_dir / 'doublet_score_histograms.png', dpi=300, bbox_inches='tight')
        print(f"\\n  Saved: {save_dir}/doublet_score_histograms.png")
        plt.close()

    # Add results to adata
    adata.obs['doublet_score'] = all_scores
    adata.obs['predicted_doublet'] = all_predictions

    # Overall summary
    total_doublets = all_predictions.sum()
    overall_rate = total_doublets / len(all_predictions) * 100

    print(f"\\n✓ Overall: {total_doublets} doublets ({overall_rate:.1f}%)")

    return adata
```

---

## 5. Filtering Order (Stage 4) - CRITICAL FIX

### Current Code (WRONG ORDER):
```python
# Remove doublets FIRST ❌
adata = adata[~adata.obs['predicted_doublet']].copy()
sc.pp.filter_cells(adata, min_genes=min_genes)
# ... more filters
sc.pp.filter_genes(adata, min_cells=GENE_FILTERS['min_cells'])
```

### Corrected Code:
```python
def filter_cells_and_genes(adata, min_genes=200, max_genes=8000, max_mt_pct=10,
                          min_counts=1000, max_counts=50000, max_ribo_pct=None):
    """Apply QC filtering in the correct order

    IMPORTANT: Filtering order matters!
    Doublets should be removed LAST after all other QC filters.
    """
    print("\\nFiltering cells and genes...")
    print(f"  Starting: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # 1. Filter cells by minimum genes
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print(f"  After min_genes filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    # 2. Filter genes by minimum cells expressing
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=GENE_FILTERS['min_cells'])
    print(f"  After min_cells filter: {adata.n_vars:,} genes ({n_genes_before - adata.n_vars:,} removed)")

    # 3. Filter cells by maximum genes
    n_before = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts <= max_genes].copy()  # Note: < in original, <= here
    print(f"  After max_genes filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    # 4. Filter cells by MT percentage
    n_before = adata.n_obs
    adata = adata[adata.obs.percent_mt <= max_mt_pct].copy()  # Note: < in original, <= here
    print(f"  After MT% filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    # 5. Optional: Filter by counts
    if min_counts is not None:
        n_before = adata.n_obs
        adata = adata[adata.obs.total_counts >= min_counts].copy()
        print(f"  After min_counts filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    if max_counts is not None:
        n_before = adata.n_obs
        adata = adata[adata.obs.total_counts <= max_counts].copy()
        print(f"  After max_counts filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    # 6. Optional: Filter by ribosomal percentage
    if max_ribo_pct is not None:
        n_before = adata.n_obs
        adata = adata[adata.obs.percent_ribo <= max_ribo_pct].copy()
        print(f"  After ribo% filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    # 7. Remove doublets LAST (most important!)
    n_before = adata.n_obs
    adata = adata[~adata.obs.predicted_doublet].copy()
    print(f"  After doublet removal: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")

    print(f"\\n  ✓ Final: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    return adata
```

---

## 6. Add Import Section at Top of Notebook

Add this cell after installation:

```python
# Import additional required libraries
import h5py
from scipy import sparse
import anndata

print("✓ Additional libraries imported")
```

---

## Summary of Changes

### Critical (Affects Results):
1. ✅ Data loading: Use custom `load_cellbender_h5()` function
2. ✅ Use `orig.ident` column instead of `sample`
3. ✅ Doublet detection: Process per-sample with proper threshold handling
4. ✅ Filtering order: Remove doublets LAST, not first

### Important (Best Practices):
5. ⚠️ Metadata: Document actual experimental design pattern
6. ℹ️ Add notes about calculation differences where acceptable

### Apply These Fixes:
Copy the corrected code blocks into the corresponding cells in the notebook.
The most critical sections are:
- Stage 1: Data Loading
- Stage 3: Doublet Detection
- Stage 4: Cell Filtering
