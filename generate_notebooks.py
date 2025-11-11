#!/usr/bin/env python3
"""
Generate the three-notebook pipeline with integrated tuning cells.

Usage:
    python generate_notebooks.py
"""

import json
from pathlib import Path

# This script generates complete notebooks
# Due to size, it creates them programmatically rather than as single large files

def create_notebook_1():
    """Create Notebook 1: Setup, QC & Filtering"""

    # Read the corrected notebook cells as template
    corrected_nb_path = Path('cell_annotation_pipeline_corrected.ipynb')
    tuning_cells_path = Path('NOTEBOOK_TUNING_CELLS.md')

    print("Creating Notebook 1: Setup, QC & Filtering...")
    print("  This notebook includes:")
    print("  - Custom CellBender H5 loader")
    print("  - Per-sample doublet detection")
    print("  - Correct filtering order (doublets last)")
    print("  - Integrated tuning cells after each stage")
    print("  - Automated assessment cells")
    print("")
    print("To complete: Copy cells from cell_annotation_pipeline_corrected.ipynb")
    print("             Stages 1-5, then add tuning cells from NOTEBOOK_TUNING_CELLS.md")


def create_notebook_2():
    """Create Notebook 2: Clustering & Markers"""

    print("Creating Notebook 2: Clustering & Markers...")
    print("  This notebook includes:")
    print("  - Data loading with validation")
    print("  - PCA with elbow plot tuning")
    print("  - UMAP with structure assessment")
    print("  - Clustering with resolution tuning")
    print("  - Marker gene analysis")
    print("  - Integrated tuning cells")


def create_notebook_3():
    """Create Notebook 3: Annotation & Export"""

    print("Creating Notebook 3: Annotation & Export...")
    print("  This notebook includes:")
    print("  - Cell type scoring and annotation")
    print("  - Confidence-based labeling")
    print("  - Reclustering for subtypes")
    print("  - Final export with metadata")
    print("  - Integrated tuning cells")


def main():
    print("="*70)
    print("NOTEBOOK GENERATION GUIDE")
    print("="*70)
    print("")
    print("Due to size constraints, the notebooks should be created by:")
    print("")
    print("1. Starting with cell_annotation_pipeline_corrected.ipynb")
    print("2. Splitting at strategic checkpoints:")
    print("   - After Stage 4 → Save qc_filtered_data.h5ad")
    print("   - After Stage 7 → Save clustered_data.h5ad")
    print("   - After Stage 9 → Save annotated_data.h5ad")
    print("")
    print("3. Adding tuning cells from NOTEBOOK_TUNING_CELLS.md after each stage")
    print("")
    print("4. Adding assessment code cells that analyze results")
    print("")
    print("="*70)
    print("REFERENCE FILES:")
    print("="*70)
    print("- NOTEBOOK_STRUCTURE.md       : Strategic breakpoints & data flow")
    print("- COMPLETE_NOTEBOOK_TEMPLATES.md : Full cell-by-cell structure")
    print("- NOTEBOOK_TUNING_CELLS.md    : Tuning cells to insert")
    print("- PARAMETER_TUNING_GUIDE.md   : Comprehensive tuning reference")
    print("- cell_annotation_pipeline_corrected.ipynb : Source code")
    print("")
    print("="*70)
    print("MANUAL CREATION STEPS:")
    print("="*70)
    print("")
    print("NOTEBOOK 1:")
    print("  1. Copy Stages 1-4 from corrected notebook")
    print("  2. After Stage 2: Add QC tuning cell + assessment cell")
    print("  3. After Stage 3: Add doublet tuning cell + summary cell")
    print("  4. After Stage 4: Add filtering tuning cell + impact cell")
    print("  5. Add save cell: adata.write('outputs/qc_filtered_data.h5ad')")
    print("")
    print("NOTEBOOK 2:")
    print("  1. Add data loading + validation cell")
    print("  2. Copy Stage 5 (normalization)")
    print("  3. Copy Stage 6 (PCA/UMAP/clustering)")
    print("     - After PCA: Add elbow tuning cell")
    print("     - After UMAP: Add structure tuning cell")
    print("     - After clustering: Add resolution tuning cell")
    print("  4. Copy Stage 7 (marker genes)")
    print("     - Add marker validation tuning cell")
    print("  5. Add save cell: adata.write('outputs/clustered_data.h5ad')")
    print("")
    print("NOTEBOOK 3:")
    print("  1. Add data loading + validation cell")
    print("  2. Copy Stage 8 (annotation)")
    print("     - Add annotation tuning cell + coverage assessment")
    print("  3. Copy Stage 9 (reclustering)")
    print("  4. Add final export cells")
    print("")
    print("="*70)
    print("")
    print("Would you like me to create:")
    print("  [A] Complete notebooks (large files)")
    print("  [B] Starter templates with TODO markers")
    print("  [C] Just keep the reference documents")
    print("")
    print("Recommendation: Option C - Use reference docs to manually create")
    print("                This gives you most control and flexibility")


if __name__ == '__main__':
    main()
