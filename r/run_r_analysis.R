#!/usr/bin/env Rscript

# Script to run the R analysis files
# Usage: Rscript run_r_analysis.R

# Set working directory
setwd("/Users/melissadu/projects/tsai-lab-scrna")

# Install required packages if not already installed
required_packages <- c(
  "Matrix", "stringr", "ggplot2", "dplyr", "tidyr", "readr",
  "Seurat", "scCustomize", "ggridges", "patchwork", "pheatmap", 
  "forcats", "purrr", "rWikiPathways"
)

bioc_packages <- c(
  "scDblFinder", "SingleCellExperiment", "edgeR", "limma",
  "fgsea", "msigdbr", "GSEABase", "AnnotationDbi", "org.Mm.eg.db"
)

# Check and install CRAN packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Check and install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Render the R Markdown files
rmarkdown::render("r/1_cellbender_qc_annotation.Rmd", 
                  output_file = "1_cellbender_qc_annotation.html",
                  output_dir = "results/")

rmarkdown::render("r/2_limma_voom_gsea.Rmd", 
                  output_file = "2_limma_voom_gsea.html",
                  output_dir = "results/")

cat("Analysis complete! Check the results/ directory for HTML reports.\n")
