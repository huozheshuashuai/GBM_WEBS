# =============================================================================
# GBM WEBS Project - Environment Setup
# =============================================================================
# R version: 4.5.0 (recommended)
# Run this first to install all required packages
# =============================================================================

# Core CRAN packages
cran_pkgs <- c(
  "Seurat", "SeuratObject", "Signac", "Slingshot",
  "glmnet", "survival", "survminer", "rms",
  "ggplot2", "dplyr", "ggpubr", "patchwork",
  "pheatmap", "RColorBrewer", "viridis",
  "fgsea", "clusterProfiler", "enrichplot",
  "igraph", "ggraph", "reshape2",
  "ggrepel", "gridExtra", "cowplot",
  "timeROC", "survIDINRI", "dcurves",
  "clinfun",  # for Jonckheere-Terpstra
  "caret", "pROC",
  "BiocManager"
)

# Bioconductor packages
bioc_pkgs <- c(
  "AUCell", "UCell", "DESeq2", "TradeRnaSeq",
  "fgsea", "ComplexHeatmap", "EnhancedVolcano",
  "AnnotationHub", "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "chromVAR", "motifmatchr", "TFBSTools"
)

# Install missing CRAN packages
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
}

install_if_missing(cran_pkgs)

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# Install from GitHub (if needed)
# remotes::install_github("satijalab/seurat", ref = "seurat5")
# remotes::install_github("satijalab/signac")

message("All packages installed successfully!")
message("R version: ", R.version.string)
