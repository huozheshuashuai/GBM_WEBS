# =============================================================================
# Corrected package installer for GBM WEBS Project (R 4.5.2 / Bioc 3.21)
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org/"))
options(timeout = 600)
options(Ncpus = 8)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.21", update = FALSE, ask = FALSE)

# ---- CRAN packages ----
cran_pkgs <- c(
  "Seurat", "SeuratObject", "Signac",
  "glmnet", "survival", "survminer", "rms",
  "ggplot2", "dplyr", "ggpubr", "patchwork",
  "pheatmap", "RColorBrewer", "viridis",
  "igraph", "ggraph", "reshape2",
  "ggrepel", "gridExtra", "cowplot",
  "timeROC", "survIDINRI", "dcurves",
  "clinfun", "caret", "pROC",
  "mgcv", "meta", "remotes", "hdf5r"
)

# ---- Bioconductor packages ----
bioc_pkgs <- c(
  "Slingshot", "SingleCellExperiment", "tradeSeq",
  "AUCell", "UCell", "DESeq2",
  "fgsea", "ComplexHeatmap", "EnhancedVolcano",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "AnnotationHub", "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "chromVAR", "motifmatchr", "TFBSTools", "JASPAR2020",
  "limma", "edgeR", "survcomp",
  "GEOquery", "TCGAbiolinks",
  "RcisTarget", "GENIE3",
  "easierData", "ExperimentHub"
)

log_result <- function(pkg, ok) {
  cat(sprintf("[%s] %s\n", ifelse(ok, "OK", "FAIL"), pkg),
      file = "results/package_install_log.txt", append = TRUE)
}

dir.create("results", showWarnings = FALSE)
cat("=== Install log ===\n", file = "results/package_install_log.txt")

for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch(install.packages(pkg, dependencies = TRUE), error = function(e) NULL)
  }
  log_result(pkg, require(pkg, character.only = TRUE, quietly = TRUE))
}

for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch(BiocManager::install(pkg, update = FALSE, ask = FALSE), error = function(e) NULL)
  }
  log_result(pkg, require(pkg, character.only = TRUE, quietly = TRUE))
}

# ---- GitHub: spacexr (RCTD) ----
if (!require("spacexr", quietly = TRUE)) {
  tryCatch(remotes::install_github("dmcable/spacexr", upgrade = "never"),
           error = function(e) message("spacexr install failed: ", e$message))
}
log_result("spacexr", require("spacexr", quietly = TRUE))

cat("\n=== Summary ===\n")
print(R.version.string)
print(BiocManager::version())
