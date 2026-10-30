# Repair pass 2: force Windows BINARY packages (avoid source compilation)
options(timeout = 600, Ncpus = 8, pkgType = "binary")
library(BiocManager)

failed <- c(
  "AnnotationDbi", "GO.db",
  "Signac", "Slingshot", "SingleCellExperiment", "tradeSeq",
  "AUCell", "UCell", "DESeq2",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "AnnotationHub", "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "chromVAR", "motifmatchr", "TFBSTools", "JASPAR2020",
  "GEOquery", "TCGAbiolinks", "RcisTarget", "GENIE3",
  "easierData", "ExperimentHub",
  "GSVA", "forestplot", "spdep", "STRINGdb"
)

for (pkg in failed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing (binary):", pkg, "\n")
    tryCatch(BiocManager::install(pkg, update = FALSE, ask = FALSE, type = "binary"),
             error = function(e) cat("  error:", e$message, "\n"))
  }
  ok <- require(pkg, character.only = TRUE, quietly = TRUE)
  cat(sprintf("[%s] %s\n", ifelse(ok, "OK", "FAIL"), pkg), "\n",
      file = "results/package_repair_log.txt", append = TRUE)
}
cat("Repair pass 2 complete.\n")
