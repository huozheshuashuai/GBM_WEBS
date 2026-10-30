# Repair pass 3: ordered install, core infrastructure first
options(timeout = 900, Ncpus = 8)
library(BiocManager)

log_ok <- function(pkg) {
  ok <- require(pkg, character.only = TRUE, quietly = TRUE)
  cat(sprintf("[%s] %s\n", ifelse(ok, "OK", "FAIL"), pkg),
      file = "results/package_repair3_log.txt", append = TRUE)
  cat(sprintf("[%s] %s\n", ifelse(ok, "OK", "FAIL"), pkg))
  ok
}

inst <- function(pkgs, type = "binary") {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(">>> Installing:", pkg, "(", type, ")\n")
      tryCatch(
        BiocManager::install(pkg, update = FALSE, ask = FALSE, type = type, force = TRUE),
        error = function(e) cat("  install error:", e$message, "\n"))
    }
    log_ok(pkg)
  }
}

unlink("results/package_repair3_log.txt")

# 1. Core infrastructure (order matters)
inst(c("GenomeInfoDbData", "GenomeInfoDb", "AnnotationDbi", "GO.db"), type = "binary")
# GO.db / GenomeInfoDbData may be source-only (data packages, no compilation needed)
inst(c("GenomeInfoDbData", "GO.db"), type = "source")

# 2. Second tier
inst(c("SingleCellExperiment", "DESeq2", "UCell", "AUCell", "Slingshot",
       "org.Hs.eg.db", "AnnotationHub", "ExperimentHub"), type = "binary")
inst(c("org.Hs.eg.db"), type = "source")

# 3. Third tier
inst(c("Signac", "tradeSeq", "clusterProfiler", "enrichplot",
       "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38",
       "chromVAR", "motifmatchr", "TFBSTools", "JASPAR2020",
       "GEOquery", "TCGAbiolinks", "RcisTarget", "GENIE3",
       "easierData", "GSVA", "forestplot", "spdep", "STRINGdb"), type = "binary")
inst(c("EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38"), type = "source")

cat("=== Repair pass 3 complete ===\n")
