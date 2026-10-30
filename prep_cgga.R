# Prepare CGGA cohorts into RDS for downstream scripts (base R only)
setwd("E:/KIMI/GBM_WEBS_Project")

prep_cgga <- function(expr_file, clin_file, out_prefix) {
  expr <- read.delim(expr_file, row.names = 1, check.names = FALSE)
  expr <- as.matrix(expr)
  # collapse duplicate gene names by max mean
  o <- order(rowMeans(expr), decreasing = TRUE)
  expr <- expr[o, ]
  expr <- expr[!duplicated(rownames(expr)), ]

  clin <- read.delim(clin_file, check.names = FALSE)
  colnames(clin)[1] <- "sample"
  censor_col <- grep("Censor", colnames(clin), value = TRUE)[1]
  out <- data.frame(
    sample = clin$sample,
    OS.time = as.numeric(clin$OS) / 30.44,  # days -> months
    OS = as.numeric(clin[[censor_col]]),
    Age = as.numeric(clin$Age),
    Gender = clin$Gender,
    Grade = clin$Grade,
    Histology = clin$Histology,
    PRS_type = clin$PRS_type,
    IDH_mutation = clin$IDH_mutation_status
  )
  out <- out[!is.na(out$OS.time) & !is.na(out$OS) & out$OS.time > 0, ]
  expr <- expr[, intersect(colnames(expr), out$sample)]
  out <- out[match(colnames(expr), out$sample), ]

  saveRDS(expr, paste0(out_prefix, "_expr.rds"))
  saveRDS(out, paste0(out_prefix, "_clinical.rds"))
  cat(out_prefix, ": ", nrow(expr), "genes x", ncol(expr), "samples,",
      sum(out$OS), "events\n")
}

dir.create("data/raw/CGGA", showWarnings = FALSE)
prep_cgga("data/raw/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",
          "data/raw/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt",
          "data/raw/CGGA/cgga693")
prep_cgga("data/raw/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
          "data/raw/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt",
          "data/raw/CGGA/cgga325")
