# Prepare LinkedOmics TCGA-GBM cohort (expression + clinical) into RDS
# Files are attribute x sample (transposed) Firehose format
setwd("E:/KIMI/GBM_WEBS_Project")
dir.create("data/raw/LinkedOmics_TCGA", showWarnings = FALSE)

read_cct <- function(f) {
  m <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
  as.matrix(m)
}

# ---- Expression: HiSeq RSEM log2 (already log2-transformed) ----
expr <- read_cct("data/raw/LinkedOmics_TCGA/Human__TCGA_GBM__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz")
# remove trailing "?" unnamed row if present
expr <- expr[rownames(expr) != "?" & !is.na(rownames(expr)), ]
o <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
expr <- expr[o, ]
expr <- expr[!duplicated(rownames(expr)), ]
colnames(expr) <- gsub("\\.", "-", colnames(expr))

# ---- Clinical: transpose attribute x sample ----
clin_raw <- read.delim("data/raw/LinkedOmics_TCGA/Human__TCGA_GBM__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
                       row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
clin_t <- as.data.frame(t(clin_raw), stringsAsFactors = FALSE)
clin_t$sample <- gsub("\\.", "-", rownames(clin_t))

clinical <- data.frame(
  sample = clin_t$sample,
  OS.time = suppressWarnings(as.numeric(clin_t$overall_survival)),  # days
  OS = suppressWarnings(as.integer(clin_t$status)),
  age = suppressWarnings(as.numeric(clin_t$years_to_birth)),
  gender = clin_t$gender,
  tumor_purity = suppressWarnings(as.numeric(clin_t$Tumor_purity)),
  stringsAsFactors = FALSE
)
clinical$OS.time <- clinical$OS.time / 30.44  # months
clinical <- clinical[!is.na(clinical$OS.time) & !is.na(clinical$OS) & clinical$OS.time > 0, ]

# Match expression & clinical
common <- intersect(colnames(expr), clinical$sample)
expr <- expr[, common]
clinical <- clinical[match(common, clinical$sample), ]

saveRDS(expr, "data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
saveRDS(clinical, "data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")

cat("TCGA-GBM (LinkedOmics):", nrow(expr), "genes x", ncol(expr), "samples,",
    sum(clinical$OS), "death events\n")
cat("Median OS:", round(median(clinical$OS.time), 1), "months\n")

# ---- RPPA proteome (for LMS analysis) ----
rppa <- read.delim("data/raw/LinkedOmics_TCGA/Human__TCGA_GBM__MDA__RPPA__MDA_RPPA__01_28_2016__BI__Gene__Firehose_RPPA.cct",
                   row.names = 1, check.names = FALSE)
rppa <- as.matrix(rppa)
colnames(rppa) <- gsub("\\.", "-", colnames(rppa))
saveRDS(rppa, "data/raw/LinkedOmics_TCGA/tcga_gbm_rppa.rds")
cat("RPPA:", nrow(rppa), "proteins x", ncol(rppa), "samples\n")
