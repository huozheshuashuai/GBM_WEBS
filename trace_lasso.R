# Trace: can LASSO reproduce the paper's six genes?
# Test: TCGA vs CGGA-693 training, lambda.min vs lambda.1se
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(glmnet); library(survival) })
set.seed(42)

CANDIDATE_GENES <- c(
  "EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP",
  "HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7",
  "HIF1A", "CA9", "VEGFA", "LDHA", "LDHB", "PKM", "SLC16A3", "ENO1", "ENO2", "GAPDH",
  "MMP2", "ITGA5", "SPP1", "FN1"
)
PAPER6 <- c("HIF1A", "SIRT2", "MMP2", "ITGA5", "HDAC1", "SPP1")

test_lasso <- function(expr_file, clin_file, label) {
  expr <- readRDS(expr_file)
  clin <- readRDS(clin_file)
  genes <- intersect(CANDIDATE_GENES, rownames(expr))
  common <- intersect(colnames(expr), clin$sample)
  x <- t(expr[genes, common])
  y <- Surv(clin$OS.time[match(common, clin$sample)], clin$OS[match(common, clin$sample)])
  cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, type.measure = "C")
  for (rule in c("lambda.min", "lambda.1se")) {
    fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = cv_fit[[rule]])
    co <- coef(fit)
    sel <- rownames(co)[co[, 1] != 0]
    cat(sprintf("%-12s %-10s n_selected=%2d  paper6 overlap: %d/6  [%s]\n",
                label, rule, length(sel), length(intersect(sel, PAPER6)),
                paste(intersect(sel, PAPER6), collapse = ",")))
  }
}

test_lasso("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds",
           "data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds", "TCGA")
test_lasso("data/raw/CGGA/cgga693_expr.rds",
           "data/raw/CGGA/cgga693_clinical.rds", "CGGA-693")
