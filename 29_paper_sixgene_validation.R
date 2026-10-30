# Validate the paper's exact six-gene model (fixed coefficients) across cohorts
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(survminer) })
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

valid <- function(expr, clin, cohort) {
  g <- intersect(names(PAPER_COEFS), rownames(expr))
  rs <- as.vector(t(expr[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  clin$risk <- rs[match(clin$sample, colnames(expr))]
  clin <- clin[!is.na(clin$risk), ]
  clin$grp <- ifelse(clin$risk > median(clin$risk), "High", "Low")
  cox <- summary(coxph(Surv(OS.time, OS) ~ risk, data = clin))
  km <- surv_pvalue(survfit(Surv(OS.time, OS) ~ grp, data = clin), data = clin)$pval
  cat(cohort, ": n=", nrow(clin),
      " HR=", round(cox$conf.int["risk", "exp(coef)"], 2),
      " cox_p=", signif(cox$coefficients["risk", "Pr(>|z|)"], 3),
      " KM_p=", signif(km, 3), "\n", sep = "")
}

e <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
c <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
valid(e, c, "TCGA-GBM")

e693 <- readRDS("data/raw/CGGA/cgga693_expr.rds")
c693 <- readRDS("data/raw/CGGA/cgga693_clinical.rds")
valid(e693, c693, "CGGA-693")

e325 <- readRDS("data/raw/CGGA/cgga325_expr.rds")
c325 <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
valid(e325, c325, "CGGA-325")

keep <- c693$Histology %in% c("GBM", "rGBM") & c693$IDH_mutation == "Wildtype"
ids <- intersect(c693$sample[keep], colnames(e693))
valid(e693[, ids], c693[match(ids, c693$sample), ], "CGGA-693-IDHwtGBM")

keep325 <- c325$Histology %in% c("GBM", "rGBM") & c325$IDH_mutation == "Wildtype"
ids325 <- intersect(c325$sample[keep325], colnames(e325))
valid(e325[, ids325], c325[match(ids325, c325$sample), ], "CGGA-325-IDHwtGBM")
