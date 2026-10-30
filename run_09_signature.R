# =============================================================================
# Driver: Six-gene LASSO-Cox signature
# =============================================================================
# Training: TCGA-GBM (per manuscript Methods, n=158; LinkedOmics n=144 here)
# Candidate pool: lactylation / hypoxia / invasion genes
# Validation: CGGA-693, CGGA-325 (+ IDH-wt GBM subsets), GEO cohorts
# Output: results/sixgene_*.csv, figures/fig_signature_*.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(glmnet)
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggplot2)
  library(dplyr)
})
set.seed(42)
source("09_sixgene_signature.R")  # prepare_training_data, run_lasso_cox, etc.

CANDIDATE_GENES <- c(
  "EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP",
  "HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7",
  "HIF1A", "CA9", "VEGFA", "LDHA", "LDHB", "PKM", "SLC16A3", "ENO1", "ENO2", "GAPDH",
  "MMP2", "ITGA5", "SPP1", "FN1"
)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- 1. Training in TCGA-GBM ----
expr <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")

genes <- intersect(CANDIDATE_GENES, rownames(expr))
common <- intersect(colnames(expr), clin$sample)
x <- t(expr[genes, common])
y <- Surv(clin$OS.time[match(common, clin$sample)], clin$OS[match(common, clin$sample)])

cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, type.measure = "C")
fit <- glmnet(x, y, family = "cox", alpha = 1, lambda = cv_fit$lambda.1se)
coefs <- as.matrix(coef(fit))
coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]

cat("=== Selected genes (lambda.1se) ===\n")
print(coefs)
cat("Paper genes: HIF1A, SIRT2, MMP2, ITGA5, HDAC1, SPP1\n")

write.csv(data.frame(gene = rownames(coefs), coefficient = coefs[, 1]),
          "results/sixgene_coefficients.csv", row.names = FALSE)

# lambda path plot
pdf("figures/fig_signature_lasso.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(cv_fit, main = "LASSO-Cox cross-validation")
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE, main = "Coefficient paths")
dev.off()

# ---- 2. Validation across cohorts ----
validate <- function(expr_v, clin_v, cohort) {
  genes_v <- intersect(rownames(coefs), rownames(expr_v))
  rs <- as.vector(t(expr_v[genes_v, , drop = FALSE]) %*% coefs[genes_v, 1])
  names(rs) <- colnames(expr_v)
  clin_v$risk <- rs[match(clin_v$sample, names(rs))]
  clin_v <- clin_v[!is.na(clin_v$risk), ]
  clin_v$risk_group <- factor(ifelse(clin_v$risk > median(clin_v$risk), "High", "Low"),
                              levels = c("Low", "High"))

  cox <- summary(coxph(Surv(OS.time, OS) ~ risk, data = clin_v))
  fit_km <- survfit(Surv(OS.time, OS) ~ risk_group, data = clin_v)
  km_p <- surv_pvalue(fit_km, data = clin_v)$pval

  p <- ggsurvplot(fit_km, data = clin_v, pval = TRUE, title = cohort,
                  palette = c("#1f77b4", "#d62728"),
                  xlab = "Time (months)", ylab = "Overall survival",
                  legend.labs = c("Low risk", "High risk"))
  pdf(paste0("figures/fig_signature_", gsub("[^A-Za-z0-9]", "_", cohort), ".pdf"),
      width = 6, height = 5)
  print(p)
  dev.off()

  data.frame(cohort = cohort, n = nrow(clin_v), events = sum(clin_v$OS),
             HR = cox$conf.int["risk", "exp(coef)"],
             HR_lower = cox$conf.int["risk", "lower .95"],
             HR_upper = cox$conf.int["risk", "upper .95"],
             cox_p = cox$coefficients["risk", "Pr(>|z|)"], km_p = km_p)
}

val_results <- list()
val_results[["TCGA-GBM (training)"]] <- validate(expr, clin, "TCGA-GBM (training)")

for (cg in c("cgga693", "cgga325")) {
  e2 <- readRDS(paste0("data/raw/CGGA/", cg, "_expr.rds"))
  c2 <- readRDS(paste0("data/raw/CGGA/", cg, "_clinical.rds"))
  val_results[[cg]] <- validate(e2, c2, cg)
  sub <- c2$Histology %in% c("GBM", "rGBM") & c2$IDH_mutation == "Wildtype"
  sub_ids <- intersect(c2$sample[sub], colnames(e2))
  if (length(sub_ids) >= 40)
    val_results[[paste0(cg, "_IDHwtGBM")]] <- validate(
      e2[, sub_ids], c2[match(sub_ids, c2$sample), ], paste0(cg, " IDHwt-GBM"))
}

sumdf <- do.call(rbind, val_results)
write.csv(sumdf, "results/sixgene_validation_summary.csv", row.names = FALSE)
print(sumdf)
message("Done.")
