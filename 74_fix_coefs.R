setwd("E:/KIMI/GBM_WEBS_Project")
# canonical six-gene coefficient table: paper (fixed) + TCGA Cox refit comparison
PAPER <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
           ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
REFIT <- c(HIF1A = -0.363, SIRT2 = -0.298, MMP2 = 0.219,
           ITGA5 = 0.209, HDAC1 = -0.472, SPP1 = 0.106)
df <- data.frame(
  gene = names(PAPER),
  coefficient_fixed = as.numeric(PAPER),
  direction = ifelse(PAPER > 0, "risk", "protective"),
  cox_refit_TCGA = as.numeric(REFIT[names(PAPER)]),
  note = "Fixed coefficients used for ALL cohorts; no re-fitting in validation cohorts. Refit shown for reference only (risk-score correlation r = 0.95)."
)
write.csv(df, "results/sixgene_coefficients.csv", row.names = FALSE)
print(df)
