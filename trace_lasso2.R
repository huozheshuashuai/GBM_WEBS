# Trace 2: does multivariable Cox refit of the 6 genes on TCGA reproduce paper coefficients?
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(survival))
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

expr <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
common <- intersect(colnames(expr), clin$sample)
clin <- clin[match(common, clin$sample), ]
for (g in names(PAPER_COEFS)) clin[[g]] <- as.numeric(expr[g, common])

f <- as.formula(paste("Surv(OS.time, OS) ~", paste(names(PAPER_COEFS), collapse = "+")))
fit <- coxph(f, data = clin)
co <- coef(fit)
cat("gene        paper    refit\n")
for (g in names(PAPER_COEFS)) cat(sprintf("%-10s %7.3f  %7.3f\n", g, PAPER_COEFS[g], co[g]))

# also test: risk score correlation between paper coefs and refit coefs
rs_paper <- as.vector(t(expr[names(PAPER_COEFS), common]) %*% PAPER_COEFS)
rs_refit <- as.vector(t(expr[names(PAPER_COEFS), common]) %*% co[names(PAPER_COEFS)])
cat("\nrisk score correlation (paper vs refit):", cor(rs_paper, rs_refit), "\n")
