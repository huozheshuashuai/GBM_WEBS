# Ablation study: leave-one-out C-index of the six-gene signature
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival) })

PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

expr <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")

cindex_for <- function(genes) {
  g <- intersect(genes, rownames(expr))
  rs <- as.vector(t(expr[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  clin$risk <- rs[match(clin$sample, colnames(expr))]
  d <- clin[!is.na(clin$risk), ]
  fit <- coxph(Surv(OS.time, OS) ~ risk, data = d)
  summary(fit)$concordance[1]
}

full <- cindex_for(names(PAPER_COEFS))
rows <- list()
rows[["Full six-gene model"]] <- data.frame(model = "Full six-gene", cindex = full,
                                            delta = 0)
for (g in names(PAPER_COEFS)) {
  ci <- cindex_for(setdiff(names(PAPER_COEFS), g))
  rows[[g]] <- data.frame(model = paste0("Minus ", g), cindex = ci,
                          delta = ci - full)
}
for (g in names(PAPER_COEFS)) {
  ci <- cindex_for(g)
  rows[[paste0(g, " only")]] <- data.frame(model = paste0(g, " only"),
                                           cindex = ci, delta = ci - full)
}

abl <- do.call(rbind, rows)
abl <- abl[order(abl$delta), ]
write.csv(abl, "results/ablation_sixgene.csv", row.names = FALSE)

sink("results/ablation_summary.txt", split = TRUE)
cat("=== Ablation: six-gene signature (TCGA-GBM, n = 144) ===\n\n")
cat("Full model C-index:", round(full, 3), "\n\n")
print(abl, row.names = FALSE)
sink()
