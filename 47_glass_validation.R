# =============================================================================
# 47. GLASS cohort validation: survival + longitudinal primary-recurrence
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(survminer); library(ggplot2); library(dplyr) })

DATA <- "E:/KIMI/data/GLASS_extracted"

# --- expression matrix (gene-level TPM) ---
tpm <- read.delim(file.path(DATA, "gene_tpm_matrix_all_samples.tsv"),
                  check.names = FALSE)
genes <- tpm[[1]]
expr <- as.matrix(tpm[, -1])
rownames(expr) <- genes
expr <- log2(expr + 1)

# --- aliquot -> case mapping (from Synapse biospecimen table w/ rna_barcode) ---
map <- read.csv("E:/KIMI/data/Job-1200313187969021782211716284.csv")
map$rna_dot <- gsub("-", ".", map$rna_barcode)
map <- map[, c("rna_dot", "sample_barcode", "case_barcode", "sample_type")]
clin_case <- read.csv(file.path(DATA, "clinical_cases.csv"))
pairs <- read.csv(file.path(DATA, "analysis_rnaseq_pairs.csv"))
pairs$tumor_a_dot <- gsub("-", ".", pairs$tumor_barcode_a)
pairs$tumor_b_dot <- gsub("-", ".", pairs$tumor_barcode_b)

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

score_samples <- function(m) {
  webs <- colMeans(m[intersect(WRITERS, rownames(m)), , drop = FALSE]) -
          colMeans(m[intersect(ERASERS, rownames(m)), , drop = FALSE])
  g <- intersect(names(PAPER_COEFS), rownames(m))
  risk <- as.vector(t(m[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  data.frame(aliquot = colnames(m), WEBS = webs, risk = risk)
}
sc <- score_samples(expr)
sc <- merge(sc, map, by.x = "aliquot", by.y = "rna_dot")
cat("Mapped", nrow(sc), "of", ncol(expr), "aliquots to cases\n")

sink("results/glass_validation.txt", split = TRUE)

# ================= 1. Survival on primary tumors =================
prim <- sc[sc$sample_type == "TP", ]
# one row per case (first aliquot)
prim <- prim[!duplicated(prim$case_barcode), ]
prim <- merge(prim, clin_case, by = "case_barcode")
prim$OS.time <- as.numeric(prim$case_overall_survival_mo)
prim$OS <- ifelse(prim$case_vital_status == "dead", 1, 0)
prim <- prim[!is.na(prim$OS.time) & prim$OS.time > 0, ]
cat("\n=== GLASS primary tumors: n =", nrow(prim), ", events =", sum(prim$OS), "===\n")

cat("\n-- Six-gene signature (Cox, continuous) --\n")
print(summary(coxph(Surv(OS.time, OS) ~ risk, data = prim))$coefficients)
prim$risk_group <- factor(ifelse(prim$risk > median(prim$risk), "High", "Low"),
                          levels = c("Low", "High"))
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = prim)
print(surv_pvalue(fit, data = prim))

cat("\n-- WEBS (Cox, continuous) --\n")
print(summary(coxph(Surv(OS.time, OS) ~ WEBS, data = prim))$coefficients)

p1 <- ggsurvplot(fit, data = prim, pval = TRUE,
                 title = paste0("GLASS primary (n=", nrow(prim), "): six-gene signature"),
                 palette = c("#1f77b4", "#d62728"),
                 legend.labs = c("Low risk", "High risk"),
                 xlab = "Time (months)", ylab = "Overall survival")$plot
ggsave("figures_final/fig_glass_km.pdf", p1, width = 6, height = 5)

# ================= 2. Longitudinal primary -> recurrence =================
lp <- pairs[pairs$comparison_type == "longitudinal", ]
lp$WEBS_TP <- sc$WEBS[match(lp$tumor_a_dot, sc$aliquot)]
lp$WEBS_R1 <- sc$WEBS[match(lp$tumor_b_dot, sc$aliquot)]
lp$risk_TP <- sc$risk[match(lp$tumor_a_dot, sc$aliquot)]
lp$risk_R1 <- sc$risk[match(lp$tumor_b_dot, sc$aliquot)]
lp <- lp[!is.na(lp$WEBS_TP) & !is.na(lp$WEBS_R1), ]
cat("\n=== Longitudinal pairs (primary -> first recurrence): n =", nrow(lp), "===\n")
wt <- wilcox.test(lp$WEBS_R1, lp$WEBS_TP, paired = TRUE)
tt <- t.test(lp$WEBS_R1, lp$WEBS_TP, paired = TRUE)
cat("WEBS primary mean =", round(mean(lp$WEBS_TP), 4),
    " recurrence mean =", round(mean(lp$WEBS_R1), 4), "\n")
cat("Paired t-test p =", signif(tt$p.value, 3),
    " Wilcoxon p =", signif(wt$p.value, 3), "\n")
cat("Pairs with WEBS increase at recurrence:", sum(lp$WEBS_R1 > lp$WEBS_TP),
    "/", nrow(lp), "\n")
wr <- wilcox.test(lp$risk_R1, lp$risk_TP, paired = TRUE)
cat("Risk score: primary mean =", round(mean(lp$risk_TP), 3),
    " recurrence mean =", round(mean(lp$risk_R1), 3),
    " Wilcoxon p =", signif(wr$p.value, 3), "\n")
sink()

# longitudinal plot
lp$direction <- ifelse(lp$WEBS_R1 > lp$WEBS_TP, "Increased", "Decreased")
plot_df <- data.frame(
  case = rep(lp$case_barcode, 2),
  stage = factor(rep(c("Primary", "Recurrence"), each = nrow(lp)),
                 levels = c("Primary", "Recurrence")),
  WEBS = c(lp$WEBS_TP, lp$WEBS_R1))
p2 <- ggplot(plot_df, aes(stage, WEBS, group = case)) +
  geom_line(alpha = 0.3, color = "grey50") +
  geom_point(aes(color = stage), size = 1.5) +
  scale_color_manual(values = c("Primary" = "#3D9B8C", "Recurrence" = "#7B4EA3")) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black",
               aes(group = 1)) +
  theme_bw(base_size = 12) +
  labs(title = paste0("GLASS longitudinal (n=", nrow(lp), " pairs, Wilcoxon p=",
                      signif(wt$p.value, 2), ")"),
       x = NULL, y = "WEBS") + theme(legend.position = "none")
ggsave("figures_final/fig_glass_longitudinal.pdf", p2, width = 5, height = 5)

write.csv(lp, "results/glass_longitudinal_pairs.csv", row.names = FALSE)
write.csv(prim, "results/glass_primary_scores.csv", row.names = FALSE)
cat("Done.\n")
