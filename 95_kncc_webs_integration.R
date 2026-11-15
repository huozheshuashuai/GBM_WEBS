# =============================================================================
# 95. Integration analysis: WEBS in the KNCC longitudinal cohort
#     (a) independent replication of WEBS-spliceosome coupling (KNCC T1, T2)
#     (b) protein WEBS longitudinal change T1 -> T2
#     (c) delta-WEBS vs delta-spliceosome coupling, cell-state adjusted
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot; delta <- b$prot_delta

writers <- c("HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT7","KAT8","CREBBP","EP300","KAT14")
erasers <- c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9",
             "HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT5","SIRT6")
w <- intersect(writers, rownames(pm)); e <- intersect(erasers, rownames(pm))
cat("writers quantified:", length(w), "/", length(writers),
    "| erasers:", length(e), "/", length(erasers), "\n")

webs <- colMeans(pm[w, , drop = FALSE], na.rm = TRUE) - colMeans(pm[e, , drop = FALSE], na.rm = TRUE)
cli$WEBS <- webs[cli$sample]

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
sp <- intersect(splice, rownames(pm))
cli$sp_score <- colMeans(pm[sp, , drop = FALSE], na.rm = TRUE)[cli$sample]

sink("results/kncc/integration_webs.txt", split = TRUE)
cat("=== WEBS x KNCC integration ===\n\n")

# (a) cross-sectional replication at T1 and T2
for (tp in c("T1", "T2")) {
  d <- cli[cli$timepoint == tp, ]
  ct <- cor.test(d$WEBS, d$sp_score, method = "spearman")
  cat("KNCC", tp, ": WEBS vs spliceosome protein score rho =", round(ct$estimate, 3),
      " p =", signif(ct$p.value, 3), " (n =", nrow(d), ")\n")
  # mRNA-Layer analog unavailable in KNCC (no RNA); CPTAC reference: protein 0.79 / mRNA -0.03
}

# (b) longitudinal change of protein WEBS
paired <- b$prot_paired
w1 <- cli$WEBS[cli$timepoint == "T1"][match(paired, cli$patient[cli$timepoint == "T1"])]
w2 <- cli$WEBS[cli$timepoint == "T2"][match(paired, cli$patient[cli$timepoint == "T2"])]
wt <- wilcox.test(w2, w1, paired = TRUE)
cat("\nprotein WEBS T1->T2: median", round(median(w1), 3), "->", round(median(w2), 3),
    " paired Wilcoxon p =", signif(wt$p.value, 3), "\n")
d_webs <- w2 - w1
cat("delta WEBS: median", round(median(d_webs), 3), "\n")

# (c) delta-WEBS vs delta-spliceosome, raw and cell-state adjusted
d_sp <- colMeans(delta[sp, , drop = FALSE], na.rm = TRUE)
prolif <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                      "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta))
neuron <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                      "GAP43","TUBB3","MAP2"), rownames(delta))
d_pro <- colMeans(delta[prolif, , drop = FALSE], na.rm = TRUE)
d_neu <- colMeans(delta[neuron, , drop = FALSE], na.rm = TRUE)

names(d_webs) <- paired
d_webs <- d_webs[names(d_sp)]
ct <- cor.test(d_webs, d_sp, method = "spearman")
cat("\ndelta WEBS vs delta spliceosome: rho =", round(ct$estimate, 3),
    " p =", signif(ct$p.value, 3), "\n")
ct2 <- cor.test(d_webs, d_pro, method = "spearman")
cat("delta WEBS vs delta proliferation: rho =", round(ct2$estimate, 3),
    " p =", signif(ct2$p.value, 3), "\n")

# adjusted: does delta WEBS predict delta spliceosome beyond cell state?
fit <- lm(d_sp ~ d_webs + d_pro + d_neu)
cf <- summary(fit)$coefficients
cat("\nmodel: d_spliceosome ~ d_WEBS + d_prolif + d_neuron\n")
cat("  beta_dWEBS =", round(cf[2, 1], 3), " p =", signif(cf[2, 4], 3), "\n")
cat("  beta_prolif =", round(cf[3, 1], 3), " p =", signif(cf[3, 4], 3), "\n")
cat("  beta_neuron =", round(cf[4, 1], 3), " p =", signif(cf[4, 4], 3), "\n")
cat("  R2 =", round(summary(fit)$r.squared, 3), "\n")

# does WEBS explain the small cell-state-independent residual from script 91?
d_imm <- colMeans(delta[intersect(c("CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC",
                                    "CSF1R","TMEM119","CX3CR1"), rownames(delta)), , drop = FALSE], na.rm = TRUE)
fit0 <- lm(d_sp ~ d_pro + d_neu + d_imm)
fit1 <- lm(d_sp ~ d_webs + d_pro + d_neu + d_imm)
an <- anova(fit0, fit1)
cat("\ndoes d_WEBS explain the residual beyond cell state? F-test p =", signif(an$`Pr(>F)`[2], 3), "\n")
sink()

# figures
t1df <- cli[cli$timepoint == "T1", ]
p1 <- ggplot(t1df, aes(WEBS, sp_score)) +
  geom_point(color = "#7B4EA3", alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(x = "Protein-level WEBS", y = "Spliceosome protein score",
       title = sprintf("KNCC primary tumors (independent cohort): rho = %.2f, p = %s",
                       cor(t1df$WEBS, t1df$sp_score, method = "spearman"),
                       signif(cor.test(t1df$WEBS, t1df$sp_score, method = "spearman")$p.value, 2)))
ggsave("figures_kncc/integration_webs_t1.pdf", p1, width = 5.5, height = 5)

df2 <- data.frame(d_webs = d_webs, d_sp = d_sp)
p2 <- ggplot(df2, aes(d_webs, d_sp)) +
  geom_point(color = "#3D9B8C", alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(x = "delta WEBS (T2 - T1)", y = "delta spliceosome score (T2 - T1)",
       title = sprintf("Longitudinal coupling: rho = %.2f, p = %s",
                       cor(d_webs, d_sp, method = "spearman"), signif(ct$p.value, 2)))
ggsave("figures_kncc/integration_delta_coupling.pdf", p2, width = 5.5, height = 5)
cat("Done.\n")
