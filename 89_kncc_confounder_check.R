# =============================================================================
# 89. KNCC module 11.8: competing-hypothesis check
#     Is the spliceosome decline tumor-intrinsic or explained by
#     proliferation drop / cell-composition shift at recurrence?
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

obj <- readRDS("results/kncc/proteome_analysis_ready.rds")
pm <- obj$pm; cli <- obj$annot; delta <- obj$delta; paired_pts <- obj$paired

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
sp <- intersect(splice, rownames(pm))

# ---- protein-marker scores (from the same matrix) ----
prolif <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                      "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(pm))
neuron <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                      "GAP43","TUBB3","MAP2"), rownames(pm))
immune <- intersect(c("CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC",
                      "CSF1R","TMEM119","CX3CR1"), rownames(pm))
cat("markers found: prolif", length(prolif), "neuron", length(neuron), "immune", length(immune), "\n")

score_delta <- function(genes) colMeans(delta[genes, , drop = FALSE], na.rm = TRUE)
d_sp    <- colMeans(delta[sp, , drop = FALSE], na.rm = TRUE)
d_pro   <- score_delta(prolif)
d_neu   <- score_delta(neuron)
d_imm   <- score_delta(immune)

# ---- correlations between score deltas across paired patients ----
ct1 <- cor.test(d_sp, d_pro, method = "spearman")
ct2 <- cor.test(d_sp, d_neu, method = "spearman")
ct3 <- cor.test(d_sp, d_imm, method = "spearman")
cat("\ndelta correlations (per-patient):\n")
cat("  spliceosome ~ proliferation delta: rho =", round(ct1$estimate, 3), " p =", signif(ct1$p.value, 3), "\n")
cat("  spliceosome ~ neuronal delta:      rho =", round(ct2$estimate, 3), " p =", signif(ct2$p.value, 3), "\n")
cat("  spliceosome ~ immune delta:        rho =", round(ct3$estimate, 3), " p =", signif(ct3$p.value, 3), "\n")

# ---- partial: does spliceosome decline survive proliferation adjustment? ----
# regress sp delta on proliferation delta; test residual against 0
fit <- lm(d_sp ~ d_pro)
resid_sp <- resid(fit)
wt_res <- suppressWarnings(wilcox.test(resid_sp, mu = 0))
cat("\nspliceosome delta after regression on proliferation delta:\n")
cat("  median raw delta:", round(median(d_sp), 3),
    " -> median residual:", round(median(resid_sp), 3),
    " (Wilcoxon p =", signif(wt_res$p.value, 3), ")\n")
cat("  variance explained by proliferation: R2 =", round(summary(fit)$r.squared, 3), "\n")

fit2 <- lm(d_sp ~ d_pro + d_neu + d_imm)
cat("  full model R2 =", round(summary(fit2)$r.squared, 3), "\n")
resid2 <- resid(fit2)
wt_res2 <- suppressWarnings(wilcox.test(resid2, mu = 0))
cat("  residual after prolif+neuronal+immune: median", round(median(resid2), 3),
    " p =", signif(wt_res2$p.value, 3), "\n")

# ---- proliferation decline itself ----
wt_pro <- suppressWarnings(wilcox.test(d_pro, mu = 0))
wt_neu <- suppressWarnings(wilcox.test(d_neu, mu = 0))
cat("\nproliferation score delta: median", round(median(d_pro), 3), " p =", signif(wt_pro$p.value, 3), "\n")
cat("neuronal score delta:      median", round(median(d_neu), 3), " p =", signif(wt_neu$p.value, 3), "\n")

# ---- per-protein: spliceosome deltas vs proliferation-linked proteins ----
prot_delta <- apply(delta, 1, median, na.rm = TRUE)
cat("\nmedian delta: spliceosome", round(median(prot_delta[sp]), 3),
    "| proliferation", round(median(prot_delta[prolif]), 3),
    "| genome", round(median(prot_delta[!names(prot_delta) %in% sp]), 3), "\n")

sink("results/kncc/module11_confounder.txt", split = TRUE)
cat("=== Module 11.8: competing hypothesis check ===\n")
cat("spliceosome delta ~ proliferation delta: rho =", round(ct1$estimate, 3), " p =", signif(ct1$p.value, 3), "\n")
cat("spliceosome delta ~ neuronal delta:      rho =", round(ct2$estimate, 3), " p =", signif(ct2$p.value, 3), "\n")
cat("spliceosome delta ~ immune delta:        rho =", round(ct3$estimate, 3), " p =", signif(ct3$p.value, 3), "\n")
cat("proliferation R2 on spliceosome delta:", round(summary(fit)$r.squared, 3), "\n")
cat("residual (prolif-adjusted) median:", round(median(resid_sp), 3), " p =", signif(wt_res$p.value, 3), "\n")
cat("full model R2:", round(summary(fit2)$r.squared, 3),
    " residual median:", round(median(resid2), 3), " p =", signif(wt_res2$p.value, 3), "\n")
sink()

# figure: scatter sp delta vs prolif delta
df <- data.frame(sp = d_sp, pro = d_pro)
p <- ggplot(df, aes(pro, sp)) +
  geom_point(color = "#7B4EA3", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(x = "Proliferation score delta (T2 - T1)", y = "Spliceosome score delta (T2 - T1)",
       title = sprintf("rho = %.2f, p = %s; R2 = %.2f",
                       ct1$estimate, signif(ct1$p.value, 2), summary(fit)$r.squared))
ggsave("figures_kncc/m11_sp_vs_prolif_delta.pdf", p, width = 5.5, height = 5)
cat("Done.\n")
