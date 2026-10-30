# Generate missing supplementary figures: benchmarking, MOFA2, spatial autocorr
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
dir.create("figures_supp_extra", showWarnings = FALSE)

# ---- S3a: method benchmarking correlation heatmap ----
bc <- read.csv("results/benchmarking_correlation.csv", row.names = 1, check.names = FALSE)
bc_long <- as.data.frame(as.table(as.matrix(bc)))
colnames(bc_long) <- c("Method1", "Method2", "rho")
p1 <- ggplot(bc_long, aes(Method1, Method2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(rho, 2)), size = 3.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "(a) Method correlation (patient-level, Spearman)",
       x = NULL, y = NULL, fill = "rho")

# ---- S3b: computational performance (already have figS_performance) ----
# keep existing figS_performance.pdf

# ---- S9: MOFA2 ----
fc <- read.csv("results/mofa2_factor_webs_corr.csv")
factors_df <- read.csv("results/mofa2_factors.csv")
p2a <- ggplot(fc, aes(x = factor(factor), y = rho)) +
  geom_bar(stat = "identity", aes(fill = q < 0.05)) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey70")) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  geom_hline(yintercept = 0) +
  labs(title = "(a) MOFA2 factor-WEBS correlations", x = "Factor", y = "Spearman rho")
p2b <- ggplot(factors_df, aes(x = WEBS, y = Factor2)) +
  geom_point(size = 2.5, color = "grey40") +
  geom_smooth(method = "lm", color = "#B2182B") +
  theme_bw(base_size = 11) +
  labs(title = "(b) Factor 2 vs WEBS (rho = -0.64, p < 1e-12)",
       x = "Protein-level WEBS", y = "MOFA2 Factor 2")
pdf("figures_supp_extra/figS_mofa2.pdf", width = 11, height = 5)
print(p2a | p2b)
dev.off()

# ---- S10: spatial autocorrelation ----
sa <- read.csv("results/spatial_autocorrelation.csv")
p3 <- ggplot(sa, aes(x = slice, y = moran_I, fill = moran_p_permutation < 0.05)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("p=%.3f", moran_p_permutation)), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey70")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Spatial autocorrelation of WEBS (Moran's I, 1,000 permutations)",
       x = NULL, y = "Moran's I")
pdf("figures_supp_extra/figS_spatial_autocorr.pdf", width = 7, height = 5)
print(p3)
dev.off()

# ---- S3 benchmarking figure ----
pdf("figures_supp_extra/figS_benchmarking.pdf", width = 6.5, height = 5.5)
print(p1)
dev.off()

message("Extra supplementary figures done")
