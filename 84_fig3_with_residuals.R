# =============================================================================
# 84. Figure 3 upgrade: add RNA-protein residual panels (e, f)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(patchwork) })
PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
webs <- read.csv("results/cptac_webs.csv")
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]

# panels a-d (rebuild from 67_gnd_figures.R logic)
lim <- read.csv("results/cptac_limma_results.csv")
if (!"is_spliceosome" %in% colnames(lim)) {
  gene_col <- colnames(lim)[sapply(lim, function(x) any(x %in% splice))][1]
  lim$is_spliceosome <- lim[[gene_col]] %in% splice
}
lim$cat <- ifelse(lim$is_spliceosome & lim$adj.P.Val < 0.05 & lim$logFC > 0, "Spliceosome up",
           ifelse(lim$is_spliceosome, "Spliceosome", "Other"))
lim <- lim[order(lim$cat != "Spliceosome up"), ]
lim$sz <- ifelse(lim$cat == "Spliceosome up", 2.0, 0.4)
lim$shp <- ifelse(lim$cat == "Spliceosome up", 17, 16)
p3a <- ggplot(lim, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = cat, size = sz, shape = shp), alpha = 0.65) +
  scale_size_identity() + scale_shape_identity() +
  scale_color_manual(values = c("Spliceosome up" = PURPLE, "Spliceosome" = TEAL, "Other" = "grey85")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = "logFC (high- vs low-WEBS)", y = "-log10(p)",
       title = "CPTAC proteome: 107/111 spliceosome up")

cs <- intersect(colnames(pm), webs$sample)
sp <- intersect(splice, rownames(pm))
score <- colMeans(pm[sp, cs, drop = FALSE], na.rm = TRUE)
df3b <- data.frame(WEBS = webs$WEBS[match(cs, webs$sample)], spliceosome = score)
p3b <- ggplot(df3b, aes(WEBS, spliceosome)) +
  geom_point(color = PURPLE, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey30", se = TRUE, formula = y ~ x) +
  theme_bw(base_size = 10) +
  labs(x = "Protein-level WEBS", y = "Spliceosome protein score",
       title = "rho = 0.79 (95% CI 0.68-0.87)")

pg <- read.csv("results/uncoupling_per_gene.csv")
p3c <- ggplot(pg, aes(rho_mrna, rho_protein)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 1.8, alpha = 0.7, color = PURPLE) +
  theme_bw(base_size = 10) +
  labs(x = "mRNA-level rho", y = "Protein-level rho",
       title = "107/109 proteins above diagonal")

td <- data.frame(cutoff = c("20%", "25%", "30%", "median"), up = c(110, 110, 111, 109))
p3d <- ggplot(td, aes(cutoff, up)) +
  geom_bar(stat = "identity", fill = PURPLE) +
  geom_hline(yintercept = 111, linetype = "dashed", color = "grey40") +
  ylim(0, 115) + theme_bw(base_size = 10) +
  labs(x = "High/low-WEBS cutoff", y = "Spliceosome proteins up",
       title = "Threshold robustness (109-111/111)")

# ---- new panels e, f: RNA-protein residuals ----
rna <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct", check.names = FALSE)
rm_ <- as.matrix(sapply(rna[, -1], as.numeric)); rownames(rm_) <- rna[[1]]
genes <- intersect(rownames(pm), rownames(rm_))
samples <- intersect(colnames(pm), intersect(colnames(rm_), webs$sample))
pm2 <- pm[genes, samples]; rm2 <- rm_[genes, samples]
w <- webs$WEBS[match(samples, webs$sample)]
grp <- factor(webs$WEBS_group[match(samples, webs$sample)], levels = c("LLS","LDTS","LHS"))
resid_m <- do.call(rbind, lapply(genes, function(g) resid(lm(pm2[g, ] ~ rm2[g, ]))))
rownames(resid_m) <- genes
splice_resid <- colMeans(resid_m[intersect(splice, genes), , drop = FALSE], na.rm = TRUE)
df3e <- data.frame(group = grp, resid = splice_resid)
kw <- kruskal.test(resid ~ group, data = df3e)
p3e <- ggplot(df3e, aes(group, resid, fill = group)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = c("LLS" = "#1f77b4", "LDTS" = "#ff7f0e", "LHS" = "#d62728")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "Spliceosome residual score",
       title = "mRNA-adjusted residuals (KW p = 1.9e-13)")

beta_df <- read.csv("results/rna_protein_residual_perprotein.csv")
beta_df$grp2 <- ifelse(beta_df$is_spliceosome, "Spliceosome", "Other")
p3f <- ggplot(beta_df, aes(grp2, beta_webs, fill = grp2)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = c("Other" = "grey75", "Spliceosome" = PURPLE)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = "WEBS beta (mRNA-adjusted)",
       title = "103/109 positive after mRNA adj. (OR = 55.9)")

fig3 <- (p3a | p3b | p3c) / (p3d | p3e | p3f) + plot_annotation(tag_levels = "a")
ggsave("figures_GND/Figure_3.pdf", fig3, width = 14, height = 9)
ggsave("figures_GND/Figure_3.png", fig3, width = 14, height = 9, dpi = 200)
cat("Figure 3 upgraded with residual panels.\n")
