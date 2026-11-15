# =============================================================================
# 104. Module co-evolution + continuous vs discrete trajectory analysis
#      Per-patient module deltas -> correlation structure, bootstrap stability,
#      PCA of delta vectors, continuous vs discrete model comparison
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta

# predefined modules: Hallmark + curated
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
keep_mods <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
               "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHIMAL_TRANSITION",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_ANGIOGENESIS",
               "HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
               "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
               "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
               "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")
keep_mods <- keep_mods[keep_mods %in% names(hm)]
mods <- lapply(hm[keep_mods], function(g) intersect(g, rownames(delta)))
mods$NEURONAL_CURATED <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                                     "GAP43","TUBB3","MAP2"), rownames(delta))
mods$PROLIF_CURATED <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                                   "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta))
mods <- mods[sapply(mods, length) >= 8]

# per-patient module delta matrix
Md <- t(sapply(mods, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)))
colnames(Md) <- colnames(delta)
Md <- Md[, colSums(is.na(Md)) == 0]
cat("module delta matrix:", nrow(Md), "modules x", ncol(Md), "patients\n")

# ---- 1. module co-evolution: correlation of per-patient deltas ----
cc <- cor(t(Md), method = "spearman")
ct_p <- matrix(NA, nrow(Md), nrow(Md), dimnames = list(rownames(Md), rownames(Md)))
for (i in 1:nrow(Md)) for (j in 1:nrow(Md))
  ct_p[i, j] <- suppressWarnings(cor.test(Md[i, ], Md[j, ], method = "spearman")$p.value)

# bootstrap stability of the strongest correlations
B <- 1000
n <- ncol(Md)
boot_cor <- function(i, j) replicate(B, cor(Md[i, sample(n, n, TRUE)], Mj <- Md[j, sample(n, n, TRUE)], method = "spearman"))
# proper paired bootstrap: resample patients jointly
boot_pair <- function(i, j) {
  replicate(B, { idx <- sample(n, n, TRUE); cor(Md[i, idx], Md[j, idx], method = "spearman") })
}

pairs_grid <- t(combn(rownames(Md), 2))
pcors <- data.frame(m1 = pairs_grid[, 1], m2 = pairs_grid[, 2],
                    rho = mapply(function(a, b) cc[a, b], pairs_grid[, 1], pairs_grid[, 2]),
                    p = mapply(function(a, b) ct_p[a, b], pairs_grid[, 1], pairs_grid[, 2]))
pcors$q <- p.adjust(pcors$p, "BH")
pcors <- pcors[order(-abs(pcors$rho)), ]

sink("results/kncc/module_coevolution.txt", split = TRUE)
cat("=== Module co-evolution (per-patient delta correlations, n =", ncol(Md), ") ===\n\n")
cat("top 15 correlated module pairs:\n")
print(head(pcors, 15), row.names = FALSE)

# bootstrap CI for top 6 pairs
cat("\nbootstrap 95% CI for top pairs:\n")
for (k in 1:min(6, nrow(pcors))) {
  i <- pcors$m1[k]; j <- pcors$m2[k]
  bs <- boot_pair(i, j)
  ci <- quantile(bs, c(0.025, 0.975))
  cat(sprintf("  %s ~ %s: rho = %.2f, boot 95%% CI [%.2f, %.2f]\n", i, j, cc[i, j], ci[1], ci[2]))
}

# ---- 2. continuous vs discrete: PCA of module delta vectors ----
pc <- prcomp(t(Md), scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
cat("\nPCA of module-delta vectors: PC1 =", round(ve[1] * 100, 1), "% PC2 =",
    round(ve[2] * 100, 1), "% PC3 =", round(ve[3] * 100, 1), "%\n")
cat("cumulative PC1-3:", round(sum(ve[1:3]) * 100, 1), "%\n")

# is PC1 just the proliferation axis?
prol_delta <- Md["PROLIF_CURATED", ]
cat("PC1 vs proliferation delta: rho =",
    round(cor(pc$x[, 1], prol_delta, method = "spearman"), 3), "\n")

# discreteness: gap statistic / silhouette for k=2..4 on scaled deltas
set.seed(1)
sil <- function(k) {
  km <- kmeans(t(Md), centers = k, nstart = 50)
  ss <- cluster::silhouette(km$cluster, dist(t(Md)))
  mean(ss[, 3])
}
sils <- sapply(2:4, sil)
cat("\nmean silhouette for k=2..4:", paste(round(sils, 3), collapse = ", "), "\n")
cat("(values >0.25 suggest structure; <0.25 => continuum, no discrete subtypes)\n")
sink()

# figures
library(pheatmap)
pdf("figures_kncc/module_coevolution_heatmap.pdf", width = 8, height = 7)
pheatmap(cc, breaks = seq(-1, 1, length.out = 101),
         color = colorRampPalette(c("#3D9B8C", "white", "#7B4EA3"))(100),
         main = "Module delta co-evolution (Spearman, per-patient)")
dev.off()

pcdf <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], prolif = prol_delta)
p2 <- ggplot(pcdf, aes(PC1, PC2, color = prolif)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "#1f77b4", mid = "grey80", high = "#d62728") +
  theme_bw(base_size = 12) +
  labs(title = "Patient evolution map (module deltas)",
       subtitle = sprintf("PC1 %.0f%% + PC2 %.0f%% variance; color = proliferation delta",
                          ve[1] * 100, ve[2] * 100),
       color = "prolif delta")
ggsave("figures_kncc/module_evolution_pca.pdf", p2, width = 6.5, height = 5.5)
cat("Done.\n")
