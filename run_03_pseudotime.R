# =============================================================================
# Driver: Slingshot pseudotime on malignant cells (paper Figure 5)
# =============================================================================
# Input: data/processed/GSE138794_malignant.rds
# Output: results/pseudotime_summary.txt, figures/fig_pseudotime.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment)
  library(ggplot2); library(dplyr); library(mgcv)
})
set.seed(42)

mal <- readRDS("data/processed/GSE138794_malignant.rds")
cat("Malignant cells:", ncol(mal), "\n")

# Re-embed malignant cells alone for a clean trajectory
mal <- NormalizeData(mal) |> FindVariableFeatures(nfeatures = 2000) |>
  ScaleData() |> RunPCA(verbose = FALSE)
mal <- FindNeighbors(mal, dims = 1:20) |> FindClusters(resolution = 0.5)
mal <- RunUMAP(mal, dims = 1:20)

# WEBS group from global IQR thresholds already stored; use as cluster labels
sce <- as.SingleCellExperiment(mal)
reducedDims(sce)$UMAP <- Embeddings(mal, "umap")
sce$WEBS_group <- mal$WEBS_group

# Start from LLS (low WEBS = putative origin per paper)
sce <- slingshot(sce, clusterLabels = "WEBS_group", reducedDim = "UMAP",
                 start.clus = "LLS", allow.breaks = FALSE)
pt <- slingPseudotime(sce)
mal$pseudotime <- pt[, 1]
mal$lineage <- factor(!is.na(pt[, 1]))

# GAM: WEBS ~ pseudotime
webs_col <- mal$WEBS
if (is.data.frame(webs_col)) webs_col <- webs_col[, 1]
df <- data.frame(WEBS = webs_col, pt = mal$pseudotime)
df <- df[!is.na(df$pt), ]
fit <- gam(WEBS ~ s(pt, bs = "cr"), data = df)
s <- summary(fit)

ct <- cor.test(df$WEBS, df$pt, method = "spearman")

sink("results/pseudotime_summary.txt")
cat("=== Pseudotime (slingshot) on malignant cells ===\n\n")
cat("Cells with pseudotime:", nrow(df), "/", ncol(mal), "\n")
cat("GAM WEBS ~ s(pseudotime): p =", signif(s$s.pv, 3),
    " dev.expl =", round(s$dev.expl, 3), "\n")
cat("Spearman WEBS vs pseudotime: rho =", round(ct$estimate, 3),
    " p =", signif(ct$p.value, 3), "\n")
cat("(Paper: pseudotime connects LLS to LHS; expect positive correlation)\n")
sink()
cat(readLines("results/pseudotime_summary.txt"), sep = "\n")

# Plots
umap_df <- data.frame(UMAP1 = Embeddings(mal, "umap")[, 1],
                      UMAP2 = Embeddings(mal, "umap")[, 2],
                      WEBS = webs_col, group = mal$WEBS_group,
                      pt = mal$pseudotime)
p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = pt)) +
  geom_point(size = 0.3) + scale_color_viridis_c(na.value = "grey80") +
  theme_bw() + ggtitle("Pseudotime (LLS origin)")
p2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = WEBS)) +
  geom_point(size = 0.3) + scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + ggtitle("WEBS")
p3 <- ggplot(df, aes(pt, WEBS)) +
  geom_point(size = 0.3, alpha = 0.3) + geom_smooth(method = "gam", color = "red") +
  theme_bw() + labs(x = "Pseudotime", y = "WEBS",
                    title = sprintf("WEBS ~ pseudotime (rho = %.2f)", ct$estimate))

dir.create("figures", showWarnings = FALSE)
pdf("figures/fig_pseudotime.pdf", width = 14, height = 5)
print(p1 + p2 + p3)
dev.off()

saveRDS(mal, "data/processed/GSE138794_malignant_pt.rds")
message("Done.")
