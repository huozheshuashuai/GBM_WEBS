# =============================================================================
# 32. Fig 9: pathway (AUCell) + TF differences LHS vs LLS in malignant cells
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(AUCell); library(ggplot2); library(dplyr)
  library(patchwork); library(pheatmap)
})
set.seed(42)
dir.create("figures_v2", showWarnings = FALSE)

mal <- readRDS("data/processed/GSE138794_malignant_pt.rds")
expr <- GetAssayData(mal, layer = "data")

# Hallmark gene sets (download once, cache)
gmt_cache <- "data/raw/msigdb_hallmark.gmt"
if (!file.exists(gmt_cache)) {
  url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
  download.file(url, gmt_cache, mode = "wb")
}
hallmark <- fgsea::gmtPathways(gmt_cache)

# subsample malignant cells for AUCell speed
cells <- sample(colnames(expr), min(6000, ncol(expr)))
expr_sub <- expr[, cells]
grp <- mal$WEBS_group[cells]

message("AUCell ranking...")
rk <- AUCell_buildRankings(as.matrix(expr_sub), plotStats = FALSE)
auc <- AUCell_calcAUC(hallmark, rk)
auc_m <- getAUC(auc)

lhs <- which(grp == "LHS"); lls <- which(grp == "LLS")
pw_res <- t(apply(auc_m, 1, function(x) {
  wt <- wilcox.test(x[lhs], x[lls])
  c(lhs = mean(x[lhs]), lls = mean(x[lls]), p = wt$p.value)
}))
pw_df <- data.frame(pathway = rownames(pw_res), pw_res, row.names = NULL)
pw_df$q <- p.adjust(pw_df$p, method = "BH")
pw_df$delta <- pw_df$lhs - pw_df$lls
write.csv(pw_df, "results/aucell_hallmark_LHS_LLS.csv", row.names = FALSE)

top <- pw_df %>% arrange(desc(abs(delta))) %>% head(20) %>%
  mutate(dir = ifelse(delta > 0, "Up in LHS", "Up in LLS"),
         pathway = gsub("HALLMARK_", "", pathway))
pA <- ggplot(top, aes(x = reorder(pathway, delta), y = delta, fill = dir)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("Up in LHS" = "#d62728", "Up in LLS" = "#1f77b4")) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "AUCell score difference (LHS - LLS)",
       title = "A. Hallmark pathways: LHS vs LLS", fill = NULL)

# TF differential expression (malignant)
tfs <- intersect(c("OLIG2","SOX2","EP300","ATF4","SNAI1","STAT3","FOS","JUN","MYC","E2F1",
                   "HIF1A","NFKB1","RELA","SP1","YY1","EGR1","KLF4","CEBPB","FOXM1","TEAD1"),
                 rownames(expr))
tf_res <- t(sapply(tfs, function(g) {
  x <- expr[g, ]
  wt <- wilcox.test(x[grp == "LHS"], x[grp == "LLS"])
  c(lfc = mean(x[grp == "LHS"]) - mean(x[grp == "LLS"]), p = wt$p.value)
}))
tf_df <- data.frame(TF = rownames(tf_res), logFC = tf_res[, "lfc"], p = tf_res[, "p"])
tf_df$q <- p.adjust(tf_df$p, method = "BH")
write.csv(tf_df, "results/tf_diff_LHS_LLS.csv", row.names = FALSE)

pB <- ggplot(tf_df, aes(x = logFC, y = -log10(q), label = TF)) +
  geom_point(aes(color = abs(logFC) > 0.25 & q < 0.05), size = 2.5) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d62728")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "B. TF differential expression (LHS vs LLS)",
       x = "Mean expression difference", y = "-log10(FDR)")

pdf("figures_v2/Fig9_pathway_TF.pdf", width = 13, height = 6)
print(pA | pB)
dev.off()
message("Fig9 done")
