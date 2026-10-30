# =============================================================================
# 30. Fix all defective figures
# =============================================================================
# (1) Survival/signature KM: single-page, no blank first page
# (2) CellChat compare: proper LLS/LHS labels
# (3) Acetylome spliceosome heatmap (missing)
# (4) Sensitivity figure: full y-range, fix platform panel
# (5) Spatial WEBS: symmetric color scale around 0
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(survival); library(survminer)
  library(ggplot2); library(dplyr); library(patchwork); library(pheatmap)
})
set.seed(42)
dir.create("figures_fixed", showWarnings = FALSE)

# ---------- (1) KM curves, single page ----------
WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")

km_plot <- function(clin, cohort, outfile) {
  q <- quantile(clin$WEBS, c(0.25, 0.75), na.rm = TRUE)
  clin$grp <- factor(ifelse(clin$WEBS <= q[1], "LLS",
                     ifelse(clin$WEBS >= q[2], "LHS", "LDTS")), levels = c("LLS","LDTS","LHS"))
  fit <- survfit(Surv(OS.time, OS) ~ grp, data = clin)
  p <- ggsurvplot(fit, data = clin, pval = TRUE, title = cohort,
                  palette = c("#1f77b4","#ff7f0e","#d62728"),
                  legend.title = "WEBS", legend.labs = c("LLS","LDTS","LHS"),
                  xlab = "Time (months)", ylab = "Overall survival")$plot
  ggsave(outfile, p, width = 6, height = 5)
}

expr <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
w <- intersect(WRITERS, rownames(expr)); e <- intersect(ERASERS, rownames(expr))
clin$WEBS <- (colMeans(expr[w, ]) - colMeans(expr[e, ]))[clin$sample]
km_plot(clin, "TCGA-GBM", "figures_fixed/fig_survival_TCGA_GBM.pdf")

for (cg in c("cgga693", "cgga325")) {
  e2 <- readRDS(paste0("data/raw/CGGA/", cg, "_expr.rds"))
  c2 <- readRDS(paste0("data/raw/CGGA/", cg, "_clinical.rds"))
  w2 <- intersect(WRITERS, rownames(e2)); ee <- intersect(ERASERS, rownames(e2))
  c2$WEBS <- (colMeans(e2[w2, ]) - colMeans(e2[ee, ]))[c2$sample]
  km_plot(c2, cg, paste0("figures_fixed/fig_survival_", cg, ".pdf"))
  keep <- c2$Histology %in% c("GBM","rGBM") & c2$IDH_mutation == "Wildtype"
  ids <- intersect(c2$sample[keep], colnames(e2))
  c3 <- c2[match(ids, c2$sample), ]
  c3$WEBS <- (colMeans(e2[w2, ids]) - colMeans(e2[ee, ids]))[ids]
  km_plot(c3, paste0(cg, " IDHwt-GBM"), paste0("figures_fixed/fig_survival_", cg, "_IDHwt_GBM.pdf"))
}

# signature KM (paper's six-gene, fixed coefficients)
PAPER_COEFS <- c(HIF1A=-0.284, SIRT2=-0.176, MMP2=0.262, ITGA5=0.230, HDAC1=-0.220, SPP1=0.146)
sig_km <- function(expr, clin, cohort, outfile) {
  g <- intersect(names(PAPER_COEFS), rownames(expr))
  rs <- as.vector(t(expr[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  clin$WEBS <- rs[match(clin$sample, colnames(expr))]  # reuse WEBS slot for plotting
  clin <- clin[!is.na(clin$WEBS), ]
  q <- median(clin$WEBS)
  clin$grp <- factor(ifelse(clin$WEBS > q, "High", "Low"), levels = c("Low","High"))
  fit <- survfit(Surv(OS.time, OS) ~ grp, data = clin)
  p <- ggsurvplot(fit, data = clin, pval = TRUE, title = paste(cohort, "(six-gene signature)"),
                  palette = c("#1f77b4","#d62728"), legend.labs = c("Low risk","High risk"),
                  xlab = "Time (months)", ylab = "Overall survival")$plot
  ggsave(outfile, p, width = 6, height = 5)
}
sig_km(expr, clin, "TCGA-GBM", "figures_fixed/fig_signature_TCGA_GBM.pdf")
e693 <- readRDS("data/raw/CGGA/cgga693_expr.rds"); c693 <- readRDS("data/raw/CGGA/cgga693_clinical.rds")
sig_km(e693, c693, "CGGA-693", "figures_fixed/fig_signature_cgga693.pdf")
e325 <- readRDS("data/raw/CGGA/cgga325_expr.rds"); c325 <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
sig_km(e325, c325, "CGGA-325", "figures_fixed/fig_signature_cgga325.pdf")

# ---------- (2) CellChat compare with proper labels ----------
suppressMessages(library(CellChat))
cc <- readRDS("results/cellchat_objects.rds")
merged <- mergeCellChat(list(LLS = cc$LLS, LHS = cc$LHS), add.names = TRUE)
pdf("figures_fixed/fig_cellchat_compare.pdf", width = 12, height = 6)
print(compareInteractions(merged, show.legend = FALSE, group = c("LLS","LHS"),
                          color.use = c("#1f77b4","#d62728")) +
        ggtitle("Differential interaction number: LLS vs LHS"))
dev.off()
pdf("figures_fixed/fig_cellchat_ranknet.pdf", width = 9, height = 7)
print(rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE,
              color.use = c("#1f77b4","#d62728")))
dev.off()

# ---------- (3) Acetylome spliceosome heatmap ----------
acetyl <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct",
                               row.names = 1, check.names = FALSE))
webs <- read.csv("results/cptac_webs.csv")
de <- read.csv("results/acetyl_sites_limma.csv")
top_sites <- de %>% filter(is_spliceosome, !is.na(adj.P.Val)) %>%
  arrange(adj.P.Val) %>% head(40) %>% pull(site)
grp <- webs$WEBS_group[match(colnames(acetyl), webs$sample)]
keep <- grp %in% c("LLS","LHS")
mat <- acetyl[top_sites, keep]
mat <- mat[complete.cases(mat), ]
ann <- data.frame(WEBS = grp[keep], row.names = colnames(mat))
ann <- ann[order(ann$WEBS), , drop = FALSE]
mat <- mat[, rownames(ann)]
pdf("figures_fixed/fig_acetyl_spliceosome_heatmap.pdf", width = 8, height = 8)
pheatmap(mat, scale = "row", annotation_col = ann, show_colnames = FALSE,
         cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main = "Top 40 spliceosome acetylation sites (ordered by WEBS group)")
dev.off()

# ---------- (4) Sensitivity figure (full range + platform from run_11) ----------
sweep <- read.csv("results/sensitivity_maxrank.csv")
loo <- read.csv("results/sensitivity_leaveout.csv")
p1 <- ggplot(sweep, aes(x = maxRank, y = spearman_vs_default)) +
  geom_point(size = 3) + geom_line() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = sweep$maxRank) +
  coord_cartesian(ylim = c(0.6, 1.0)) +
  theme_bw(base_size = 12) +
  labs(x = "UCell maxRank (fraction of genes)", y = "Spearman rho vs default WEBS",
       title = "(a) WEBS sensitivity to maxRank (fixed at 1500 = 0.075)")
p2 <- ggplot(loo, aes(x = spearman_vs_default, fill = type)) +
  geom_histogram(bins = 25, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(x = "Spearman rho vs default WEBS", y = "Count",
       title = "(b) Leave-one/two-out robustness (n = 57)")
p3 <- ggplot() + theme_bw(base_size = 12) +
  annotate("text", x = 0.5, y = 0.5, size = 4.5,
           label = "(c) Cross-platform (GSE131928):\npatient-level WEBS rho = 0.67\n(Smart-seq2 vs 10x, n = 10)") +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())
pdf("figures_fixed/figS_sensitivity.pdf", width = 15, height = 5)
print(p1 + p2 + p3)
dev.off()

# ---------- (5) Spatial with symmetric scale ----------
for (f in list.files("data/processed", pattern = "^spatial_.*\\.rds$", full.names = TRUE)) {
  sid <- sub("^spatial_(.*)\\.rds$", "\\1", basename(f))
  if (grepl("rctd", sid)) next
  sobj <- readRDS(f)
  lim <- max(abs(sobj$WEBS), na.rm = TRUE)
  p <- SpatialFeaturePlot(sobj, features = "WEBS", pt.size.factor = 1.8) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim)) +
    ggtitle(paste(sid, "WEBS"))
  ggsave(paste0("figures_fixed/fig_spatial_", sid, ".pdf"), p, width = 6, height = 6)
}

message("All fixed figures saved to figures_fixed/")
