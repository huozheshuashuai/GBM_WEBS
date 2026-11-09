# =============================================================================
# 77. Immune infiltration deconvolution (xCell): TCGA-GBM + CPTAC-GBM
#     LHS vs LLS differences and WEBS correlations
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(xCell); library(dplyr); library(ggplot2); library(tidyr) })

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
webs_bulk <- function(m) colMeans(m[intersect(WRITERS, rownames(m)), , drop = FALSE]) -
                           colMeans(m[intersect(ERASERS, rownames(m)), , drop = FALSE])

# ---------- cohort 1: TCGA-GBM (Xena star_tpm, back-transform to TPM) ----------
tpm <- read.delim(gzfile("E:/KIMI/data/Xena_GDC/TCGA-GBM.star_tpm.tsv.gz"), check.names = FALSE)
ens <- sub("\\..*$", "", tpm[[1]])
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
mmap <- mapIds(org.Hs.eg.db, keys = unique(ens), column = "SYMBOL",
               keytype = "ENSEMBL", multiVals = "first")
sym <- unname(mmap[ens])
keep <- !is.na(sym)
m <- as.matrix(tpm[keep, -1]); rownames(m) <- sym[keep]
m <- m[!duplicated(rownames(m)), ]
m_tpm <- 2^m - 0.001            # Xena star_tpm is log2(TPM+0.001)
m_tpm[m_tpm < 0] <- 0
webs1 <- webs_bulk(m)           # WEBS on the same (log) scale as before
grp1 <- ifelse(webs1 <= quantile(webs1, 0.25), "LLS",
               ifelse(webs1 >= quantile(webs1, 0.75), "LHS", "LDTS"))

xc1 <- xCellAnalysis(m_tpm)
saveRDS(xc1, "results/xcell_tcga_gbm.rds")

# ---------- cohort 2: CPTAC-GBM RNA-seq (log2 FPKM -> back-transform) ----------
rna <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct", check.names = FALSE)
m2 <- as.matrix(sapply(rna[, -1], as.numeric)); rownames(m2) <- rna[[1]]
m2 <- m2[!is.na(rownames(m2)), ]
webs2 <- webs_bulk(m2)
grp2 <- ifelse(webs2 <= quantile(webs2, 0.25), "LLS",
               ifelse(webs2 >= quantile(webs2, 0.75), "LHS", "LDTS"))
m2_fpkm <- 2^m2
xc2 <- xCellAnalysis(m2_fpkm)
saveRDS(xc2, "results/xcell_cptac_gbm.rds")

# ---------- compare LHS vs LLS per cell type ----------
compare <- function(xc, webs, grp, label) {
  sel <- grp %in% c("LLS", "LHS")
  out <- lapply(rownames(xc), function(ct) {
    x <- xc[ct, sel]
    g <- factor(grp[sel], levels = c("LLS", "LHS"))
    wt <- suppressWarnings(wilcox.test(x ~ g))
    ct2 <- suppressWarnings(cor.test(x, webs[sel], method = "spearman"))
    data.frame(cell_type = ct, rho = unname(ct2$estimate), p_rho = ct2$p.value,
               p_group = wt$p.value,
               mean_LLS = mean(x[g == "LLS"]), mean_LHS = mean(x[g == "LHS"]))
  })
  df <- do.call(rbind, out)
  df$q_group <- p.adjust(df$p_group, "BH")
  df$q_rho <- p.adjust(df$p_rho, "BH")
  df$cohort <- label
  df
}
r1 <- compare(xc1, webs1, grp1, "TCGA-GBM")
r2 <- compare(xc2, webs2, grp2, "CPTAC-GBM")
res <- rbind(r1, r2)
write.csv(res, "results/immune_infiltration_xcell.csv", row.names = FALSE)

sink("results/immune_infiltration_summary.txt", split = TRUE)
cat("=== xCell immune infiltration: LHS vs LLS ===\n\n")
for (lab in c("TCGA-GBM", "CPTAC-GBM")) {
  sub <- res[res$cohort == lab, ]
  sig <- sub[sub$q_group < 0.05, ]
  cat(lab, ": significant cell-type differences (FDR<0.05):", nrow(sig), "\n")
  print(sig %>% arrange(q_group) %>% head(12) %>%
          dplyr::select(cell_type, mean_LLS, mean_LHS, p_group, q_group))
  cat("\n")
}
# concordance across cohorts
mrg <- merge(r1[, c("cell_type","rho","p_group")], r2[, c("cell_type","rho","p_group")],
             by = "cell_type", suffixes = c("_tcga", "_cptac"))
ct <- cor.test(mrg$rho_tcga, mrg$rho_cptac, method = "spearman")
cat("WEBS-correlation concordance across cohorts: rho =", round(ct$estimate, 3),
    " p =", signif(ct$p.value, 3), "\n")
both <- mrg[abs(mrg$rho_tcga) > 0.2 & abs(mrg$rho_cptac) > 0.2 &
              sign(mrg$rho_tcga) == sign(mrg$rho_cptac), ]
cat("Cell types with |rho|>0.2 and same direction in BOTH cohorts:\n")
print(both[order(-abs(both$rho_tcga)), ])
sink()
cat(readLines("results/immune_infiltration_summary.txt"), sep = "\n")

# ---------- figures ----------
top <- res %>% filter(cohort == "TCGA-GBM") %>% arrange(q_group) %>% head(12)
top$dir <- ifelse(top$mean_LHS > top$mean_LLS, "Higher in LHS", "Higher in LLS")
p1 <- ggplot(top, aes(reorder(cell_type, -q_group), -log10(p_group), fill = dir)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_fill_manual(values = c("Higher in LHS" = "#7B4EA3", "Higher in LLS" = "#3D9B8C"), name = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "-log10(Wilcoxon p)",
       title = "xCell infiltration differences: LHS vs LLS (TCGA-GBM)")
ggsave("figures_final/fig_immune_xcell.pdf", p1, width = 8, height = 5)

key_cells <- intersect(c("Macrophages", "Macrophages M1", "Macrophages M2", "Microglia",
                         "Monocytes", "CD4+ T-cells", "CD8+ T-cells", "Neutrophils"),
                       rownames(xc1))
df <- lapply(key_cells, function(ct) {
  data.frame(cell_type = ct, score = xc1[ct, grp1 %in% c("LLS","LHS")],
             group = grp1[grp1 %in% c("LLS","LHS")])
}) %>% bind_rows()
df$group <- factor(df$group, levels = c("LLS","LHS"))
p2 <- ggplot(df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.6) + facet_wrap(~cell_type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("LLS" = "#1f77b4", "LHS" = "#d62728")) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "xCell score",
       title = "Key immune populations by WEBS group (TCGA-GBM)")
ggsave("figures_final/fig_immune_boxplots.pdf", p2, width = 10, height = 6)
cat("Done.\n")
