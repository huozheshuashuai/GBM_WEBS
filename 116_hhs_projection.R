# =============================================================================
# 116. HHS frozen-signature projection (advisor V6 task 2)
#      Freeze KNCC module definitions + PC2 loadings -> project onto HHS
#      3 validation levels: module direction | axis projection | compartment
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

# ---------- KNCC frozen objects ----------
b <- readRDS("results/kncc/dual_omics_bundle.rds")
kdelta <- b$prot_delta   # KNCC delta matrix (proteins x 105)

gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm_all <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm_all) <- sub("\t.*", "", gmt)

keep_mods <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
               "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
               "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
               "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
               "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")
mods_k <- lapply(hm_all[keep_mods], function(g) intersect(g, rownames(kdelta)))
mods_k$NEURONAL_CURATED <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                                       "GAP43","TUBB3","MAP2"), rownames(kdelta))
mods_k$PROLIF_CURATED <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                                     "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(kdelta))
mods_k <- mods_k[sapply(mods_k, length) >= 8]

# KNCC module delta matrix + PC2 loadings (frozen)
Md <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
pc <- prcomp(t(Md), scale. = TRUE)
pc2_load <- pc$rotation[, 2]   # frozen loadings
cat("frozen PC2 loadings (top +/-):\n")
print(round(sort(pc2_load)[c(1:5, (length(pc2_load)-4):length(pc2_load))], 3))

# ---------- HHS projection ----------
h <- readRDS("results/hhs/hhs_analysis_ready.rds")
hd <- h$delta   # proteins x 40 patients

# shared module definitions restricted to genes present in BOTH cohorts
mods_h <- lapply(names(mods_k), function(mn) intersect(mods_k[[mn]], rownames(hd)))
names(mods_h) <- names(mods_k)
usable <- sapply(mods_h, length) >= 5
cat("\nmodules with >=5 shared proteins:", sum(usable), "/", length(mods_h), "\n")
mods_h <- mods_h[usable]
mods_k <- mods_k[names(mods_h)]

Md_h <- t(sapply(mods_h, function(g) colMeans(hd[g, , drop = FALSE], na.rm = TRUE)))
Md_k <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))

# ---- Level 1: module direction consistency ----
l1 <- data.frame(
  module = names(mods_h),
  n_shared = sapply(mods_h, length),
  kncc_med = apply(Md_k, 1, median, na.rm = TRUE),
  hhs_med = apply(Md_h, 1, median, na.rm = TRUE),
  kncc_dir_up = apply(Md_k, 1, function(x) mean(x > 0, na.rm = TRUE)),
  hhs_dir_up = apply(Md_h, 1, function(x) mean(x > 0, na.rm = TRUE)),
  hhs_p = apply(Md_h, 1, function(x) suppressWarnings(wilcox.test(x, mu = 0)$p.value))
)
l1$hhs_q <- p.adjust(l1$hhs_p, "BH")
l1$same_sign <- sign(l1$kncc_med) == sign(l1$hhs_med)
ct_l1 <- cor.test(l1$kncc_med, l1$hhs_med, method = "spearman")

sink("results/hhs/hhs_projection.txt", split = TRUE)
cat("=== HHS frozen-signature projection ===\n\n")
cat("--- Level 1: module-level delta consistency (", nrow(l1), " shared modules) ---\n")
cat("KNCC vs HHS module median deltas: rho =", round(ct_l1$estimate, 3),
    " p =", signif(ct_l1$p.value, 3), "\n")
cat("sign concordance:", round(mean(l1$same_sign) * 100, 1), "%\n\n")
print(l1[order(-l1$kncc_med), c("module","n_shared","kncc_med","hhs_med","hhs_q","same_sign")],
      row.names = FALSE)

# ---- Level 2: frozen PC2 axis projection ----
load2 <- pc2_load[names(mods_h)]
pc2_k <- as.numeric(t(load2) %*% Md_k[names(mods_h), ])
pc2_h <- as.numeric(t(load2) %*% Md_h[names(mods_h), ])
wt_h <- suppressWarnings(wilcox.test(pc2_h, mu = 0))
wt_k <- suppressWarnings(wilcox.test(pc2_k, mu = 0))
cat("\n--- Level 2: frozen KNCC PC2 axis projected onto HHS ---\n")
cat("KNCC PC2 score: median", round(median(pc2_k), 3), " p =", signif(wt_k$p.value, 3), "\n")
cat("HHS PC2 score:  median", round(median(pc2_h), 3), " p =", signif(wt_h$p.value, 3), "\n")
cat("HHS direction: ", round(mean(pc2_h > 0) * 100, 1), "% patients positive\n")

# ---- Level 3: compartment proxies in HHS (T1 cross-sectional) ----
proxies <- list(
  Microglia = c("PTPRC","CD14","CSF1R","C1QA","AIF1","TMEM119","P2RY12","CX3CR3" ,
                "CX3CR1","CD68","CD163"),
  MES_state = c("CD44","VIM","S100A11","S100A10","ANXA1","LGALS1","MT2A","C1S","NAMPT","EMP3"),
  Endothelial = c("PECAM1","VWF","CLDN5","CDH5","KDR","ESAM","FLT1","ENG"),
  Neuron = c("RBFOX3","SNAP25","SYP","SYT1","STMN2","TUBB3","MAP2","NEFL")
)
proxies <- lapply(proxies, function(g) intersect(g, rownames(h$mat)))
pm1 <- h$mat[, h$ptype == "Primary"]
second <- c("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_IL2_STAT5_SIGNALING",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA")
cat("\n--- Level 3: compartment association in HHS primary tumors ---\n")
for (mn in second) {
  g <- intersect(mods_k[[mn]], rownames(pm1))
  sc <- colMeans(pm1[g, , drop = FALSE], na.rm = TRUE)
  line <- sapply(names(proxies), function(pn) {
    pg <- proxies[[pn]]
    if (length(pg) < 3) return(NA)
    pcv <- colMeans(pm1[pg, , drop = FALSE], na.rm = TRUE)
    r <- suppressWarnings(cor(sc, pcv, method = "spearman", use = "complete.obs"))
    paste0(pn, " rho=", round(r, 2))
  })
  cat(mn, ":", paste(line, collapse = " | "), "\n")
}
sink()

# figures
l1$module <- factor(l1$module, levels = l1$module[order(l1$kncc_med)])
p1 <- ggplot(l1, aes(kncc_med, hhs_med)) +
  geom_point(aes(color = same_sign), size = 2.5) +
  ggrepel::geom_text_repel(aes(label = sub("HALLMARK_", "", module)), size = 2.5, max.overlaps = 25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  scale_color_manual(values = c("TRUE" = "#7B4EA3", "FALSE" = "#d62728")) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(x = "KNCC module delta (median)", y = "HHS module delta (median)",
       title = sprintf("Module-level cross-cohort consistency (rho = %.2f, p = %s)",
                       ct_l1$estimate, signif(ct_l1$p.value, 2)))
ggsave("figures_kncc/hhs_module_consistency.pdf", p1, width = 6.5, height = 5.5)

df2 <- data.frame(cohort = rep(c("KNCC (n=105)", "HHS (n=40)"), c(length(pc2_k), length(pc2_h))),
                  pc2 = c(pc2_k, pc2_h))
p2 <- ggplot(df2, aes(cohort, pc2, fill = cohort)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#7B4EA3", "#3D9B8C")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(x = NULL, y = "Frozen KNCC PC2 axis score (projected delta)",
       title = sprintf("PC2 axis projection: KNCC p = %s | HHS p = %s",
                       signif(wt_k$p.value, 2), signif(wt_h$p.value, 2)))
ggsave("figures_kncc/hhs_pc2_projection.pdf", p2, width = 5.5, height = 5)
cat("Done.\n")
