# =============================================================================
# 35. Paper-style figures (matching original v16 style)
# =============================================================================
# Style: theme_bw, colored violin + white boxplot, on-plot KW p, panel letters,
# original color scheme (LLS light blue / LDTS gold / LHS salmon)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
  library(ggpubr); library(ggrepel)
})
set.seed(42)
dir.create("figures_paper", showWarnings = FALSE)

C_GROUP <- c("LLS" = "#7BAFD4", "LDTS" = "#F0B323", "LHS" = "#E8685E")
C_REGION <- c("Normal-like" = "#6FBF8B", "Core-like" = "#D96C6C", "Margin-like" = "#F0B323")

obj <- readRDS("data/processed/GSE138794_final.rds")
mal <- readRDS("data/processed/GSE138794_malignant_pt.rds")

# helpers -------------------------------------------------------------------
violin_box <- function(df, x, y, fill_vals, ylab = "WEBS", test = c("kw", "wilcox"), title = NULL) {
  test <- match.arg(test)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.9) +
    geom_boxplot(width = 0.13, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = fill_vals) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = ylab)
  if (!is.null(title)) p <- p + ggtitle(title)
  if (test == "kw") {
    kt <- kruskal.test(df[[y]] ~ df[[x]])
    lab <- ifelse(kt$p.value < 1e-300, "p < 1e-300", sprintf("Kruskal-Wallis p = %s", formatC(kt$p.value, format = "e", digits = 1)))
    p <- p + labs(subtitle = lab)
  }
  p
}

# ===========================================================================
# Fig 2: WEBS landscape (original: cell-type UMAP, WEBS UMAP, violin by group)
# ===========================================================================
obj$cc_label <- obj$cell_type
has <- !is.na(obj$original_type)
obj$cc_label[has] <- as.character(obj$original_type[has])
obj$cc_label <- dplyr::recode(obj$cc_label,
  "Neoplastic_cell" = "Malignant", "Myeloid" = "Macrophage",
  "Endothelial_cell" = "Endothelial")
obj$cc_label[!obj$cc_label %in% c("Astrocyte","Endothelial","Macrophage","Malignant",
                                  "Neuron","Oligodendrocyte","Pericyte")] <- "Unassigned"

pA <- DimPlot(obj, group.by = "cc_label", label = TRUE, repel = TRUE, label.size = 3.5) +
  ggtitle("(a)  Cell Type Annotation (Marker-Based)") +
  theme_bw(base_size = 11) +
  theme(legend.text = element_text(size = 8), legend.key.size = unit(0.35, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2.5)))

pB <- FeaturePlot(obj, features = "WEBS", pt.size = 0.3) +
  scale_color_gradientn(colors = c("#3B3B6D", "#9C9C9C", "#8B1E3F"),
                        limits = c(-1, 1), oob = scales::squish) +
  ggtitle("(b)  Writer-Eraser Balance Score (WEBS)") +
  labs(subtitle = "Writers - Erasers") +
  theme_bw(base_size = 11)

md <- obj@meta.data
pC <- violin_box(md, "WEBS_group", "WEBS", C_GROUP,
                 ylab = "Writer - Eraser (WEBS proxy)",
                 title = "(c)  WEBS by group") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

pdf("figures_paper/Fig2_WEBS_landscape.pdf", width = 15, height = 5.5)
print(pA + pB + pC)
dev.off()

# ===========================================================================
# Fig 3: Neftel states (violin per module + KW)
# ===========================================================================
neftel_cols <- grep("^Neftel_", colnames(mal@meta.data), value = TRUE)
long <- mal@meta.data %>%
  select(WEBS_group, all_of(neftel_cols)) %>%
  tidyr::pivot_longer(all_of(neftel_cols), names_to = "module", values_to = "score") %>%
  mutate(module = sub("Neftel_", "", module))
pD <- ggplot(long, aes(x = WEBS_group, y = score, fill = WEBS_group)) +
  geom_violin(scale = "width", alpha = 0.9) +
  geom_boxplot(width = 0.13, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
  facet_wrap(~module, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = C_GROUP) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  ggtitle("(a)  Neftel module scores by WEBS group (malignant cells)") +
  labs(x = "WEBS group", y = "Module score")

pE <- ggplot(mal@meta.data, aes(x = WEBS, y = Neftel_MES)) +
  geom_point(size = 0.15, alpha = 0.25, color = "grey40") +
  geom_smooth(method = "gam", color = "#8B1E3F", se = TRUE) +
  theme_bw(base_size = 11) +
  ggtitle("(b)  WEBS vs Mesenchymal module") +
  labs(x = "WEBS", y = "MES score")
ct_mes <- cor.test(mal$WEBS, mal$Neftel_MES, method = "spearman")
pE <- pE + annotate("text", x = min(mal$WEBS), y = max(mal$Neftel_MES),
                    hjust = 0, vjust = 1, size = 4,
                    label = sprintf("Spearman rho = %.2f, p = %s",
                                    ct_mes$estimate, formatC(ct_mes$p.value, format = "e", digits = 1)))

pdf("figures_paper/Fig3_Neftel_WEBS.pdf", width = 15, height = 5.5)
print(pD + pE + plot_layout(widths = c(2, 1)))
dev.off()

# ===========================================================================
# Fig 4: Region + hypoxia decoupling (matches original 4-panel layout)
# ===========================================================================
md_all <- obj@meta.data
p4a <- violin_box(md_all[!is.na(md_all$Region_Type), ], "Region_Type", "WEBS",
                  C_REGION, ylab = "WEBS",
                  title = "(a)  WEBS across invasion region types")

# patient-paired margin vs core (malignant cells, median per patient-region)
pr <- md_all %>% filter(Region_Type %in% c("Core-like", "Margin-like")) %>%
  group_by(sample, Region_Type) %>%
  summarise(m = median(WEBS), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Region_Type, values_from = m) %>%
  filter(!is.na(`Core-like`), !is.na(`Margin-like`))
wt <- wilcox.test(pr$`Core-like`, pr$`Margin-like`, paired = TRUE)
pr_long <- pr %>% tidyr::pivot_longer(-sample, names_to = "region", values_to = "webs")
p4b <- ggplot(pr_long, aes(x = region, y = webs)) +
  geom_boxplot(aes(fill = region), width = 0.4, outlier.shape = NA, alpha = 0.85) +
  geom_line(aes(group = sample), color = "grey70") +
  geom_point(aes(color = region), size = 2) +
  scale_fill_manual(values = C_REGION) + scale_color_manual(values = C_REGION) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  ggtitle(sprintf("(b)  Patient-level paired Margin vs Core WEBS\npaired Wilcoxon, n=%d, p=%s",
                  nrow(pr), formatC(wt$p.value, format = "e", digits = 1))) +
  labs(x = "Region type", y = "Patient median WEBS")

p4c <- ggplot(md_all, aes(x = CA9_VEGFA, y = WEBS)) +
  geom_point(size = 0.15, alpha = 0.25, color = "grey40") +
  geom_smooth(method = "lm", color = "#E8685E", se = FALSE) +
  theme_bw(base_size = 11) +
  ggtitle("(c)  CA9 + VEGFA vs WEBS") +
  labs(subtitle = sprintf("All cells (n = %s)", format(nrow(md_all), big.mark = ",")),
       x = "CA9 + VEGFA score", y = "WEBS")
ct_cv <- cor.test(md_all$CA9_VEGFA, md_all$WEBS, method = "spearman")
p4c <- p4c + annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 4,
                      label = sprintf("Spearman rho = %.3f", ct_cv$estimate))

p4d <- violin_box(md_all, "WEBS_group", "Hallmark_Hypoxia", C_GROUP,
                  ylab = "Hallmark Hypoxia score",
                  title = "(d)  Hallmark Hypoxia score by WEBS group")

pdf("figures_paper/Fig4_hypoxia_decoupling.pdf", width = 12, height = 10)
print((p4a + p4b) / (p4c + p4d))
dev.off()

# ===========================================================================
# Fig 9: functional characterization (original layout a-f)
# ===========================================================================
# (a) MES violin by group
p9a <- violin_box(mal@meta.data, "WEBS_group", "Neftel_MES", C_GROUP,
                  ylab = "Neftel Mesenchymal Score") +
  ggtitle("(a)  Mesenchymal score by WEBS group")

# (b) AUCell pathway bar (from 32 results)
pw <- read.csv("results/aucell_hallmark_LHS_LLS.csv")
pw_sel <- pw %>% filter(pathway %in% paste0("HALLMARK_", c(
  "DNA_REPAIR","EPITHELIAL_MESENCHYMAL_TRANSITION","INFLAMMATORY_RESPONSE","HYPOXIA",
  "APOPTOSIS","OXIDATIVE_STRESS","GLYCOLYSIS","MYC_TARGETS_V1","MYC_TARGETS_V2",
  "E2F_TARGETS","G2M_CHECKPOINT","MTORC1_SIGNALING","PI3K_AKT_MTOR_SIGNALING",
  "UNFOLDED_PROTEIN_RESPONSE","TNFA_SIGNALING_VIA_NFKB"))) %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway))
pw_sel$dir <- ifelse(pw_sel$delta > 0, "Up in LHS", "Up in LLS")
p9b <- ggplot(pw_sel, aes(x = reorder(pathway, delta), y = delta, fill = dir)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("Up in LHS" = "#E8685E", "Up in LLS" = "#4B6CB7")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  ggtitle("(b)  AUCell pathway activity (LHS - LLS)") +
  labs(x = NULL, y = "AUC Diff (LHS-LLS)")

# (c) CytoTRACE-like stemness (gene-count proxy)
mal$stemness <- scale(mal$nFeature_RNA)[, 1]
p9c <- violin_box(mal@meta.data, "WEBS_group", "stemness", C_GROUP,
                  ylab = "Stemness Score") +
  ggtitle("(c)  Stemness by WEBS group")

# (d) perturbation sensitivity: leave-one-out delta WEBS (from sensitivity csv)
loo <- read.csv("results/sensitivity_leaveout.csv")
loo1 <- loo %>% filter(type == "leave-one-out") %>%
  mutate(gene = sub("minus_", "", perturbation),
         delta = spearman_vs_default - 1,
         type_g = ifelse(gene %in% c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP"),
                         "Writer", "Eraser"))
p9d <- ggplot(loo1, aes(x = reorder(gene, delta), y = delta, fill = type_g)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("Writer" = "#E8685E", "Eraser" = "#4B6CB7")) +
  theme_bw(base_size = 10) +
  ggtitle("(d)  Perturbation sensitivity of WEBS (leave-one-out)") +
  labs(x = NULL, y = "Delta correlation vs default", fill = "Type")

# (e) top differential TFs — full 67-TF set, proper stats
TF67 <- c("OLIG2","SOX2","ATF4","SNAI1","STAT3","FOS","JUN","MYC","E2F1","HIF1A",
          "NFKB1","RELA","SP1","YY1","EGR1","KLF4","CEBPB","FOXM1","TEAD1","ASCL1",
          "ID2","TCF12","NRF1","ZEB1","ZEB2","TWIST1","SNAI2","ETS1","ELK1","JUNB",
          "JUND","FOSB","MXI1","MAX","USF1","USF2","XBP1","ATF6","ATF3","DDIT3",
          "CEBPA","CEBPD","CEBPG","KLF6","KLF9","EGR2","EGR3","FOSB","FOSL1","FOSL2",
          "JUN","JUNB","BACH1","NFE2L2","MAFF","MAFG","MAFK","CREB1","CREM","ATF1",
          "ATF2","ATF5","ATF7","CREB3","CREB3L1","CREB5","HES1","HEY1","SOX4")
TF67 <- unique(TF67)
expr <- GetAssayData(mal, layer = "data")
grp <- mal$WEBS_group
tf_use <- intersect(TF67, rownames(expr))
tf_res <- t(sapply(tf_use, function(g) {
  x <- expr[g, ]
  c(lfc = mean(x[grp == "LHS"]) - mean(x[grp == "LLS"]),
    p = wilcox.test(x[grp == "LHS"], x[grp == "LLS"])$p.value)
}))
tf_df <- data.frame(TF = rownames(tf_res), logFC = tf_res[, 1], p = tf_res[, 2])
tf_df$q <- p.adjust(tf_df$p, method = "BH")
tf_df <- tf_df %>% filter(q < 0.05) %>% arrange(desc(abs(logFC))) %>% head(12) %>%
  arrange(logFC) %>% mutate(TF = factor(TF, levels = TF))
write.csv(tf_df, "results/tf_diff_67_LHS_LLS.csv", row.names = FALSE)
p9e <- ggplot(tf_df, aes(x = TF, y = logFC, fill = logFC > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E8685E", "FALSE" = "#4B6CB7")) +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  ggtitle("(e)  Top differential TFs (FDR < 0.05)") +
  labs(x = NULL, y = "log2 Fold Change (LHS vs LLS)")

pdf("figures_paper/Fig9_functional.pdf", width = 14, height = 10)
print((p9a + p9b + p9c) / (p9d + p9e))
dev.off()

message("Paper-style figures done")
