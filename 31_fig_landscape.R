# =============================================================================
# 31. Main figures: WEBS landscape (Fig2), Neftel modules (Fig3), hypoxia (Fig4)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork); library(UCell)
})
dir.create("figures_v2", showWarnings = FALSE)
obj <- readRDS("data/processed/GSE138794_final.rds")

COL_CT <- c("Malignant"="#d62728","Macrophage_Microglia"="#ff7f0e","Oligodendrocyte"="#2ca02c",
            "Astrocyte"="#9467bd","Neuron"="#8c564b","Endothelial"="#e377c2","Pericyte"="#7f7f7f")
COL_WEBS <- c("LLS"="#1f77b4","LDTS"="#ff7f0e","LHS"="#d62728")

# unified label
obj$cc_label <- obj$cell_type
has <- !is.na(obj$original_type)
obj$cc_label[has] <- as.character(obj$original_type[has])
obj$cc_label[obj$cc_label == "Neoplastic_cell"] <- "Malignant"
obj$cc_label[obj$cc_label == "Myeloid"] <- "Macrophage_Microglia"
obj$cc_label[obj$cc_label == "Endothelial_cell"] <- "Endothelial"
obj$cc_label[!obj$cc_label %in% names(COL_CT)] <- NA

# ---- Fig 2: landscape ----
pA <- DimPlot(obj, group.by = "cc_label", cols = COL_CT, label = TRUE, repel = TRUE) +
  ggtitle("A. Cell types") + theme(legend.position = "none")
pB <- FeaturePlot(obj, features = "WEBS") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  ggtitle("B. WEBS")
pC <- VlnPlot(obj, features = "WEBS", group.by = "cc_label", cols = COL_CT, pt.size = 0) +
  ggtitle("C. WEBS by cell type") + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                          legend.position = "none")
pD <- DimPlot(obj, group.by = "WEBS_group", cols = COL_WEBS) +
  ggtitle("D. WEBS groups (IQR)")
pdf("figures_v2/Fig2_WEBS_landscape.pdf", width = 13, height = 11)
print((pA | pB) / (pC | pD))
dev.off()

# ---- Fig 3: Neftel modules vs WEBS (malignant only) ----
mal <- readRDS("data/processed/GSE138794_malignant_pt.rds")
neftel_cols <- grep("^Neftel_", colnames(mal@meta.data), value = TRUE)
long <- mal@meta.data %>%
  select(WEBS, WEBS_group, all_of(neftel_cols)) %>%
  tidyr::pivot_longer(cols = all_of(neftel_cols), names_to = "module", values_to = "score") %>%
  mutate(module = sub("Neftel_", "", module))
pE <- ggplot(long, aes(x = WEBS_group, y = score, fill = WEBS_group)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.1, outlier.size = 0) +
  facet_wrap(~module, scales = "free_y") +
  scale_fill_manual(values = COL_WEBS) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "A. Neftel modules by WEBS group (malignant cells)", x = NULL, y = "Module score")
pF <- ggplot(mal@meta.data, aes(x = WEBS, y = Neftel_MES)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth(method = "gam", color = "#8c564b") +
  theme_bw(base_size = 12) +
  labs(title = "B. WEBS vs mesenchymal module", x = "WEBS", y = "MES score")
pdf("figures_v2/Fig3_Neftel_WEBS.pdf", width = 13, height = 6)
print(pE | pF)
dev.off()

# ---- Fig 4: hypoxia decoupling ----
md <- obj@meta.data
md$is_mal <- (!is.na(obj$original_type) & obj$original_type == "Neoplastic_cell") |
             (is.na(obj$original_type) & obj$cell_type == "Malignant")
md_mal <- md[md$is_mal, ]
pG <- ggplot(md_mal, aes(x = Hallmark_Hypoxia, y = WEBS)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 12) +
  labs(title = "A. WEBS vs hypoxia (malignant cells, rho = -0.05)",
       x = "Hallmark hypoxia score", y = "WEBS")
pH <- ggplot(md_mal[!is.na(md_mal$Region_Type), ],
             aes(x = Region_Type, y = WEBS, fill = Region_Type)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("Core-like"="#d62728","Margin-like"="#ff7f0e","Normal-like"="#2ca02c")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "B. WEBS by region (malignant cells)", x = NULL, y = "WEBS")
pdf("figures_v2/Fig4_hypoxia_decoupling.pdf", width = 12, height = 5.5)
print(pG | pH)
dev.off()

message("Fig2-4 done")
