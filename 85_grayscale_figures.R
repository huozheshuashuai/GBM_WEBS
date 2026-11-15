# =============================================================================
# 85. Grayscale (BJC print-safe) versions of main figures
#     F2B and F4a need luminance-monotone regeneration; others convert directly
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(ggplot2); library(dplyr); library(patchwork) })

dir.create("figures_BW", showWarnings = FALSE)

obj <- readRDS("data/processed/GSE138794_final.rds")
COL_CT <- c("Malignant" = "grey15", "Macrophage_Microglia" = "grey45", "Oligodendrocyte" = "grey70",
            "Astrocyte" = "grey55", "Neuron" = "grey35", "Endothelial" = "grey80", "Pericyte" = "grey60")
obj$cc_label <- obj$cell_type
has <- !is.na(obj$original_type)
obj$cc_label[has] <- as.character(obj$original_type[has])
obj$cc_label[obj$cc_label == "Neoplastic_cell"] <- "Malignant"
obj$cc_label[obj$cc_label == "Myeloid"] <- "Macrophage_Microglia"
obj$cc_label[obj$cc_label == "Endothelial_cell"] <- "Endothelial"

pA <- DimPlot(obj, group.by = "cc_label", cols = COL_CT, label = TRUE, repel = TRUE) +
  ggtitle("A. Cell types") + theme_bw(base_size = 10)

# F2B grayscale: luminance-monotone (light = low WEBS, dark = high WEBS)
um <- Embeddings(obj, "umap")
dfB <- data.frame(umap_1 = um[, 1], umap_2 = um[, 2], WEBS = obj$WEBS)
pB <- ggplot(dfB, aes(umap_1, umap_2, color = WEBS)) +
  geom_point(size = 0.3) +
  scale_color_gradient(low = "grey92", high = "grey5", name = "WEBS") +
  ggtitle("B. WEBS") + theme_bw(base_size = 10) +
  theme(legend.position = "right")

# F2C violins with grey fills
pC <- VlnPlot(obj, features = "WEBS", group.by = "cc_label", pt.size = 0) +
  scale_fill_grey(start = 0.85, end = 0.25) +
  ggtitle("C. WEBS by cell type") + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# F2D: shapes already encode groups; use grey fills too
pD <- DimPlot(obj, group.by = "WEBS_group", shape.by = "WEBS_group") +
  scale_color_grey(start = 0.8, end = 0.1, name = "WEBS group") +
  scale_shape_manual(values = c("LLS" = 16, "LDTS" = 17, "LHS" = 15), name = "WEBS group") +
  ggtitle("D. WEBS groups (IQR)") + theme_bw(base_size = 10)

fig2_bw <- (pA | pB) / (pC | pD)
ggsave("figures_BW/Figure_2.pdf", fig2_bw, width = 12, height = 10)
ggsave("figures_BW/Figure_2.png", fig2_bw, width = 12, height = 10, dpi = 200)

# ---- F4a heatmap in luminance-monotone grey (acetylation: white=low, black=high) ----
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
webs <- read.csv("results/cptac_webs.csv")
tt <- read.csv("results/acetyl_sites_limma.csv")
tt$gene <- sub("-.*$", "", tt$site)
top <- tt %>% filter(gene %in% splice, adj.P.Val < 0.05) %>% arrange(desc(abs(t))) %>% head(30)
ac <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct", check.names = FALSE)
sub <- ac[match(top$site, ac[[1]]), ]
am <- as.matrix(sapply(sub[, -1], as.numeric)); rownames(am) <- top$site
common <- intersect(colnames(am), webs$sample)
am <- am[, common]; w2 <- webs[match(common, webs$sample), ]
am <- am[, order(w2$WEBS)]
grp_ord <- w2$WEBS_group[order(w2$WEBS)]
suppressMessages(library(tidyr))
long <- pivot_longer(as.data.frame(am) %>% mutate(site = rownames(am)), -site,
                     names_to = "sample", values_to = "val")
long$site <- factor(long$site, levels = rev(top$site))
long$sample <- factor(long$sample, levels = colnames(am))
ann <- data.frame(sample = colnames(am), group = grp_ord)
ann$sample <- factor(ann$sample, levels = colnames(am))
p4a_main <- ggplot(long, aes(sample, site, fill = val)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "grey5", na.value = "grey85", name = "Acetyl.") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6.5)) +
  labs(x = NULL, y = NULL, title = "Spliceosome acetylation")
p_ann <- ggplot(ann, aes(sample, 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = c("LLS" = "grey85", "LDTS" = "grey50", "LHS" = "grey10"), name = "WEBS group") +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))
p4a <- (p4a_main / p_ann + plot_layout(heights = c(20, 1)))

# phospho volcano (grey): sig spliceosome = black triangles, others grey
ph <- read.csv("results/phospho_limma_LHS_vs_LLS.csv")
ph$cat <- ifelse(ph$adj.P.Val < 0.05 & ph$is_spliceosome & ph$logFC > 0, "Spliceosome sig (up)",
          ifelse(ph$adj.P.Val < 0.05 & ph$is_spliceosome, "Spliceosome sig (down)",
          ifelse(ph$adj.P.Val < 0.05, "Other sig",
          ifelse(ph$is_spliceosome, "Spliceosome NS", "NS"))))
ph <- ph[order(!grepl("Spliceosome sig", ph$cat)), ]
ph$sz <- ifelse(grepl("Spliceosome sig", ph$cat), 1.8, 0.4)
ph$shp <- ifelse(grepl("Spliceosome sig", ph$cat), 17, 16)
ph$gcol <- ifelse(ph$cat == "Spliceosome sig (up)", "grey10",
           ifelse(ph$cat == "Spliceosome sig (down)", "grey40",
           ifelse(ph$cat == "Other sig", "grey70",
           ifelse(ph$cat == "Spliceosome NS", "grey50", "grey90"))))
p4b <- ggplot(ph, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = gcol, size = sz, shape = shp), alpha = 0.55) +
  scale_size_identity() + scale_shape_identity() + scale_color_identity() +
  theme_bw(base_size = 10) +
  labs(x = "logFC (high- vs low-WEBS)", y = "-log10(p)",
       title = "Phosphosites (triangles = spliceosome sig)")

# kinases bar (grey shades by significance)
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
kk <- intersect(colnames(pm), webs$sample)
kinases <- intersect(c("SRPK1","SRPK2","CLK1","CLK2","CLK3","CLK4","CDK7","CDK9","CDK12","CDK13"), rownames(pm))
kdf <- data.frame(kinase = kinases, rho = NA, p = NA)
for (k in kinases) {
  ct <- suppressWarnings(cor.test(pm[k, kk], webs$WEBS[match(kk, webs$sample)], method = "spearman"))
  kdf$rho[kdf$kinase == k] <- ct$estimate; kdf$p[kdf$kinase == k] <- ct$p.value
}
kdf$sig <- ifelse(kdf$p < 0.05, "p < 0.05", "n.s.")
p4c <- ggplot(kdf, aes(reorder(kinase, rho), rho, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("p < 0.05" = "grey20", "n.s." = "grey75"), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Spearman rho vs WEBS", title = "Splicing-kinase abundance vs WEBS")

# GENIE3 (grey: writers dark, erasers mid)
ko <- read.csv("results/virtual_ko_abs_summary.csv")
ko$type <- ifelse(ko$ko %in% c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP"), "Writers", "Erasers")
p4d <- ggplot(ko, aes(reorder(ko, ratio), ratio, fill = type)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("Writers" = "grey15", "Erasers" = "grey60"), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Relative |shift| (spliceosome / other)",
       title = "GENIE3 connectivity sensitivity")

fig4_bw <- ((p4a | p4b) / (p4c | p4d)) + plot_annotation(tag_levels = list(c("a","","b","c","d")))
ggsave("figures_BW/Figure_4.pdf", fig4_bw, width = 12, height = 10)
ggsave("figures_BW/Figure_4.png", fig4_bw, width = 12, height = 10, dpi = 200)
cat("BW F2, F4 done.\n")
