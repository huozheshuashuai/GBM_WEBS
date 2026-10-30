# =============================================================================
# Figure Generation
# =============================================================================
# Standardized color schemes and plotting functions for all main figures
# =============================================================================

library(ggplot2)
library(patchwork)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# ---- Color Schemes ----
COLORS <- list(
  webs_group = c("LLS" = "#1f77b4", "LDTS" = "#ff7f0e", "LHS" = "#d62728"),
  neftel = c("MES1" = "#8c564b", "MES2" = "#c49c94", 
             "AC" = "#2ca02c", "NPC1" = "#9467bd", "NPC2" = "#bcbd22", "OPC" = "#17becf"),
  cell_type = c("Malignant" = "#d62728", "Macrophage_Microglia" = "#ff7f0e",
                "Oligodendrocyte" = "#2ca02c", "Astrocyte" = "#9467bd",
                "Neuron" = "#8c564b", "Endothelial" = "#e377c2", "Pericyte" = "#7f7f7f"),
  region = c("Core" = "#d62728", "Invasion" = "#ff7f0e", 
             "Normal" = "#2ca02c", "Intermediate" = "#9467bd")
)

# ---- Figure 1: Study Design ----
plot_fig1_design <- function() {
  # Conceptual diagram - typically created in Illustrator/PowerPoint
  message("Figure 1 is a conceptual diagram. Create in Illustrator with:")
  message("  - Data sources (scRNA, snATAC, Spatial, Bulk, CPTAC, DepMap)")
  message("  - WEBS computation flow")
  message("  - Three-tier analysis (Discovery, External, Pharmacogenomic)")
}

# ---- Figure 2: scRNA-seq WEBS Landscape ----
plot_fig2_landscape <- function(obj) {
  p1 <- DimPlot(obj, reduction = "umap", group.by = "cell_type", 
                cols = COLORS$cell_type, label = TRUE) + 
    ggtitle("Cell Types")
  
  p2 <- FeaturePlot(obj, features = "WEBS", cols = c("blue", "white", "red")) +
    ggtitle("WEBS Score")
  
  p3 <- VlnPlot(obj, features = "WEBS", group.by = "cell_type", 
                cols = COLORS$cell_type) + 
    ggtitle("WEBS by Cell Type") + NoLegend()
  
  p4 <- DimPlot(obj, reduction = "umap", group.by = "WEBS_group",
                cols = COLORS$webs_group) +
    ggtitle("WEBS Groups (IQR)")
  
  (p1 + p2) / (p3 + p4)
}

# ---- Figure 3: Neftel Module vs WEBS ----
plot_fig3_neftel <- function(obj) {
  neftel_features <- c("Neftel_MES1", "Neftel_MES2", "Neftel_AC", 
                       "Neftel_NPC1", "Neftel_NPC2", "Neftel_OPC")
  
  VlnPlot(obj, features = neftel_features, group.by = "WEBS_group",
          cols = COLORS$webs_group, pt.size = 0) +
    plot_layout(ncol = 3) +
    ggtitle("Neftel Modules by WEBS Group")
}

# ---- Figure 4: AUCell Pathway Heatmap ----
plot_fig4_pathway_heatmap <- function(pathway_results, top_n = 30) {
  
  # Select top pathways by effect size
  top <- pathway_results %>%
    arrange(desc(abs(effect_size_r))) %>%
    head(top_n)
  
  # Prepare matrix for heatmap
  mat <- matrix(c(top$lls_mean, top$lhs_mean), ncol = 2)
  colnames(mat) <- c("LLS", "LHS")
  rownames(mat) <- top$pathway
  
  pheatmap(mat, scale = "row", color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Top Pathways: LHS vs LLS",
           cluster_cols = FALSE)
}

# ---- Figure 5: Survival K-M Curves ----
plot_fig5_survival <- function(clinical_data, cohort_name) {
  plot_km(clinical_data, cohort_name)
}

# ---- Figure 6: Six-Gene Signature ----
plot_fig6_signature <- function(cox_results, cindex_results) {
  # Forest plot of HR across cohorts
  library(forestplot)
  
  forest_data <- data.frame(
    cohort = cox_results$cohort,
    HR = cox_results$HR,
    lower = cox_results$lower_95,
    upper = cox_results$upper_95,
    p = cox_results$p_value
  )
  
  forestplot(
    labeltext = c(forest_data$cohort, format(forest_data$p, digits = 2)),
    mean = forest_data$HR,
    lower = forest_data$lower,
    upper = forest_data$upper,
    zero = 1,
    boxsize = 0.3,
    lineheight = unit(8, "mm"),
    col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue")
  )
}

# ---- Figure 7: CPTAC Volcano + Spliceosome ----
plot_fig7_cptac <- function(limma_results, spliceosome_results) {
  
  # Volcano plot
  p1 <- ggplot(limma_results$results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = direction), alpha = 0.5) +
    geom_point(data = spliceosome_results$splice_results,
               aes(x = logFC, y = -log10(adj.P.Val)), color = "red", size = 2) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    ggtitle("CPTAC: LHS vs LLS Differential Proteins")
  
  # Spliceosome bar plot
  splice_df <- data.frame(
    Protein = rownames(spliceosome_results$splice_results),
    logFC = spliceosome_results$splice_results$logFC
  ) %>% arrange(logFC)
  
  p2 <- ggplot(splice_df, aes(x = reorder(Protein, logFC), y = logFC)) +
    geom_bar(stat = "identity", aes(fill = logFC > 0)) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    theme_bw() +
    ggtitle("Spliceosome Protein Changes")
  
  p1 + p2
}

# ---- Figure 8: DepMap CRISPR ----
plot_fig8_depmap <- function(depmap_results, top_genes = 20) {
  
  top <- depmap_results %>%
    arrange(q_value) %>%
    head(top_genes)
  
  ggplot(top, aes(x = reorder(gene, spearman_r), y = spearman_r)) +
    geom_bar(stat = "identity", aes(fill = spearman_r < 0)) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                      labels = c("Negative" = "Synthetic Lethal", "Positive" = "Not SL")) +
    labs(x = "Gene", y = "Spearman r (WEBS vs Dependency)",
         title = "Top DepMap Synthetic Lethal Candidates") +
    theme_bw()
}

# ---- Figure 9: GDSC2 Drug Screen ----
plot_fig9_gdsc2 <- function(gdsc_results) {
  
  ggplot(gdsc_results %>% filter(q_value < 0.1), 
         aes(x = reorder(DRUG_NAME, spearman_r), y = spearman_r)) +
    geom_bar(stat = "identity", aes(fill = spearman_r < 0)) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    labs(x = "Drug", y = "Spearman r (WEBS vs LN_IC50)",
         title = "GDSC2: WEBS-Associated Drug Sensitivity") +
    theme_bw()
}

# ---- Figure 10: Calibration + Nomogram ----
plot_fig10_calibration <- function(clinical_data, time_point = 24) {
  cal <- plot_calibration(clinical_data, time_point)
  return(cal)
}

# ---- Figure 11: Benchmarking ----
plot_fig11_benchmarking <- function(benchmark_results) {
  
  library(reshape2)
  df <- melt(benchmark_results, id.vars = "method")
  colnames(df) <- c("Method", "Time", "AUC")
  
  ggplot(df, aes(x = Method, y = AUC, fill = Method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Time) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    coord_cartesian(ylim = c(0.3, 0.9)) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Method Comparison: Survival Prediction AUC")
}

# ---- Figure 12: Working Model ----
plot_fig12_model <- function() {
  message("Figure 12 is a conceptual working model.")
  message("Create in Illustrator with:")
  message("  - Left: Epigenetic Writer-Eraser Balance axis (histone lactylation)")
  message("  - Center: Spliceosome activation in LHS cells")
  message("  - Right: Therapeutic implications (HDAC inhibitors, splicing modulators)")
  message("  - Key genes: SIRT2, HDAC1, MMP2, EP300, KAT2A, HDAC2")
}
