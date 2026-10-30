# =============================================================================
# Pathway Analysis (AUCell) & TF Regulon Analysis (SCENIC / pySCENIC)
# =============================================================================
# AUCell: Per-cell pathway enrichment via ranking-based AUC
# SCENIC: RcisTarget + GRNBoost2 for TF regulon inference
# =============================================================================

library(AUCell)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RCisTarget)
library(GENIE3)

# ---- 1. AUCell Pathway Scoring ----
# Paper: 10,507 pathways tested, leading-edge = 0.375

score_pathways_aucell <- function(obj, pathway_list, aucMaxRank = 0.375) {
  
  # Build ranking for each cell
  expr_matrix <- GetAssayData(obj, layer = "data")
  cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
  
  # Score pathways
  cells_AUC <- AUCell_calcAUC(
    geneSets = pathway_list,
    rankings = cells_rankings,
    aucMaxRank = ceiling(aucMaxRank * nrow(expr_matrix)),
    normAUC = TRUE
  )
  
  # Extract AUC matrix
  auc_matrix <- getAUC(cells_AUC)
  
  # Add to metadata
  for (pw in rownames(auc_matrix)) {
    obj[[pw]] <- auc_matrix[pw, colnames(obj)]
  }
  
  return(list(obj = obj, auc = cells_AUC))
}

# ---- 2. Pathway Comparison: LHS vs LLS ----
compare_pathways <- function(obj, pathways, group_col = "WEBS_group") {
  
  lhs_cells <- WhichCells(obj, expression = !!sym(group_col) == "LHS")
  lls_cells <- WhichCells(obj, expression = !!sym(group_col) == "LLS")
  
  results <- lapply(pathways, function(pw) {
    lhs_scores <- obj[[pw]][lhs_cells, 1]
    lls_scores <- obj[[pw]][lls_cells, 1]
    
    # Wilcoxon rank-sum test (Mann-Whitney U)
    test <- wilcox.test(lhs_scores, lls_scores)
    
    # Effect size: rank-biserial correlation
    # r = Z / sqrt(N)
    z <- qnorm(test$p.value / 2) * sign(mean(lhs_scores) - mean(lls_scores))
    n <- length(lhs_cells) + length(lls_cells)
    r <- z / sqrt(n)
    
    data.frame(
      pathway = pw,
      lhs_mean = mean(lhs_scores),
      lls_mean = mean(lls_scores),
      p_value = test$p.value,
      effect_size_r = r
    )
  })
  
  results <- do.call(rbind, results)
  results$q_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)
}

# ---- 3. GSEA on Bulk RNA-seq ----
run_gsea_bulk <- function(ranked_genes, gmt_file, nPerm = 1000) {
  library(fgsea)
  
  # ranked_genes: named numeric vector (log2FC), names = gene symbols
  pathways <- gmtPathways(gmt_file)
  
  fgsea_results <- fgsea(
    pathways = pathways,
    stats = ranked_genes,
    nperm = nPerm,
    minSize = 15,
    maxSize = 500
  )
  
  fgsea_results <- fgsea_results[order(fgsea_results$padj), ]
  return(fgsea_results)
}

# ---- 4. KEGG Enrichment (CPTAC Proteins) ----
run_kegg_proteins <- function(protein_list, organism = "hsa", pvalueCutoff = 0.05) {
  
  # Convert protein symbols to Entrez IDs
  entrez_ids <- bitr(protein_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  kegg <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = organism,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  return(kegg)
}

# ---- 5. SCENIC / pySCENIC Pipeline (R wrapper) ----
# Note: pySCENIC is Python-based. Here is the R equivalent workflow.

run_scenic_r <- function(obj, tf_list, nCores = 4) {
  
  expr_matrix <- as.matrix(GetAssayData(obj, layer = "data"))
  
  # Step 1: GRNBoost2 equivalent (GENIE3 in R)
  # Use only variable features for speed
  var_genes <- VariableFeatures(obj)
  tf_var <- intersect(tf_list, var_genes)
  
  weightMatrix <- GENIE3(expr_matrix[var_genes, ], 
                          regulators = tf_var,
                          nCores = nCores)
  
  # Step 2: RcisTarget motif enrichment
  # Requires RcisTarget databases (download separately)
  
  # Step 3: AUCell regulon scoring
  # Similar to AUCell pathway scoring above
  
  message("SCENIC R workflow completed. For full pySCENIC, use Python pipeline.")
  return(weightMatrix)
}

# ---- 6. High-Confidence TF Identification ----
# Criteria:
#   a) Significantly associated with WEBS (Kruskal-Wallis)
#   b) AUC > 0.75 in ROC (distinguishing LHS vs LLS)
#   c) Positive correlation with survival HR

identify_high_confidence_tfs <- function(tf_auc_matrix, webs_group, survival_hr) {
  
  # a) Kruskal-Wallis across WEBS groups
  kw_results <- apply(tf_auc_matrix, 1, function(x) {
    kruskal.test(x ~ webs_group)$p.value
  })
  
  # b) AUC for LHS vs LLS
  roc_results <- apply(tf_auc_matrix, 1, function(x) {
    library(pROC)
    roc_obj <- roc(webs_group == "LHS", x)
    auc(roc_obj)
  })
  
  # c) Correlation with survival HR
  hr_corr <- apply(tf_auc_matrix, 1, function(x) {
    cor(x, survival_hr, use = "pairwise.complete.obs", method = "spearman")
  })
  
  results <- data.frame(
    TF = rownames(tf_auc_matrix),
    KW_pvalue = kw_results,
    ROC_AUC = roc_results,
    HR_correlation = hr_corr
  )
  
  # Filter: p < 0.05, AUC > 0.75, positive HR correlation
  high_conf <- results %>%
    filter(KW_pvalue < 0.05, ROC_AUC > 0.75, HR_correlation > 0) %>%
    arrange(desc(ROC_AUC))
  
  return(list(all = results, high_confidence = high_conf))
}
