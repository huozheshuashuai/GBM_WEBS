# =============================================================================
# Method Comparison: WEBS vs ssGSEA / AUCell / PLAGE / All-gene Mean
# =============================================================================
# Tested on: GSE131928 external cohort (28 patients, 24,131 cells)
# Metric: AUC for patient-level survival stratification (1-year, 2-year, 3-year)
# =============================================================================

library(UCell)
library(GSVA)
library(AUCell)
library(dplyr)

# ---- 1. Compute Scores with Different Methods ----

compute_ssgsea <- function(obj, writer_genes, eraser_genes) {
  # ssGSEA: Single-sample GSEA
  expr <- as.matrix(GetAssayData(obj, layer = "data"))
  
  writer_score <- gsva(expr, gset.idx.list = list(Writers = writer_genes),
                       method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
  eraser_score <- gsva(expr, gset.idx.list = list(Erasers = eraser_genes),
                       method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
  
  score <- writer_score[1, ] - eraser_score[1, ]
  return(score)
}

compute_aucell_method <- function(obj, writer_genes, eraser_genes) {
  # AUCell method (different from WEBS's UCell)
  expr <- GetAssayData(obj, layer = "data")
  cells_rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
  
  writer_auc <- AUCell_calcAUC(geneSets = list(Writers = writer_genes),
                                rankings = cells_rankings)
  eraser_auc <- AUCell_calcAUC(geneSets = list(Erasers = eraser_genes),
                                rankings = cells_rankings)
  
  score <- getAUC(writer_auc)["Writers", ] - getAUC(eraser_auc)["Erasers", ]
  return(score)
}

compute_plage <- function(obj, writer_genes, eraser_genes) {
  # PLAGE: Pathway Level Analysis of Gene Expression
  expr <- as.matrix(GetAssayData(obj, layer = "data"))
  
  writer_score <- gsva(expr, gset.idx.list = list(Writers = writer_genes),
                       method = "plage", verbose = FALSE)
  eraser_score <- gsva(expr, gset.idx.list = list(Erasers = eraser_genes),
                       method = "plage", verbose = FALSE)
  
  score <- writer_score[1, ] - eraser_score[1, ]
  return(score)
}

compute_allgene_mean <- function(obj, writer_genes, eraser_genes) {
  # All-gene mean: simple mean expression difference
  expr <- GetAssayData(obj, layer = "data")
  
  writer_mean <- Matrix::colMeans(expr[intersect(writer_genes, rownames(expr)), , drop = FALSE])
  eraser_mean <- Matrix::colMeans(expr[intersect(eraser_genes, rownames(expr)), , drop = FALSE])
  
  score <- writer_mean - eraser_mean
  return(score)
}

# ---- 2. Patient-Level Aggregation ----
aggregate_to_patient <- function(scores, patient_ids) {
  # Median score per patient
  patient_scores <- tapply(scores, patient_ids, median, na.rm = TRUE)
  return(patient_scores)
}

# ---- 3. Compare Methods ----
compare_methods <- function(obj, clinical_data, patient_id_col = "patient") {
  
  writer_genes <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
  eraser_genes <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")
  
  # Compute all scores
  message("Computing ssGSEA...")
  ssgsea <- compute_ssgsea(obj, writer_genes, eraser_genes)
  
  message("Computing AUCell...")
  aucell <- compute_aucell_method(obj, writer_genes, eraser_genes)
  
  message("Computing PLAGE...")
  plage <- compute_plage(obj, writer_genes, eraser_genes)
  
  message("Computing all-gene mean...")
  allgene <- compute_allgene_mean(obj, writer_genes, eraser_genes)
  
  message("Computing WEBS (UCell)...")
  source("tools/WEBS/webs.R")
  obj <- webs_compute(obj, writer_genes, eraser_genes)
  webs <- obj$WEBS
  
  # Aggregate to patient level
  patient_ids <- obj[[patient_id_col]][, 1]
  
  methods <- list(
    WEBS = aggregate_to_patient(webs, patient_ids),
    ssGSEA = aggregate_to_patient(ssgsea, patient_ids),
    AUCell = aggregate_to_patient(aucell, patient_ids),
    PLAGE = aggregate_to_patient(plage, patient_ids),
    AllGeneMean = aggregate_to_patient(allgene, patient_ids)
  )
  
  # Compare AUC for survival prediction
  library(timeROC)
  
  results <- lapply(names(methods), function(m) {
    scores <- methods[[m]]
    
    # Match to clinical
    common_pts <- intersect(names(scores), clinical_data$patient)
    scores_sub <- scores[common_pts]
    clin_sub <- clinical_data[match(common_pts, clinical_data$patient), ]
    
    aucs <- sapply(c(12, 24, 36), function(t) {
      roc <- timeROC(T = clin_sub$OS.time, delta = clin_sub$OS,
                     marker = scores_sub, cause = 1, times = t, iid = FALSE)
      roc$AUC[2]
    })
    
    data.frame(
      method = m,
      AUC_1yr = aucs[1],
      AUC_2yr = aucs[2],
      AUC_3yr = aucs[3]
    )
  })
  
  do.call(rbind, results)
}
