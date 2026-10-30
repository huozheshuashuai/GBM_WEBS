# =============================================================================
# Pharmacogenomics: DepMap CRISPR + GDSC2 Drug Screening
# =============================================================================
# DepMap: Avana CRISPR dependency (Chronos scores, negative = lethal)
# GDSC2: Pan-drug screening (IC50, AUC)
# Correlation: WEBS vs gene dependency, WEBS vs drug sensitivity
# =============================================================================

library(dplyr)
library(ggplot2)

# ---- 1. Load DepMap Data ----
load_depmap <- function(effects_file, model_file, expression_file) {
  
  crispr <- read.csv(effects_file, row.names = 1, check.names = FALSE)
  models <- read.csv(model_file)
  expression <- read.csv(expression_file, row.names = 1, check.names = FALSE)
  
  # Filter to CNS/GBM cell lines
  gbm_lines <- models$ModelID[grepl("CNS|BRAIN|GLIOMA|ASTROCYTOMA", 
                                     models$OncotreePrimaryDisease, ignore.case = TRUE)]
  
  crispr_gbm <- crispr[, intersect(colnames(crispr), gbm_lines)]
  expr_gbm <- expression[, intersect(colnames(expression), gbm_lines)]
  
  return(list(crispr = crispr_gbm, expression = expr_gbm, models = models))
}

# ---- 2. Compute WEBS for Cell Lines ----
compute_depmap_webs <- function(expression_matrix) {
  source("tools/WEBS/webs.R")
  
  webs_df <- webs_compute_bulk(
    expression_matrix,
    writer_genes = c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP"),
    eraser_genes = c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")
  )
  
  return(webs_df)
}

# ---- 3. Correlation: WEBS vs Gene Dependency ----
correlate_webs_dependency <- function(webs_df, crispr_matrix) {
  
  common_lines <- intersect(webs_df$sample, colnames(crispr_matrix))
  webs <- webs_df$WEBS[match(common_lines, webs_df$sample)]
  
  correlations <- apply(crispr_matrix[, common_lines], 1, function(dep) {
    cor(webs, dep, method = "spearman", use = "pairwise.complete.obs")
  })
  
  pvalues <- apply(crispr_matrix[, common_lines], 1, function(dep) {
    cor.test(webs, dep, method = "spearman")$p.value
  })
  
  results <- data.frame(
    gene = rownames(crispr_matrix),
    spearman_r = correlations,
    p_value = pvalues,
    q_value = p.adjust(pvalues, method = "BH")
  ) %>%
    arrange(p_value)
  
  return(results)
}

# ---- 4. Find Synthetic Lethal Candidates ----
find_synthetic_lethal <- function(cor_results, top_n = 20) {
  
  # Positive correlation: high WEBS = less dependent (less negative Chronos)
  # Negative correlation: high WEBS = more dependent (more negative Chronos) = synthetic lethal
  
  sl_candidates <- cor_results %>%
    filter(spearman_r < -0.3, q_value < 0.05) %>%
    arrange(spearman_r) %>%
    head(top_n)
  
  return(sl_candidates)
}

# ---- 5. Load GDSC2 Data ----
load_gdsc2 <- function(response_file, expression_file) {
  
  response <- read.delim(response_file)
  expression <- read.delim(expression_file)
  
  return(list(response = response, expression = expression))
}

# ---- 6. Correlation: WEBS vs Drug Sensitivity ----
correlate_webs_drug <- function(webs_df, gdsc_response, drug_name = NULL) {
  
  # Filter to specific drug or all drugs
  if (!is.null(drug_name)) {
    gdsc_response <- gdsc_response[gdsc_response$DRUG_NAME == drug_name, ]
  }
  
  common_lines <- intersect(webs_df$sample, gdsc_response$CELL_LINE_NAME)
  webs <- webs_df$WEBS[match(common_lines, webs_df$sample)]
  
  # Correlation per drug
  drug_cor <- gdsc_response %>%
    group_by(DRUG_NAME) %>%
    filter(n() >= 5) %>%  # At least 5 cell lines
    summarise(
      n_lines = n(),
      spearman_r = cor(LN_IC50, webs[match(CELL_LINE_NAME, common_lines)], 
                       method = "spearman", use = "pairwise.complete.obs"),
      p_value = cor.test(LN_IC50, webs[match(CELL_LINE_NAME, common_lines)], 
                         method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  return(drug_cor)
}

# ---- 7. HDAC Inhibitor Analysis ----
analyze_hdac_inhibitors <- function(gdsc_response, webs_df) {
  
  hdac_drugs <- gdsc_response %>%
    filter(grepl("HDAC|VORINOSTAT|SAHA|ROMIDEPSIN|PANOBINOSTAT|BELINOSTAT|entinostat",
                 DRUG_NAME, ignore.case = TRUE))
  
  results <- correlate_webs_drug(webs_df, hdac_drugs)
  
  return(results)
}

# ---- 8. Wilcoxon: LHS vs LLS Drug Sensitivity ----
compare_drug_sensitivity <- function(gdsc_response, webs_df, drug_name) {
  
  drug_data <- gdsc_response[gdsc_response$DRUG_NAME == drug_name, ]
  
  # Merge WEBS
  drug_data$WEBS <- webs_df$WEBS[match(drug_data$CELL_LINE_NAME, webs_df$sample)]
  drug_data <- drug_data[!is.na(drug_data$WEBS), ]
  
  # Stratify
  q1 <- quantile(drug_data$WEBS, 0.25)
  q3 <- quantile(drug_data$WEBS, 0.75)
  drug_data$group <- ifelse(drug_data$WEBS <= q1, "LLS",
                            ifelse(drug_data$WEBS >= q3, "LHS", "LDTS"))
  
  # Compare LHS vs LLS
  lhs_ic50 <- drug_data$LN_IC50[drug_data$group == "LHS"]
  lls_ic50 <- drug_data$LN_IC50[drug_data$group == "LLS"]
  
  test <- wilcox.test(lhs_ic50, lls_ic50)
  
  return(list(
    test = test,
    lhs_median = median(lhs_ic50, na.rm = TRUE),
    lls_median = median(lls_ic50, na.rm = TRUE),
    data = drug_data
  ))
}
