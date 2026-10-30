# =============================================================================
# CPTAC Proteomics Analysis
# =============================================================================
# Data: Gene-level proteome (11,141 proteins x 99 samples)
#       + Site-level acetylome (13,556 sites x 99 samples)
# Methods: Limma for differential expression, GSEA, KEGG
# =============================================================================

library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)

# ---- 1. Load CPTAC Data ----
load_cptac <- function(proteome_file, acetylome_file = NULL, clinical_file = NULL) {
  
  # Proteome: rows = proteins, cols = samples
  proteome <- readRDS(proteome_file)  # or read.table/read.csv
  
  # Acetylome (optional)
  if (!is.null(acetylome_file)) {
    acetylome <- readRDS(acetylome_file)
  } else {
    acetylome <- NULL
  }
  
  # Clinical metadata with WEBS group assignments
  if (!is.null(clinical_file)) {
    clinical <- readRDS(clinical_file)
  } else {
    clinical <- NULL
  }
  
  return(list(proteome = proteome, acetylome = acetylome, clinical = clinical))
}

# ---- 2. WEBS Stratification at Protein Level ----
stratify_cptac <- function(proteome, clinical = NULL) {
  source("tools/WEBS/webs.R")
  
  webs_df <- webs_compute_protein(
    proteome,
    writer_proteins = c("HAT1", "KAT2A", "KAT2B", "KAT5", "KAT6A",
                        "KAT7", "KAT8", "CREBBP", "EP300", "KAT14"),
    eraser_proteins = c("HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5",
                        "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10",
                        "HDAC11", "SIRT1", "SIRT2", "SIRT3", "SIRT5", "SIRT6")
  )
  
  # IQR-based stratification
  q1 <- quantile(webs_df$WEBS, 0.25)
  q3 <- quantile(webs_df$WEBS, 0.75)
  
  webs_df$WEBS_group <- ifelse(
    webs_df$WEBS <= q1, "LLS",
    ifelse(webs_df$WEBS >= q3, "LHS", "LDTS")
  )
  
  return(webs_df)
}

# ---- 3. Differential Protein Expression (LHS vs LLS) ----
run_limma_proteomics <- function(proteome, group_vector, 
                                  lhs_label = "LHS", lls_label = "LLS") {
  
  # Subset to LHS and LLS only
  keep <- group_vector %in% c(lhs_label, lls_label)
  proteome_sub <- proteome[, keep]
  group_sub <- factor(group_vector[keep], levels = c(lls_label, lhs_label))
  
  # Design matrix
  design <- model.matrix(~ 0 + group_sub)
  colnames(design) <- levels(group_sub)
  
  # Limma pipeline
  fit <- lmFit(proteome_sub, design)
  
  contrast <- makeContrasts(LHSvsLLS = LHS - LLS, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Results
  results <- topTable(fit2, number = Inf, sort.by = "p", adjust.method = "BH")
  results$protein <- rownames(results)
  
  # Add direction
  results$direction <- ifelse(results$logFC > 0, "Up", "Down")
  
  return(list(results = results, fit = fit2))
}

# ---- 4. Spliceosome Protein Analysis ----
analyze_spliceosome_proteins <- function(limma_results, spliceosome_genes) {
  
  # Drop NA fit rows, then filter to spliceosome proteins
  res_all <- limma_results$results[!is.na(limma_results$results$logFC) &
                                   !is.na(limma_results$results$adj.P.Val), ]
  splice_results <- res_all[rownames(res_all) %in% spliceosome_genes, ]

  n_up <- sum(splice_results$logFC > 0 & splice_results$adj.P.Val < 0.05)
  n_down <- sum(splice_results$logFC < 0 & splice_results$adj.P.Val < 0.05)
  n_total <- nrow(splice_results)

  # Fisher's exact test for enrichment
  all_up <- sum(res_all$logFC > 0 & res_all$adj.P.Val < 0.05)
  all_down <- sum(res_all$logFC < 0 & res_all$adj.P.Val < 0.05)
  all_total <- nrow(res_all)
  
  contingency <- matrix(
    c(n_up, all_up - n_up,
      n_down, all_down - n_down),
    nrow = 2
  )
  
  fisher <- fisher.test(contingency, alternative = "two.sided")
  
  return(list(
    splice_results = splice_results,
    summary = data.frame(
      n_spliceosome = n_total,
      n_up = n_up,
      n_down = n_down,
      odds_ratio = fisher$estimate,
      p_value = fisher$p.value
    )
  ))
}

# ---- 5. Transcript-Protein Correlation (Cross-Layer) ----
compute_cross_layer_corr <- function(rna_matrix, protein_matrix, shared_samples) {
  
  # Match genes
  shared_genes <- intersect(rownames(rna_matrix), rownames(protein_matrix))
  
  rna_sub <- rna_matrix[shared_genes, shared_samples]
  protein_sub <- protein_matrix[shared_genes, shared_samples]
  
  # Spearman correlation per gene
  correlations <- sapply(shared_genes, function(g) {
    cor(rna_sub[g, ], protein_sub[g, ], method = "spearman", use = "pairwise.complete.obs")
  })
  
  # Focus on spliceosome genes
  spliceosome_genes <- intersect(shared_genes, get_spliceosome_genes())
  splice_corr <- correlations[spliceosome_genes]
  
  return(list(
    all_genes = correlations,
    spliceosome = splice_corr,
    mean_corr_all = mean(correlations, na.rm = TRUE),
    mean_corr_splice = mean(splice_corr, na.rm = TRUE)
  ))
}

# ---- 6. Helper: Spliceosome Gene List ----
get_spliceosome_genes <- function() {
  # Major spliceosome components
  c("SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5", "SF3B6",
    "U2AF1", "U2AF2", "U2AF1L4", "U2AF1L5",
    "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF31", "PRPF38A", "PRPF38B", "PRPF39", "PRPF40A", "PRPF40B",
    "SNRNP70", "SNRPA", "SNRPB", "SNRPC", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG",
    "SNRPA1", "SNRPB2", "SNRNP200", "SNRNP25", "SNRNP27", "SNRNP35", "SNRNP40", "SNRNP48",
    "EFTUD2", "SNU13", "LSM2", "LSM3", "LSM4", "LSM5", "LSM6", "LSM7", "LSM8",
    "SMNDC1", "GEMIN2", "GEMIN4", "GEMIN5", "GEMIN6", "GEMIN7", "GEMIN8",
    "DDX39B", "DDX42", "DDX46",
    "TCERG1", "RBM22", "RBM25", "RBM39", "RBM45", "SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7", "SRSF8", "SRSF9", "SRSF10", "SRSF11",
    "TRA2A", "TRA2B", "HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPK", "HNRNPL", "HNRNPM",
    "PHF5A", "PPIE", "PPIH", "CWC15", "CWC22", "CWC25", "CWC27",
    "AQR", "ISY1", "CDC40", "CDC5L", "PLRG1", "BCAS2",
    "DHX8", "DHX15", "DHX16", "DHX35", "DHX38",
    "SYF1", "SYF2", "CRNKL1", "BUD31", "SART1",
    "SNW1", "SKIV2L2", "THOC1", "THOC2", "THOC5", "THOC6", "THOC7",
    "ALYREF", "NXF1", "NXT1", "UAP56", "URH49")
}
