# =============================================================================
# Driver: CPTAC-GBM proteomics - the central spliceosome analysis
# =============================================================================
# Reproduces: protein-level WEBS, limma LHS vs LLS, spliceosome enrichment
# (paper: 117/129 up, OR = 39.2, FDR = 6.9e-64), cross-layer mRNA-protein corr
# Input: data/raw/CPTAC/HS_CPTAC_GBM_*.cct
# Output: results/cptac_*.csv, results/cptac_limma_results.csv
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
})
source("tools/WEBS/webs.R")
source("08_cptac_proteomics.R")   # for get_spliceosome_genes()

# ---- 1. Load CPTAC data (gene x sample format) ----
proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))
acetylome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct",
                                  row.names = 1, check.names = FALSE))
rnaseq <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct",
                               row.names = 1, check.names = FALSE))
# drop all-NA rows/cols
proteome <- proteome[rowSums(!is.na(proteome)) > 0, colSums(!is.na(proteome)) > 0]
rnaseq <- rnaseq[rowSums(!is.na(rnaseq)) > 0, ]
message("Proteome: ", nrow(proteome), " x ", ncol(proteome))
message("Acetylome: ", nrow(acetylome), " x ", ncol(acetylome))
message("RNAseq: ", nrow(rnaseq), " x ", ncol(rnaseq))

# ---- 2. Protein-level WEBS and stratification ----
webs_df <- webs_compute_protein(proteome)
q <- quantile(webs_df$WEBS, c(0.25, 0.75), na.rm = TRUE)
webs_df$WEBS_group <- ifelse(webs_df$WEBS <= q[1], "LLS",
                      ifelse(webs_df$WEBS >= q[2], "LHS", "LDTS"))
write.csv(webs_df, "results/cptac_webs.csv", row.names = FALSE)
print(table(webs_df$WEBS_group))

# ---- 3. Limma differential abundance: LHS vs LLS ----
limma_res <- run_limma_proteomics(proteome, webs_df$WEBS_group)
de <- limma_res$results
de$is_spliceosome <- rownames(de) %in% get_spliceosome_genes()
write.csv(de, "results/cptac_limma_results.csv")

# ---- 4. Spliceosome enrichment (paper: 117/129 up, OR = 39.2) ----
sp <- analyze_spliceosome_proteins(limma_res, get_spliceosome_genes())
write.csv(sp$summary, "results/cptac_spliceosome_summary.csv", row.names = FALSE)
write.csv(sp$splice_results, "results/cptac_spliceosome_proteins.csv")

# ---- 5. Histone acetylation correlation ----
# Overall histone acetylation = mean across histone H1/H2A/H2B/H3/H4 sites
acetyl_genes <- sub("\\..*$", "", rownames(acetylome))  # site ids like GENE-K123
histone_pat <- "^H1|^H2A|^H2B|^H3|^H4|HIST"
hist_rows <- grepl(histone_pat, acetyl_genes)
if (sum(hist_rows) > 10) {
  hist_ac <- colMeans(acetylome[hist_rows, , drop = FALSE], na.rm = TRUE)
  common <- intersect(webs_df$sample, names(hist_ac))
  ct <- cor.test(webs_df$WEBS[match(common, webs_df$sample)], hist_ac[common],
                 method = "pearson")
  message("Protein WEBS vs overall histone acetylation: r = ",
          round(ct$estimate, 3), ", p = ", signif(ct$p.value, 3))
  acetyl_res <- data.frame(metric = "WEBS_vs_histone_acetylation",
                           r = ct$estimate, p = ct$p.value)
  write.csv(acetyl_res, "results/cptac_acetylation_corr.csv", row.names = FALSE)
}

# ---- 6. Cross-layer correlation (mRNA vs protein for spliceosome) ----
common_samples <- intersect(colnames(rnaseq), colnames(proteome))
shared_genes <- intersect(rownames(rnaseq), rownames(proteome))
splice_genes <- intersect(get_spliceosome_genes(), shared_genes)

rna_sub <- rnaseq[splice_genes, common_samples]
prot_sub <- proteome[splice_genes, common_samples]
splice_corr <- sapply(splice_genes, function(g) {
  cor(rna_sub[g, ], prot_sub[g, ], method = "spearman", use = "pairwise.complete.obs")
})

# also: correlation of spliceosome score with WEBS at mRNA vs protein level
splice_score_prot <- colMeans(prot_sub, na.rm = TRUE)
splice_score_rna <- colMeans(rna_sub, na.rm = TRUE)
webs_common <- webs_df$WEBS[match(common_samples, webs_df$sample)]
# mRNA-level WEBS must be computed from the RNA-seq matrix (transcriptomic gene sets)
rna_webs <- webs_compute_bulk(rnaseq)
webs_rna_common <- rna_webs$WEBS[match(common_samples, rna_webs$sample)]
ct_prot <- cor.test(webs_common, splice_score_prot, method = "spearman")
ct_rna <- cor.test(webs_rna_common, splice_score_rna, method = "spearman")

sink("results/cptac_summary.txt")
cat("=== CPTAC proteomics summary ===\n\n")
cat("Samples:", ncol(proteome), "| Proteins:", nrow(proteome), "\n\n")
cat("-- Spliceosome enrichment (LHS vs LLS) --\n")
cat("Quantified spliceosome proteins:", sp$summary$n_spliceosome, "\n")
cat("Up-regulated:", sp$summary$n_up, " Down:", sp$summary$n_down, "\n")
cat("Odds ratio:", round(sp$summary$odds_ratio, 1), " p:", sp$summary$p_value, "\n")
cat("(Paper: 117/129 up, OR = 39.2, FDR = 6.9e-64)\n\n")
cat("-- Spliceosome machinery vs WEBS across samples --\n")
cat("Protein level: Spearman rho =", round(ct_prot$estimate, 3),
    " p =", signif(ct_prot$p.value, 3), "\n")
cat("mRNA level:    Spearman rho =", round(ct_rna$estimate, 3),
    " p =", signif(ct_rna$p.value, 3), "\n")
cat("(Paper: protein rho = 0.44, mRNA rho = -0.29)\n\n")
cat("-- Per-gene mRNA-protein correlation (spliceosome) --\n")
cat("Median rho:", round(median(splice_corr, na.rm = TRUE), 3), "\n")
sink()

cat(readLines("results/cptac_summary.txt"), sep = "\n")
message("Done.")
