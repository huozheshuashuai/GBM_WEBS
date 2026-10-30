# =============================================================================
# 20. Site-level acetylation of spliceosome proteins (mechanism deep-dive)
# =============================================================================
# Hypothesis: writer enzymes acetylate splicing factors; lysine acetylation
# competitively blocks ubiquitination -> protein stabilization -> spliceosome
# accumulation despite lower mRNA (explains the mRNA-protein discordance).
# Tests: (a) spliceosome acetylation sites up in LHS vs LLS (limma)
#        (b) per-site acetylation ~ WEBS correlation
#        (c) acetylation-protein abundance coupling per spliceosome protein
# Input: data/raw/CPTAC/ acetylome + proteome + results/cptac_webs.csv
# Output: results/acetyl_*.csv, figures/fig_acetyl_*.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})
source("08_cptac_proteomics.R")

# ---- Load ----
acetyl <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct",
                               row.names = 1, check.names = FALSE))
proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))
webs <- read.csv("results/cptac_webs.csv", stringsAsFactors = FALSE)

# Parse site IDs: expect formats like "GENE-K123" or "GENE_K123" — inspect
site_ids <- rownames(acetyl)
cat("Example site IDs:", paste(head(site_ids, 3), collapse = ", "), "\n")
site_gene <- sub("[-_\\.][Kk][0-9]+.*$", "", site_ids)  # gene symbol part before K-site

splice_genes <- get_spliceosome_genes()
is_splice_site <- site_gene %in% splice_genes
cat("Acetylation sites on spliceosome proteins:", sum(is_splice_site), "\n")

# ---- (a) Limma LHS vs LLS on all acetylation sites ----
grp <- webs$WEBS_group[match(colnames(acetyl), webs$sample)]
keep <- grp %in% c("LHS", "LLS")
A <- acetyl[, keep]
g <- factor(grp[keep], levels = c("LLS", "LHS"))
# drop sites with too many NAs
ok <- rowSums(!is.na(A)) >= 10
A <- A[ok, ]
design <- model.matrix(~ 0 + g)
colnames(design) <- levels(g)
fit <- eBayes(contrasts.fit(lmFit(A, design), makeContrasts(LHSvsLLS = LHS - LLS, levels = design)))
de <- topTable(fit, number = Inf, sort.by = "p", adjust.method = "BH")
de$site <- rownames(de)
de$gene <- site_gene[ok][match(de$site, rownames(A))]
de$is_spliceosome <- de$gene %in% splice_genes
write.csv(de, "results/acetyl_sites_limma.csv", row.names = FALSE)

# spliceosome site enrichment among significant up-acetylated sites
sig_up <- !is.na(de$adj.P.Val) & de$adj.P.Val < 0.05 & de$logFC > 0
tab <- matrix(c(
  sum(sig_up & de$is_spliceosome), sum(sig_up & !de$is_spliceosome),
  sum(!sig_up & de$is_spliceosome), sum(!sig_up & !de$is_spliceosome)), nrow = 2)
ft <- fisher.test(tab)

# ---- (b) Per-site correlation with WEBS ----
w <- webs$WEBS[match(colnames(A), webs$sample)]
cor_res <- t(sapply(seq_len(nrow(A)), function(i) {
  x <- A[i, ]
  ok2 <- !is.na(x) & !is.na(w)
  if (sum(ok2) < 15) return(c(r = NA, p = NA))
  ct <- suppressWarnings(cor.test(w[ok2], x[ok2], method = "spearman"))
  c(r = unname(ct$estimate), p = ct$p.value)
}))
site_corr <- data.frame(site = rownames(A), gene = de$gene[match(rownames(A), de$site)],
                        spearman_r = cor_res[, 1], p = cor_res[, 2])
site_corr$q <- p.adjust(site_corr$p, method = "BH")
site_corr$is_spliceosome <- site_corr$gene %in% splice_genes
write.csv(site_corr, "results/acetyl_site_webs_corr.csv", row.names = FALSE)

splice_corr <- site_corr[site_corr$is_spliceosome & !is.na(site_corr$spearman_r), ]
n_sig_pos <- sum(splice_corr$q < 0.05 & splice_corr$spearman_r > 0, na.rm = TRUE)

# ---- (c) Acetylation-abundance coupling for top spliceosome proteins ----
# mean acetylation per protein vs protein abundance (same samples)
mean_acetyl_by_gene <- t(sapply(split(seq_len(nrow(A)), site_gene[ok]), function(idx)
  colMeans(A[idx, , drop = FALSE], na.rm = TRUE)))
common_g <- intersect(rownames(mean_acetyl_by_gene), rownames(proteome))
common_g <- intersect(common_g, splice_genes)
coupling <- sapply(common_g, function(gg) {
  s <- intersect(colnames(proteome), colnames(mean_acetyl_by_gene))
  cor(proteome[gg, s], mean_acetyl_by_gene[gg, s],
      method = "spearman", use = "pairwise.complete.obs")
})
coupling_df <- data.frame(gene = common_g, acetyl_protein_coupling = coupling)
write.csv(coupling_df, "results/acetyl_protein_coupling.csv", row.names = FALSE)

# ---- Summary ----
sink("results/acetyl_summary.txt")
cat("=== Site-level acetylation: spliceosome mechanism test ===\n\n")
cat("Total sites tested:", nrow(A), "| spliceosome sites:", sum(de$is_spliceosome), "\n")
cat("Significantly up-acetylated in LHS (FDR<0.05):", sum(sig_up), "\n")
cat("  of which spliceosome:", sum(sig_up & de$is_spliceosome), "\n")
cat("Fisher enrichment OR:", round(ft$estimate, 2), " p:", signif(ft$p.value, 3), "\n\n")
cat("Spliceosome sites positively correlating with WEBS (FDR<0.05):",
    n_sig_pos, "/", nrow(splice_corr), "\n")
cat("Median site-WEBS rho (spliceosome):", round(median(splice_corr$spearman_r, na.rm = TRUE), 3), "\n\n")
cat("Acetylation-abundance coupling (spliceosome proteins, median rho):",
    round(median(coupling, na.rm = TRUE), 3), "\n")
cat("Top coupled proteins:\n")
print(head(coupling_df[order(-coupling_df$acetyl_protein_coupling), ], 10))
sink()
cat(readLines("results/acetyl_summary.txt"), sep = "\n")

# ---- Figure: heatmap of top spliceosome sites ----
dir.create("figures", showWarnings = FALSE)
top_sites <- de[de$is_spliceosome, ] %>% filter(!is.na(adj.P.Val)) %>%
  arrange(adj.P.Val) %>% head(40) %>% pull(site)
if (length(top_sites) >= 5) {
  mat <- A[top_sites, ]
  mat <- mat[complete.cases(mat), ]  # hclust cannot handle NA
  ann <- data.frame(WEBS_group = g, row.names = colnames(A))
  pdf("figures/fig_acetyl_spliceosome_heatmap.pdf", width = 9, height = 8)
  pheatmap(mat, scale = "row", annotation_col = ann,
           show_colnames = FALSE, cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Top 40 spliceosome acetylation sites (LHS vs LLS)")
  dev.off()
}
message("Done.")
