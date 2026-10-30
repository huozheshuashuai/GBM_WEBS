# =============================================================================
# 28. Specificity control: is the mRNA-protein discordance spliceosome-specific?
# =============================================================================
# Compare layer-discordance (protein rho - mRNA rho vs WEBS) for:
#   - spliceosome (111 proteins)
#   - 10 random gene sets of the same size (from quantified proteome)
#   - KEGG proteasome (PSM*) as degradation-pathway control
# Output: results/specificity_control.csv, figures/figS_specificity.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr) })
set.seed(42)
source("08_cptac_proteomics.R")

proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))
rnaseq <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct",
                               row.names = 1, check.names = FALSE))
webs <- read.csv("results/cptac_webs.csv")

common_samples <- intersect(colnames(proteome), intersect(colnames(rnaseq), webs$sample))
prot_w <- webs$WEBS[match(common_samples, webs$sample)]
# mRNA-layer WEBS from CPTAC RNAseq (paper's mRNA-layer configuration)
source("tools/WEBS/webs.R")
rna_webs <- webs_compute_bulk(rnaseq)
mrna_w <- rna_webs$WEBS[match(common_samples, rna_webs$sample)]

layer_cor <- function(genes) {
  g_p <- intersect(genes, rownames(proteome))
  g_r <- intersect(genes, rownames(rnaseq))
  score_p <- colMeans(proteome[g_p, common_samples, drop = FALSE], na.rm = TRUE)
  score_r <- colMeans(rnaseq[g_r, common_samples, drop = FALSE], na.rm = TRUE)
  c(prot = cor(prot_w, score_p, method = "spearman", use = "pairwise.complete.obs"),
    mrna = cor(mrna_w, score_r, method = "spearman", use = "pairwise.complete.obs"),
    n = length(g_p))
}

splice_genes <- get_spliceosome_genes()
splice_res <- layer_cor(splice_genes)

# 100 random gene sets of equal size (stable empirical p)
all_genes <- intersect(rownames(proteome), rownames(rnaseq))
rand_res <- t(sapply(1:100, function(i) {
  layer_cor(sample(all_genes, splice_res["n"]))
}))

# Proteasome (degradation) control
proteasome <- grep("^PSM[ABCD]", all_genes, value = TRUE)
prot_res <- layer_cor(proteasome)

res <- data.frame(
  set = c("Spliceosome", paste0("Random_", 1:100), "Proteasome"),
  rho_protein = c(splice_res["prot"], rand_res[, "prot"], prot_res["prot"]),
  rho_mrna = c(splice_res["mrna"], rand_res[, "mrna"], prot_res["mrna"]),
  n = c(splice_res["n"], rand_res[, "n"], prot_res["n"])
)
res$discordance <- res$rho_protein - res$rho_mrna
write.csv(res, "results/specificity_control.csv", row.names = FALSE)

sink("results/specificity_summary.txt", split = TRUE)
cat("=== Spliceosome specificity control ===\n\n")
cat(sprintf("Spliceosome:  protein rho = %.2f | mRNA rho = %.2f | discordance = %.2f\n",
            splice_res["prot"], splice_res["mrna"], splice_res["prot"] - splice_res["mrna"]))
cat(sprintf("Random sets:  discordance mean = %.2f (range %.2f to %.2f)\n",
            mean(res$discordance[2:101]), min(res$discordance[2:101]), max(res$discordance[2:101])))
cat(sprintf("Proteasome:   protein rho = %.2f | mRNA rho = %.2f | discordance = %.2f\n",
            prot_res["prot"], prot_res["mrna"], prot_res["prot"] - prot_res["mrna"]))
# empirical p: spliceosome discordance vs random
p_emp <- (sum(res$discordance[2:101] >= res$discordance[1]) + 1) / 101
cat(sprintf("\nEmpirical p (spliceosome discordance >= random): %.3f\n", p_emp))
sink()

# Plot
res$type <- c("Spliceosome", rep("Random", 100), "Proteasome")
p <- ggplot(res, aes(x = type, y = discordance, fill = type)) +
  geom_point(size = 4, shape = 21, color = "black") +
  scale_fill_manual(values = c("Spliceosome" = "#d62728", "Random" = "grey70",
                               "Proteasome" = "#1f77b4")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(y = "Layer discordance (protein rho - mRNA rho vs WEBS)",
       title = "Spliceosome mRNA-protein discordance is not generic",
       caption = sprintf("Empirical p = %.3f (spliceosome vs 10 random sets)", p_emp))

pdf("figures/figS_specificity.pdf", width = 6, height = 5)
print(p)
dev.off()
message("Done.")
