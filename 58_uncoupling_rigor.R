# =============================================================================
# 58. RNA-protein uncoupling rigor:
#     bootstrap CI of delta-rho, leave-one-protein-out, per-gene scatter,
#     formal comparison, and control for proliferation/immune content
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr) })
set.seed(42)
source("08_cptac_proteomics.R")
source("tools/WEBS/webs.R")

proteome <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct",
                                 row.names = 1, check.names = FALSE))
rnaseq <- as.matrix(read.delim("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct",
                               row.names = 1, check.names = FALSE))
webs <- read.csv("results/cptac_webs.csv")
common_samples <- intersect(colnames(proteome), intersect(colnames(rnaseq), webs$sample))
prot_w <- webs$WEBS[match(common_samples, webs$sample)]
rna_webs <- webs_compute_bulk(rnaseq)
mrna_w <- rna_webs$WEBS[match(common_samples, rna_webs$sample)]

splice_in <- intersect(get_spliceosome_genes(),
                       intersect(rownames(proteome), rownames(rnaseq)))
n_s <- length(splice_in)
cat("Spliceosome genes with both layers:", n_s, "\n")

score_p <- colMeans(proteome[splice_in, common_samples], na.rm = TRUE)
score_r <- colMeans(rnaseq[splice_in, common_samples], na.rm = TRUE)
rho_p <- cor(prot_w, score_p, method = "spearman")
rho_m <- cor(mrna_w, score_r, method = "spearman")

# --- per-gene correlations ---
per_gene <- data.frame(
  gene = splice_in,
  rho_protein = sapply(splice_in, function(g)
    cor(prot_w, proteome[g, common_samples], method = "spearman", use = "pairwise.complete.obs")),
  rho_mrna = sapply(splice_in, function(g)
    cor(mrna_w, rnaseq[g, common_samples], method = "spearman", use = "pairwise.complete.obs")))
write.csv(per_gene, "results/uncoupling_per_gene.csv", row.names = FALSE)

# --- bootstrap CI of rho_p, rho_m, delta ---
B <- 2000
n <- length(common_samples)
boot <- t(sapply(seq_len(B), function(i) {
  idx <- sample(n, replace = TRUE)
  c(cor(prot_w[idx], score_p[idx], method = "spearman"),
    cor(mrna_w[idx], score_r[idx], method = "spearman"))
}))
boot_d <- boot[, 1] - boot[, 2]
ci_p <- quantile(boot[, 1], c(0.025, 0.975))
ci_m <- quantile(boot[, 2], c(0.025, 0.975))
ci_d <- quantile(boot_d, c(0.025, 0.975))
p_delta <- 2 * min(mean(boot_d <= 0), mean(boot_d >= 0))

# --- leave-one-protein-out ---
loo <- sapply(seq_len(n_s), function(i) {
  sc <- colMeans(proteome[splice_in[-i], common_samples], na.rm = TRUE)
  cor(prot_w, sc, method = "spearman")
})

# --- partial correlation controlling proliferation & immune protein scores ---
suppressMessages(library(ppcor))
prolif <- colMeans(proteome[intersect(c("MKI67","TOP2A","PCNA","MCM2"), rownames(proteome)),
                            common_samples, drop = FALSE], na.rm = TRUE)
immune <- colMeans(proteome[intersect(c("PTPRC","CD68","CD14"), rownames(proteome)),
                            common_samples, drop = FALSE], na.rm = TRUE)
pc_full <- pcor.test(prot_w, score_p, cbind(prolif, immune), method = "spearman")
pc_mrna <- pcor.test(mrna_w, score_r, cbind(prolif, immune), method = "spearman")

sink("results/uncoupling_rigor.txt", split = TRUE)
cat("=== RNA-protein uncoupling rigor ===\n\n")
cat(sprintf("Pathway-level: rho_protein = %.3f (95%% CI %.3f to %.3f)\n", rho_p, ci_p[1], ci_p[2]))
cat(sprintf("                 rho_mRNA    = %.3f (95%% CI %.3f to %.3f)\n", rho_m, ci_m[1], ci_m[2]))
cat(sprintf("                 delta       = %.3f (95%% CI %.3f to %.3f), bootstrap p = %.4f\n",
            rho_p - rho_m, ci_d[1], ci_d[2], p_delta))
cat(sprintf("\nLeave-one-protein-out: rho_protein range %.3f to %.3f (min drop %.3f)\n",
            min(loo), max(loo), rho_p - min(loo)))
cat(sprintf("\nPartial correlation controlling proliferation + immune protein scores:\n"))
cat(sprintf("  protein: rho = %.3f, p = %.3g\n", pc_full$estimate, pc_full$p.value))
cat(sprintf("  mRNA:    rho = %.3f, p = %.3g\n", pc_mrna$estimate, pc_mrna$p.value))
cat(sprintf("\nPer-gene: protein rho > mRNA rho in %d/%d genes; median protein rho = %.2f, median mRNA rho = %.2f\n",
            sum(per_gene$rho_protein > per_gene$rho_mrna), n_s,
            median(per_gene$rho_protein), median(per_gene$rho_mrna)))
sink()
cat(readLines("results/uncoupling_rigor.txt"), sep = "\n")

# figure: per-gene scatter
p <- ggplot(per_gene, aes(rho_mrna, rho_protein)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, alpha = 0.7, color = "#7B4EA3") +
  annotate("text", x = min(per_gene$rho_mrna), y = max(per_gene$rho_protein),
           hjust = 0, label = paste0(sum(per_gene$rho_protein > per_gene$rho_mrna),
                                     "/", n_s, " proteins above diagonal"),
           color = "#d62728") +
  theme_bw(base_size = 12) +
  labs(title = "Per-gene RNA vs protein correlation with WEBS (spliceosome)",
       x = "mRNA-level rho", y = "Protein-level rho")
ggsave("figures_final/fig_uncoupling_pergene.pdf", p, width = 6.5, height = 5.5)
cat("Done.\n")
