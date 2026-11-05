# =============================================================================
# 28b. Specificity control UPGRADE: 10,000 expression-matched random gene sets
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
source("tools/WEBS/webs.R")
rna_webs <- webs_compute_bulk(rnaseq)
mrna_w <- rna_webs$WEBS[match(common_samples, rna_webs$sample)]

# expression strata for matching (proteome mean abundance deciles)
all_genes <- intersect(rownames(proteome), rownames(rnaseq))
expr_mean <- rowMeans(proteome[all_genes, common_samples], na.rm = TRUE)
decile <- cut(expr_mean, breaks = quantile(expr_mean, probs = 0:10/10, na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
names(decile) <- all_genes

splice_genes <- get_spliceosome_genes()
splice_in <- intersect(splice_genes, all_genes)

layer_cor <- function(genes) {
  score_p <- colMeans(proteome[genes, common_samples, drop = FALSE], na.rm = TRUE)
  score_r <- colMeans(rnaseq[genes, common_samples, drop = FALSE], na.rm = TRUE)
  c(prot = cor(prot_w, score_p, method = "spearman", use = "pairwise.complete.obs"),
    mrna = cor(mrna_w, score_r, method = "spearman", use = "pairwise.complete.obs"))
}
splice_res <- layer_cor(splice_in)
splice_disc <- splice_res["prot"] - splice_res["mrna"]

# 10,000 expression-matched random sets: for each spliceosome gene, sample a gene
# from the same proteome-abundance decile
B <- 10000
by_decile <- split(all_genes, decile)
splice_dec <- decile[splice_in]
rand_disc <- numeric(B)
for (i in seq_len(B)) {
  gs <- sapply(seq_along(splice_in), function(j) {
    pool <- by_decile[[as.character(splice_dec[j])]]
    sample(pool, 1)
  })
  r <- layer_cor(gs)
  rand_disc[i] <- r["prot"] - r["mrna"]
}
p_emp <- (sum(rand_disc >= splice_disc) + 1) / (B + 1)

sink("results/specificity_10k_summary.txt", split = TRUE)
cat("=== Specificity control: 10,000 expression-matched random gene sets ===\n\n")
cat(sprintf("Spliceosome discordance: %.3f (protein rho = %.2f, mRNA rho = %.2f)\n",
            splice_disc, splice_res["prot"], splice_res["mrna"]))
cat(sprintf("Random (expression-matched): mean = %.3f, SD = %.3f, max = %.3f\n",
            mean(rand_disc), sd(rand_disc), max(rand_disc)))
cat(sprintf("Empirical p = %.5f  (n = %d permutations)\n", p_emp, B))
cat(sprintf("Z-score = %.2f\n", (splice_disc - mean(rand_disc)) / sd(rand_disc)))
sink()
cat(readLines("results/specificity_10k_summary.txt"), sep = "\n")

df <- data.frame(discordance = rand_disc)
p <- ggplot(df, aes(discordance)) +
  geom_histogram(bins = 60, fill = "grey75", color = "white") +
  geom_vline(xintercept = splice_disc, color = "#d62728", linewidth = 1.2) +
  annotate("text", x = splice_disc, y = Inf, vjust = 1.5, hjust = -0.1,
           label = paste0("Spliceosome (p = ", signif(p_emp, 2), ")"), color = "#d62728") +
  theme_bw(base_size = 12) +
  labs(title = "Spliceosome layer discordance vs 10,000 expression-matched random gene sets",
       x = "Layer discordance (protein rho - mRNA rho vs WEBS)", y = "Count")
ggsave("figures_final/figS_specificity_10k.pdf", p, width = 7, height = 5)
cat("Done.\n")
