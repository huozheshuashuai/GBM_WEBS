# =============================================================================
# 41. SNF (Similarity Network Fusion) — true multi-omics integration on CPTAC
# =============================================================================
# Views: proteome / acetylome (site->gene mean) / phosphoproteome / RNA-seq
# Question: do consensus multi-omic clusters recapitulate WEBS groups?
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(SNFtool); library(dplyr); library(ggplot2) })
set.seed(42)

webs <- read.csv("results/cptac_webs.csv")

load_view <- function(f, collapse = FALSE) {
  m <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  m <- m[rowSums(!is.na(m)) > ncol(m) * 0.5, ]
  if (collapse) {
    g <- sub("[-_\\.][A-Za-z]*[0-9]+.*$", "", rownames(m))
    m <- rowsum(m, g, reorder = FALSE, na.rm = TRUE) /
         pmax(1, as.vector(table(g)[rownames(rowsum(m, g, reorder = FALSE))]))
  }
  m
}

prot <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct")
acet <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_acetylome_mssm_site.cct", collapse = TRUE)
phos <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_phosphoproteome_mssm_site.cct", collapse = TRUE)
rna  <- load_view("data/raw/CPTAC/HS_CPTAC_GBM_rnaseq_fpkm_uq_log2.cct")

topv <- function(m, n = 2000) m[order(apply(m, 1, var, na.rm = TRUE), decreasing = TRUE)[1:min(n, nrow(m))], ]
views <- list(proteome = topv(prot), acetylome = topv(acet),
              phospho = topv(phos), rnaseq = topv(rna))
samples <- Reduce(intersect, lapply(views, colnames))
views <- lapply(views, function(m) { m <- m[, samples]; m[complete.cases(m), ] })
cat("Views x", length(samples), "samples\n")

# affinity matrices + SNF
make_W <- function(m) {
  d <- dist(scale(t(m)))
  affinityMatrix(as.matrix(d), K = 20, sigma = 0.5)
}
Ws <- lapply(views, make_W)
W <- SNF(Ws, K = 20, t = 20)

# spectral clustering into 2-4 clusters, pick best by eigengap
est <- tryCatch(estimateNumberOfClustersGivenGraph(W, NUMC = 2:5), error = function(e) NULL)
k_best <- 3
if (!is.null(est)) {
  vals <- sapply(est, function(x) x[1])
  nm <- names(vals)
  if (!is.null(nm) && !any(is.na(vals))) {
    k_best <- as.integer(nm[which.max(vals)])
    if (is.na(k_best)) k_best <- 3
  }
}
cat("SNF clusters:", k_best, "\n")
clusters <- spectralClustering(W, K = k_best)

# association with WEBS groups
snf_df <- data.frame(sample = samples, cluster = factor(clusters))
snf_df$WEBS_group <- webs$WEBS_group[match(samples, webs$sample)]
snf_df$WEBS <- webs$WEBS[match(samples, webs$sample)]

tab <- table(snf_df$cluster, snf_df$WEBS_group)
chi <- chisq.test(tab)
kw <- kruskal.test(WEBS ~ cluster, data = snf_df)

sink("results/snf_summary.txt", split = TRUE)
cat("=== SNF multi-omics integration ===\n")
cat("Optimal clusters:", k_best, "\n\n")
cat("Cluster x WEBS_group contingency:\n"); print(tab)
cat("Chi-square p:", signif(chi$p.value, 3), "\n")
cat("KW WEBS ~ cluster p:", signif(kw$p.value, 3), "\n")
sink()

write.csv(snf_df, "results/snf_clusters.csv", row.names = FALSE)

# heatmap of fused network
ord <- order(clusters)
pdf("figures_final/fig_snf_heatmap.pdf", width = 6.5, height = 6)
image(W[ord, ord], col = colorRampPalette(c("white", "#6B4C9A"))(50),
      axes = FALSE, main = sprintf("SNF fused similarity (k = %d, chi2 p = %.2g)", k_best, chi$p.value))
dev.off()

# WEBS by cluster boxplot
p <- ggplot(snf_df, aes(x = cluster, y = WEBS, fill = cluster)) +
  geom_boxplot(outlier.size = 0.8) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = sprintf("SNF consensus clusters vs WEBS (KW p = %.2g)", kw$p.value),
       x = "SNF cluster", y = "Protein-level WEBS")
ggsave("figures_final/fig_snf_webs.pdf", p, width = 5.5, height = 5)
cat("Done.\n")
