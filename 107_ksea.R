# =============================================================================
# 107. Kinase activity inference (KSEA-style) on KNCC longitudinal phospho
#      Gene-level kinase-substrate sets (OmniPath phospho databases)
#      on raw and parent-protein-corrected site deltas
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

ks <- read.delim("data/raw/kinase_substrate.tsv")
site <- read.csv("results/kncc/phosphosite_paired.csv")

# gene-level deltas: median of per-site deltas (raw and regression-corrected)
site_ok <- site[!is.na(site$med), ]
gene_raw <- tapply(site_ok$med, site_ok$gene, median)
site_okc <- site[!is.na(site$med_reg), ]
gene_cor <- tapply(site_okc$med_reg, site_okc$gene, median)
cat("genes with raw site deltas:", length(gene_raw), "| corrected:", length(gene_cor), "\n")

ks <- ks[ks$substrate %in% names(gene_cor), ]
cat("kinase-substrate edges with corrected coverage:", nrow(ks), "\n")
cat("kinases with >=5 covered substrates:", sum(table(ks$kinase) >= 5), "\n")

kinases <- names(table(ks$kinase)[table(ks$kinase) >= 5])
ksea_test <- function(kin, vals) {
  sub <- intersect(ks$substrate[ks$kinase == kin], names(vals))
  bg <- setdiff(names(vals), sub)
  if (length(sub) < 5) return(c(n = length(sub), stat = NA, p = NA))
  wt <- suppressWarnings(wilcox.test(vals[sub], vals[bg]))
  c(n = length(sub), stat = median(vals[sub]) - median(vals[bg]), p = wt$p.value)
}
rr_raw <- t(sapply(kinases, ksea_test, vals = gene_raw))
rr_cor <- t(sapply(kinases, ksea_test, vals = gene_cor))
res <- data.frame(kinase = kinases,
                  n_sub = rr_raw[, 1],
                  raw_stat = rr_raw[, 2], raw_p = rr_raw[, 3],
                  cor_stat = rr_cor[, 2], cor_p = rr_cor[, 3])
res$raw_q <- p.adjust(res$raw_p, "BH")
res$cor_q <- p.adjust(res$cor_p, "BH")
res <- res[order(res$cor_p), ]
write.csv(res, "results/kncc/ksea_kinases.csv", row.names = FALSE)

sink("results/kncc/ksea.txt", split = TRUE)
cat("=== KSEA-style kinase activity (KNCC paired T2-T1) ===\n")
cat("kinases tested (>=5 covered substrates):", nrow(res), "\n")
cat("RAW deltas: FDR<0.05:", sum(res$raw_q < 0.05, na.rm = TRUE), "\n")
print(res %>% filter(raw_q < 0.05) %>% dplyr::select(kinase, n_sub, raw_stat, raw_q) , row.names = FALSE)
cat("\nCORRECTED (parent-protein regression): FDR<0.05:", sum(res$cor_q < 0.05, na.rm = TRUE), "\n")
print(head(res %>% dplyr::select(kinase, n_sub, cor_stat, cor_p, cor_q), 15), row.names = FALSE)

# single-substrate dominance check for top corrected kinases
cat("\nleave-one-substrate-out for top 5 corrected kinases:\n")
for (kin in head(res$kinase, 5)) {
  sub <- intersect(ks$substrate[ks$kinase == kin], names(gene_cor))
  loo_p <- sapply(sub, function(s) {
    ss <- setdiff(sub, s)
    if (length(ss) < 4) return(NA)
    bg <- setdiff(names(gene_cor), sub)
    suppressWarnings(wilcox.test(gene_cor[ss], gene_cor[bg]))$p.value
  })
  cat(" ", kin, ": LOO p range [", signif(min(loo_p, na.rm = TRUE), 2), ",",
      signif(max(loo_p, na.rm = TRUE), 2), "]",
      if (max(loo_p, na.rm = TRUE) > 0.05) "(single-substrate driven!)" else "(robust)", "\n")
}
sink()
cat("Done.\n")
