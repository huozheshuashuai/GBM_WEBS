# =============================================================================
# 94. KNCC phospho module-level competing-hypothesis screen
#     per-site regression-corrected deltas (site ~ parent protein across patients)
#     -> gene-level mean -> Hallmark module scores -> same adjudication model
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
ph_delta <- b$phos_delta          # sites x 89 patients
ph_gene <- b$phos_gene
prot_delta <- b$prot_delta        # proteins x 105 patients
pts <- intersect(colnames(ph_delta), colnames(prot_delta))
cat("dual-paired patients used:", length(pts), "\n")
ph_d <- ph_delta[, pts]
pd <- prot_delta[, pts]

# ---- per-site regression correction vs parent protein delta ----
parent <- pd[match(ph_gene, rownames(pd)), ]
site_test <- function(v, min_n = 25) {
  v <- v[!is.na(v)]
  if (length(v) < min_n) return(NA_real_)
  median(v)
}
resid_med <- numeric(nrow(ph_d))
for (i in seq_len(nrow(ph_d))) {
  ok <- !is.na(ph_d[i, ]) & !is.na(parent[i, ])
  if (sum(ok) < 25) { resid_med[i] <- NA; next }
  fit <- tryCatch(lm(ph_d[i, ok] ~ parent[i, ok]), error = function(e) NULL)
  resid_med[i] <- if (is.null(fit)) NA else median(resid(fit))
}
gene_corr <- tapply(resid_med, ph_gene, median, na.rm = TRUE)
gene_corr <- gene_corr[!is.na(gene_corr)]
cat("genes with corrected phospho signal:", length(gene_corr), "\n")

# ---- Hallmark module scores on corrected phospho ----
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, names(gene_corr)))
hm <- hm[sapply(hm, length) >= 10]
cat("modules with >=10 phospho-covered genes:", length(hm), "\n")

# predictors from PROTEIN deltas (same adjudication model)
prolif_core <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                           "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(pd))
neuron <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                      "GAP43","TUBB3","MAP2"), rownames(pd))
immune <- intersect(c("CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC",
                      "CSF1R","TMEM119","CX3CR1"), rownames(pd))
d_p0 <- colMeans(pd[prolif_core, ], na.rm = TRUE)
d_n0 <- colMeans(pd[neuron, ], na.rm = TRUE)
d_i0 <- colMeans(pd[immune, ], na.rm = TRUE)

# module score per patient: mean corrected site delta over module genes' sites
# (build per-gene corrected delta matrix first for patient-level module scores)
gene_mat <- rowsum(ifelse(is.na(resid_med), 0, 0), ph_gene)  # placeholder to keep order
# per-patient per-gene corrected delta matrix (median of site residuals per gene)
genes_u <- unique(ph_gene[!is.na(ph_gene)])
gm <- matrix(NA_real_, length(genes_u), length(pts), dimnames = list(genes_u, pts))
resid_full <- matrix(NA_real_, nrow(ph_d), ncol(ph_d))
for (i in seq_len(nrow(ph_d))) {
  ok <- !is.na(ph_d[i, ]) & !is.na(parent[i, ])
  if (sum(ok) < 25) next
  fit <- tryCatch(lm(ph_d[i, ok] ~ parent[i, ok]), error = function(e) NULL)
  if (!is.null(fit)) resid_full[i, ok] <- resid(fit)
}
for (g in genes_u) {
  idx <- which(ph_gene == g)
  gm[g, ] <- colMeans(resid_full[idx, , drop = FALSE], na.rm = TRUE)
}

out <- lapply(names(hm), function(mn) {
  g <- intersect(hm[[mn]], rownames(gm))
  d_m <- colMeans(gm[g, , drop = FALSE], na.rm = TRUE)
  ok <- !is.na(d_m)
  if (sum(ok) < 50) return(NULL)
  fit <- lm(d_m[ok] ~ d_p0[ok] + d_n0[ok] + d_i0[ok])
  cf <- summary(fit)$coefficients
  ci <- confint(fit)
  wt <- suppressWarnings(wilcox.test(d_m[ok], mu = 0))
  data.frame(module = mn, n_genes = length(g),
             raw_median = median(d_m[ok]), raw_p = wt$p.value,
             adj_intercept = cf[1, 1], adj_lo = ci[1, 1], adj_hi = ci[1, 2], adj_p = cf[1, 4],
             r2 = summary(fit)$r.squared)
})
res <- do.call(rbind, out)
res$adj_q <- p.adjust(res$adj_p, "BH")
res$raw_q <- p.adjust(res$raw_p, "BH")
res$class <- ifelse(res$raw_q < 0.05 & res$adj_q < 0.05, "B (independent residual)",
             ifelse(res$raw_q < 0.05, "A (explained by cell state)", "C (n.s.)"))
res <- res[order(res$adj_intercept), ]
write.csv(res, "results/kncc/phospho_module_adjudication.csv", row.names = FALSE)

sink("results/kncc/phospho_module_adjudication.txt", split = TRUE)
cat("=== Phospho module adjudication (regression-corrected, n =", length(pts), "patients) ===\n\n")
cat("class A:", sum(res$class == "A (explained by cell state)"),
    "| B:", sum(res$class == "B (independent residual)"),
    "| C:", sum(res$class == "C (n.s.)"), "\n\n")
print(res[, c("module", "n_genes", "raw_median", "raw_q", "adj_intercept", "adj_q", "r2", "class")],
      row.names = FALSE)
sink()
cat("Done.\n")
