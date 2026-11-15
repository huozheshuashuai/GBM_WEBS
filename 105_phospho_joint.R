# =============================================================================
# 105. Phospho effective-coverage audit + patient-level protein-phospho
#      joint trajectory analysis
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
ph_delta <- b$phos_delta; ph_gene <- b$phos_gene
prot_delta <- b$prot_delta
pts <- intersect(colnames(ph_delta), colnames(prot_delta))
ph_d <- ph_delta[, pts]; pd <- prot_delta[, pts]
cat("dual-paired patients:", length(pts), "\n")

# ---- 1. effective coverage audit ----
pairs_per_site <- rowSums(!is.na(ph_d))
sink("results/kncc/phospho_coverage.txt", split = TRUE)
cat("=== phospho effective coverage audit ===\n")
cat("total sites:", nrow(ph_d), "\n")
cat("sites with >=25 pairs:", sum(pairs_per_site >= 25), "\n")
cat("sites with >=40 pairs:", sum(pairs_per_site >= 40), "\n")
cat("sites with >=60 pairs:", sum(pairs_per_site >= 60), "\n")
cat("median pairs per site:", median(pairs_per_site), "\n")
sink()

# ---- 2. per-site regression-corrected deltas (recompute; as in 90/94) ----
keep <- pairs_per_site >= 40
ph_k <- ph_d[keep, ]
parent <- pd[match(ph_gene[keep], rownames(pd)), ]
cat("sites with >=40 pairs:", nrow(ph_k), "\n")

resid_mat <- matrix(NA_real_, nrow(ph_k), ncol(ph_k), dimnames = dimnames(ph_k))
for (i in seq_len(nrow(ph_k))) {
  ok <- !is.na(ph_k[i, ]) & !is.na(parent[i, ])
  if (sum(ok) < 40) next
  fit <- tryCatch(lm(ph_k[i, ok] ~ parent[i, ok]), error = function(e) NULL)
  if (!is.null(fit)) resid_mat[i, ok] <- resid(fit)
}
valid <- rowSums(!is.na(resid_mat)) >= 40
resid_mat <- resid_mat[valid, ]
genes_v <- ph_gene[keep][valid]
cat("sites passing correction:", nrow(resid_mat), "\n")

# ---- 3. gene-level corrected signal -> module joint trajectories ----
# per-patient per-gene corrected phospho delta
genes_u <- unique(genes_v)
gm <- matrix(NA_real_, length(genes_u), ncol(resid_mat), dimnames = list(genes_u, colnames(resid_mat)))
for (g in genes_u) {
  idx <- which(genes_v == g)
  gm[g, ] <- colMeans(resid_mat[idx, , drop = FALSE], na.rm = TRUE)
}
gm[gm == "NaN"] <- NA

# Hallmark modules on corrected phospho gene matrix
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(gm)))
hm <- hm[sapply(hm, length) >= 8]

ph_mod <- t(sapply(hm, function(g) colMeans(gm[g, , drop = FALSE], na.rm = TRUE)))
pr_mod <- t(sapply(hm, function(g) colMeans(pd[intersect(g, rownames(pd)), , drop = FALSE], na.rm = TRUE)))

# ---- 4. protein-phospho discordance at module level per patient ----
cc <- sapply(rownames(ph_mod), function(mn)
  cor(ph_mod[mn, ], pr_mod[mn, ], use = "pairwise.complete.obs", method = "spearman"))
sink("results/kncc/phospho_joint.txt", split = TRUE)
cat("\n=== module-level protein vs corrected-phospho delta correlation ===\n")
print(round(sort(cc), 2))
cat("\nmodules where phospho moves but protein does not (or vice versa):\n")
# per-module: phospho delta sig but protein not, patient-level
res <- lapply(rownames(ph_mod), function(mn) {
  x <- ph_mod[mn, ]; y <- pr_mod[mn, ]
  ok <- !is.na(x) & !is.na(y)
  wtp <- suppressWarnings(wilcox.test(x[ok], mu = 0)$p.value)
  wtr <- suppressWarnings(wilcox.test(y[ok], mu = 0)$p.value)
  data.frame(module = mn, phos_med = median(x[ok]), prot_med = median(y[ok]),
             phos_p = wtp, prot_p = wtr)
})
res <- do.call(rbind, res)
res$phos_q <- p.adjust(res$phos_p, "BH")
res$prot_q <- p.adjust(res$prot_p, "BH")
res$quadrant <- ifelse(res$phos_q < 0.05 & res$prot_q >= 0.05, "phospho-only",
                ifelse(res$phos_q >= 0.05 & res$prot_q < 0.05, "protein-only",
                ifelse(res$phos_q < 0.05 & res$prot_q < 0.05, "both", "neither")))
print(res[order(res$phos_p), ], row.names = FALSE)
sink()

# figure: module protein vs phospho medians
res$lab <- ifelse(res$quadrant == "phospho-only" | res$quadrant == "protein-only",
                  sub("HALLMARK_", "", res$module), NA)
p <- ggplot(res, aes(prot_med, phos_med, color = quadrant)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.8, max.overlaps = 30) +
  scale_color_manual(values = c("both" = "grey50", "neither" = "grey85",
                                "phospho-only" = "#7B4EA3", "protein-only" = "#3D9B8C")) +
  theme_bw(base_size = 11) +
  labs(x = "Protein module delta (median)", y = "Corrected phospho module delta (median)",
       title = "Protein vs phosphosignaling evolution (module level)",
       color = "Quadrant")
ggsave("figures_kncc/phospho_joint_quadrant.pdf", p, width = 7, height = 5.5)
saveRDS(list(resid_mat = resid_mat, genes_v = genes_v, gm = gm, ph_mod = ph_mod,
             pr_mod = pr_mod, pts = pts), "results/kncc/phospho_corrected_bundle.rds")
cat("Done.\n")
