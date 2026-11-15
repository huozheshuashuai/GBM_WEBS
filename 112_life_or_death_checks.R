# =============================================================================
# 112. Life-or-death checks (advisor V5):
#      A. DIRECT module-level technical variance from duplicate aliquots
#         (module score on each duplicate copy; dup-diff variance estimates
#          delta technical variance incl. plex effects; CI + LODO)
#      B. Module overlap audit for co-evolution structure
#         (Jaccard matrix, shared-protein-removed correlations, strict
#          non-overlapping modules PCA, leave-one-patient-out PCA)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

# ---------- load raw proteome keeping duplicate aliquots separate ----------
f <- "E:/KIMI/data/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log KNCC_GBM", hdr))
cc <- rep("NULL", length(hdr)); cc[1] <- "character"; cc[log_idx] <- "numeric"
pr <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc, check.names = FALSE)
pr <- pr[!pr$Gene %in% c("Mean","Median","StdDev") & !is.na(pr$Gene), ]
raw_mat <- as.matrix(pr[, -1]); rownames(raw_mat) <- pr$Gene
colnames(raw_mat) <- sub(":.*", "", sub("^Log ", "", colnames(raw_mat)))

dup_names <- unique(colnames(raw_mat)[duplicated(colnames(raw_mat))])
cat("duplicate aliquot pairs:", length(dup_names), "\n")

# frozen 7 modules
b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(raw_mat)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(raw_mat)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(raw_mat)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)

# ---------- A. direct module-level technical variance ----------
# module score on each duplicate copy; diff between copies = technical delta error
# (spanning plexes, so it includes plex effects -> conservative overestimate)
tech <- lapply(dup_names, function(s) {
  idx <- which(colnames(raw_mat) == s)
  m1 <- raw_mat[, idx[1]]; m2 <- raw_mat[, idx[2]]
  sapply(mods, function(g) mean(m1[g] - m2[g], na.rm = TRUE))
})
tech_mat <- do.call(cbind, tech)   # modules x duplicate-pairs
tech_var_mod <- apply(tech_mat, 1, var, na.rm = TRUE)

# patient delta variance per module (biological + technical)
pat_var <- sapply(mods, function(g) var(colMeans(delta[g, , drop = FALSE], na.rm = TRUE), na.rm = TRUE))

# ratio with bootstrap CI over duplicate pairs AND patients
B <- 2000
boot_ratio <- array(NA, c(B, length(mods)), dimnames = list(NULL, names(mods)))
ndup <- ncol(tech_mat); npts <- ncol(delta)
pd_mat <- sapply(mods, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE))
for (i in seq_len(B)) {
  bi <- sample(ndup, ndup, TRUE); bj <- sample(npts, npts, TRUE)
  tv <- apply(tech_mat[, bi, drop = FALSE], 1, var, na.rm = TRUE)
  pv <- apply(pd_mat[bj, , drop = FALSE], 2, var, na.rm = TRUE)
  boot_ratio[i, ] <- pv / tv
}
ratio_ci <- apply(boot_ratio, 2, quantile, c(0.025, 0.975), na.rm = TRUE)

sink("results/kncc/techvariance_direct.txt", split = TRUE)
cat("=== A. DIRECT module-level technical variance (V5 check 1) ===\n\n")
cat("estimator: variance of module-score differences between duplicate aliquots\n")
cat("(duplicates span plexes -> includes plex effects -> conservative, i.e.\n")
cat(" overestimates pure measurement error -> ratios are conservative)\n")
cat("delta technical variance of T2-T1 = variance of duplicate-pair difference\n")
cat("(two independent measurements, same logic as two timepoints)\n\n")
res <- data.frame(module = names(mods),
                  n_proteins = sapply(mods, length),
                  tech_var_delta = round(tech_var_mod, 5),
                  patient_var_delta = round(pat_var, 4),
                  ratio = round(pat_var / tech_var_mod, 1),
                  ratio_ci_lo = round(ratio_ci[1, ], 1),
                  ratio_ci_hi = round(ratio_ci[2, ], 1))
print(res, row.names = FALSE)

# leave-one-duplicate-out
cat("\nleave-one-duplicate-out ratio range:\n")
for (j in seq_len(ndup)) {
  tv <- apply(tech_mat[, -j, drop = FALSE], 1, var, na.rm = TRUE)
  cat("  drop", dup_names[j], ": ratio range [",
      round(min(pat_var / tv), 1), ",", round(max(pat_var / tv), 1), "]\n")
}

# single-protein sanity: single-protein delta tech var vs module
cat("\nsingle-protein delta technical variance (dup-diff):",
    round(mean(sapply(dup_names, function(s) {
      idx <- which(colnames(raw_mat) == s)
      var(raw_mat[, idx[1]] - raw_mat[, idx[2]], na.rm = TRUE)
    })), 4), "\n")
sink()

# ---------- B. module overlap audit ----------
all23 <- c(hm[c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
                "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
                "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
                "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")],
           list(NEURONAL_CURATED = mods$NEURO, PROLIF_CURATED = mods$PROLIF))
all23 <- all23[sapply(all23, length) >= 8]

# Jaccard matrix
nm <- names(all23)
jac <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
for (i in seq_along(nm)) for (j in seq_along(nm))
  jac[i, j] <- length(intersect(all23[[i]], all23[[j]])) / length(union(all23[[i]], all23[[j]]))

sink("results/kncc/techvariance_direct.txt", append = TRUE, split = TRUE)
cat("\n=== B. module overlap audit (V5 check 2) ===\n\n")
cat("max Jaccard between module pairs:", round(max(jac[upper.tri(jac)]), 3), "\n")
top_overlap <- which(jac > 0.05 & upper.tri(jac), arr.ind = TRUE)
if (nrow(top_overlap)) {
  for (k in seq_len(nrow(top_overlap)))
    cat("  ", nm[top_overlap[k, 1]], "~", nm[top_overlap[k, 2]],
        "Jaccard =", round(jac[top_overlap[k, 1], top_overlap[k, 2]], 3), "\n")
} else cat("  (no pair with Jaccard > 0.05)\n")

# delta correlations with shared proteins removed (for top correlated pairs)
Md <- t(sapply(all23, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
cc_full <- cor(t(Md), method = "spearman")
pairs_grid <- t(combn(nm, 2))
pc <- data.frame(m1 = pairs_grid[, 1], m2 = pairs_grid[, 2],
                 rho_full = mapply(function(a, b) cc_full[a, b], pairs_grid[, 1], pairs_grid[, 2]))
pc <- pc[order(-pc$rho_full), ]
cat("\ntop pairs: full vs shared-proteins-removed delta correlation:\n")
for (k in 1:min(8, nrow(pc))) {
  a <- pc$m1[k]; bb <- pc$m2[k]
  ga <- setdiff(all23[[a]], all23[[bb]]); gb <- setdiff(all23[[bb]], all23[[a]])
  if (length(ga) < 5 || length(gb) < 5) { cat("  ", a, "~", bb, ": too few unique proteins\n"); next }
  da <- colMeans(delta[ga, , drop = FALSE], na.rm = TRUE)
  db <- colMeans(delta[gb, , drop = FALSE], na.rm = TRUE)
  ok <- !is.na(da) & !is.na(db)
  rho_u <- cor(da[ok], db[ok], method = "spearman")
  cat("  ", a, "~", bb, ": full rho =", round(pc$rho_full[k], 3),
      " -> unique-only rho =", round(rho_u, 3),
      " (shared:", length(intersect(all23[[a]], all23[[bb]])), ")\n")
}

# strict PCA on non-overlapping module definitions
# greedy dedup: for each shared protein assign to the module where it has higher |loading| proxy
# simpler: iteratively remove proteins appearing in >1 module from all but the largest module
mod_strict <- all23
seen <- character(0)
for (mn in nm[order(-sapply(all23, length))]) {
  mod_strict[[mn]] <- setdiff(all23[[mn]], seen)
  seen <- c(seen, all23[[mn]])
}
mod_strict <- mod_strict[sapply(mod_strict, length) >= 8]
cat("\nstrict non-overlapping modules retained:", length(mod_strict), "\n")
Ms <- t(sapply(mod_strict, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)))
Ms <- Ms[, colSums(is.na(Ms)) == 0]
pc_s <- prcomp(t(Ms), scale. = TRUE)
ve_s <- pc_s$sdev^2 / sum(pc_s$sdev^2)
cat("strict-module PCA: PC1 =", round(ve_s[1] * 100, 1), "% PC2 =", round(ve_s[2] * 100, 1), "%\n")

# leave-one-patient-out influence on PC1/PC2 variance shares
pc_full <- prcomp(t(Md), scale. = TRUE)
ve_full <- pc_full$sdev^2 / sum(pc_full$sdev^2)
loo_ve <- t(sapply(seq_len(ncol(Md)), function(j) {
  v <- prcomp(t(Md[, -j]), scale. = TRUE)$sdev
  v^2 / sum(v^2)
}))
cat("LOPO PC1 variance share range: [", round(min(loo_ve[, 1]) * 100, 1), ",",
    round(max(loo_ve[, 1]) * 100, 1), "]% | PC2 range: [", round(min(loo_ve[, 2]) * 100, 1), ",",
    round(max(loo_ve[, 2]) * 100, 1), "]%\n")
cat("(PCA used correlation matrix, scale. = TRUE)\n")
sink()

# figures
res$module <- factor(res$module, levels = res$module[order(-res$ratio)])
p1 <- ggplot(res, aes(module, ratio, ymin = ratio_ci_lo, ymax = ratio_ci_hi)) +
  geom_pointrange(color = "#7B4EA3") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() + theme_bw(base_size = 11) +
  labs(y = "Patient delta variance / duplicate-estimated technical delta variance",
       x = NULL, title = "Heterogeneity vs DIRECT module-level technical variance (95% CI)")
ggsave("figures_kncc/techvariance_direct.pdf", p1, width = 7.5, height = 5)

pdf("figures_kncc/module_jaccard.pdf", width = 8, height = 7)
pheatmap::pheatmap(jac, breaks = seq(0, max(jac) + 0.01, length.out = 101),
                   color = colorRampPalette(c("white", "#3D9B8C"))(100),
                   main = "Module overlap (Jaccard index)")
dev.off()
cat("Done.\n")
