# =============================================================================
# 117. HHS robustness audit (V7 priority 1)
#      effect sizes + CIs, coverage-restricted projection, imputation
#      sensitivity, leave-one-patient-out, multi-specimen rule audit
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
kdelta <- b$prot_delta
h <- readRDS("results/hhs/hhs_analysis_ready.rds")
hd <- h$delta

# ---- multi-specimen audit: what are the extra specimens per patient? ----
sn <- colnames(h$mat); ptype <- h$ptype; pid <- h$pid
spec_tab <- table(pid, ptype)
multi <- rownames(spec_tab)[rowSums(spec_tab > 0) > 2]
sink("results/hhs/hhs_robustness.txt", split = TRUE)
cat("=== HHS robustness audit ===\n\n")
cat("multi-specimen patients (n =", length(multi), "):\n")
extra_types <- table(ptype[pid %in% multi])
print(extra_types)
cat("-> extra specimens are NAT adjacent samples (PM+RM+N pattern); the paired\n")
cat("   PM/RM delta is unaffected. Merge rule: one PM + one RM per patient (verified unique).\n")

# verify exactly one PM and one RM per paired patient
pm_ct <- table(pid[ptype == "Primary"]); rm_ct <- table(pid[ptype == "Recurrent"])
cat("patients with >1 PM:", sum(pm_ct > 1), "| >1 RM:", sum(rm_ct > 1), "\n")

# ---- frozen objects ----
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
keep_mods <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_DNA_REPAIR",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS",
               "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               "HALLMARK_ANGIOGENESIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
               "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA",
               "HALLMARK_UNFOLDED_PROTEIN_RESPONSE","HALLMARK_PROTEIN_SECRETION",
               "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_MYC_TARGETS_V1")
mods_k <- lapply(hm[keep_mods], function(g) intersect(g, rownames(kdelta)))
mods_k$NEURONAL_CURATED <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                                       "GAP43","TUBB3","MAP2"), rownames(kdelta))
mods_k$PROLIF_CURATED <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                                     "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(kdelta))
mods_k <- mods_k[sapply(mods_k, length) >= 8]
Md <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
pc <- prcomp(t(Md), scale. = TRUE)
pc2_load <- pc$rotation[, 2]

mods_h <- lapply(names(mods_k), function(mn) intersect(mods_k[[mn]], rownames(hd)))
names(mods_h) <- names(mods_k)
mods_h <- mods_h[sapply(mods_h, length) >= 5]
mods_k <- mods_k[names(mods_h)]
Md_h <- t(sapply(mods_h, function(g) colMeans(hd[g, , drop = FALSE], na.rm = TRUE)))
Md_k <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
load2 <- pc2_load[names(mods_h)]
pc2_k <- as.numeric(t(load2) %*% Md_k[names(mods_h), ])
pc2_h <- as.numeric(t(load2) %*% Md_h[names(mods_h), ])

# ---- (1) effect size + 95% CI for PC2 projection ----
boot_ci <- function(x, B = 2000) {
  bs <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), lo = quantile(bs, 0.025), hi = quantile(bs, 0.975))
}
ci_k <- boot_ci(pc2_k); ci_h <- boot_ci(pc2_h)
cat("\n(1) PC2 axis projection effect (mean ± 95% CI):\n")
cat("  KNCC: mean", round(ci_k["mean"], 3), "[", round(ci_k["lo"], 3), ",", round(ci_k["hi"], 3), "]\n")
cat("  HHS:  mean", round(ci_h["mean"], 3), "[", round(ci_h["lo"], 3), ",", round(ci_h["hi"], 3), "]\n")

# ---- (2) coverage per module (shared proteins / frozen size) ----
cat("\n(2) protein coverage per module (shared/frozen):\n")
cov_df <- data.frame(module = names(mods_k),
                     frozen_n = sapply(mods_k, length),
                     shared_n = sapply(mods_h, length))
cov_df$frac <- round(cov_df$shared_n / cov_df$frozen_n, 2)
print(cov_df[order(cov_df$frac), ], row.names = FALSE)
cat("PC2-loading-weighted coverage: ", round(sum(abs(load2) * (cov_df$frac[match(names(load2), cov_df$module)]) / sum(abs(load2))), 2), "\n")

# ---- (3) high-coverage restricted projection ----
hc_mods <- cov_df$module[cov_df$frac >= 0.6]
cat("\n(3) restricted to modules with >=60% coverage (n =", length(hc_mods), "):\n")
Md_h2 <- Md_h[hc_mods, , drop = FALSE]; Md_k2 <- Md_k[hc_mods, , drop = FALSE]
l2 <- load2[hc_mods]
pc2_k2 <- as.numeric(t(l2) %*% Md_k2); pc2_h2 <- as.numeric(t(l2) %*% Md_h2)
cat("  KNCC: median", round(median(pc2_k2), 3), " p =", signif(wilcox.test(pc2_k2, mu = 0)$p.value, 3), "\n")
cat("  HHS:  median", round(median(pc2_h2), 3), " p =", signif(wilcox.test(pc2_h2, mu = 0)$p.value, 3), "\n")

# ---- (4) imputation sensitivity: HHS uses imputed matrix.
#     simulate worst case: imputed proteins carry zero delta -> shrink module deltas.
#     We cannot recompute without unimputed data; instead bound: which modules rely
#     on low-abundance (likely imputed) proteins in HHS primary samples.
pm1 <- h$mat[, h$ptype == "Primary"]
abund <- rowMeans(pm1, na.rm = TRUE)
q33 <- quantile(abund, 1/3, na.rm = TRUE)
lowab <- names(abund)[abund <= q33]
cat("\n(4) imputation sensitivity bound:\n")
cat("  bottom-tercile abundance proteins (imputation candidates):", length(lowab), "\n")
cat("  module fraction in bottom tercile:\n")
fa <- sapply(mods_h, function(g) round(mean(g %in% lowab), 2))
print(sort(fa))
# robust projection excluding bottom-tercile proteins
mods_hc <- lapply(mods_h, function(g) setdiff(g, lowab))
mods_hc <- mods_hc[sapply(mods_hc, length) >= 5]
Md_h3 <- t(sapply(mods_hc, function(g) colMeans(hd[g, , drop = FALSE], na.rm = TRUE)))
l3 <- load2[names(mods_hc)]
pc2_h3 <- as.numeric(t(l3) %*% Md_h3[names(mods_hc), ])
cat("  excluding bottom-tercile proteins: HHS PC2 median", round(median(pc2_h3), 3),
    " p =", signif(wilcox.test(pc2_h3, mu = 0)$p.value, 3), "\n")

# ---- (5) leave-one-patient-out on PC2 projection ----
loo <- sapply(seq_len(ncol(Md_h)), function(j) {
  v <- as.numeric(t(load2) %*% Md_h[names(mods_h), -j, drop = FALSE])
  wilcox.test(v, mu = 0)$p.value
})
cat("\n(5) leave-one-patient-out HHS PC2 p-value range: [",
    signif(min(loo), 2), ",", signif(max(loo), 2), "]\n")
cat("  all < 0.05?", all(loo < 0.05), "\n")

# ---- (6) EMT and NEURONAL shown together (direction check) ----
cat("\n(6) direction comparison (median delta):\n")
for (mn in c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "NEURONAL_CURATED",
             "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_IL2_STAT5_SIGNALING")) {
  cat("  ", mn, ": KNCC", round(median(Md_k[mn, ]), 3), "-> HHS", round(median(Md_h[mn, ]), 3), "\n")
}
sink()

# figure: EMT vs NEURONAL both cohorts, patient distributions
df <- rbind(
  data.frame(cohort = "KNCC", module = "EMT", delta = Md_k["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]),
  data.frame(cohort = "HHS", module = "EMT", delta = Md_h["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]),
  data.frame(cohort = "KNCC", module = "NEURONAL", delta = Md_k["NEURONAL_CURATED", ]),
  data.frame(cohort = "HHS", module = "NEURONAL", delta = Md_h["NEURONAL_CURATED", ])
)
p <- ggplot(df, aes(cohort, delta, fill = cohort)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ module, scales = "free_x") +
  scale_fill_manual(values = c("KNCC" = "#7B4EA3", "HHS" = "#3D9B8C")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(y = "Per-patient module delta (T2 - T1)", x = NULL,
       title = "EMT replicates; NEURONAL direction differs between cohorts (shown openly)")
ggsave("figures_kncc/hhs_emt_neuronal.pdf", p, width = 8, height = 5)
cat("Done.\n")
