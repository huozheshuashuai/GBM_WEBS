# =============================================================================
# 122. HHS patient-level co-variation structure validation
#      (distinguish: mean-direction replication vs patient-structure replication)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
h <- readRDS("results/hhs/hhs_analysis_ready.rds")
kdelta <- b$prot_delta; hd <- h$delta

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
mods_h <- lapply(names(mods_k), function(mn) intersect(mods_k[[mn]], rownames(hd)))
names(mods_h) <- names(mods_k)
mods_h <- mods_h[sapply(mods_h, length) >= 5]
mods_k <- mods_k[names(mods_h)]

Md_k <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
Md_h <- t(sapply(mods_h, function(g) colMeans(hd[g, , drop = FALSE], na.rm = TRUE)))

# module-delta co-variation matrices (per-patient Spearman)
cc_k <- cor(t(Md_k), method = "spearman")
cc_h <- cor(t(Md_h), method = "spearman")

pairs <- t(combn(names(mods_h), 2))
cmp <- data.frame(m1 = pairs[, 1], m2 = pairs[, 2],
                  kncc = mapply(function(a, b2) cc_k[a, b2], pairs[, 1], pairs[, 2]),
                  hhs  = mapply(function(a, b2) cc_h[a, b2], pairs[, 1], pairs[, 2]))
ct <- cor.test(cmp$kncc, cmp$hhs, method = "spearman")

# bootstrap CI of the structure-consistency correlation (resample HHS patients)
B <- 2000; n_h <- ncol(Md_h)
bs <- replicate(B, {
  idx <- sample(n_h, n_h, TRUE)
  ch <- cor(t(Md_h[, idx]), method = "spearman")
  cor(mapply(function(a, b2) cc_k[a, b2], pairs[, 1], pairs[, 2]),
      mapply(function(a, b2) ch[a, b2], pairs[, 1], pairs[, 2]), method = "spearman")
})
ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)

sink("results/hhs/hhs_structure_validation.txt", split = TRUE)
cat("=== HHS patient-level co-variation structure validation ===\n\n")
cat("question: does the KNCC module co-variation STRUCTURE (not just mean\n")
cat("direction) replicate in HHS longitudinal deltas (n = 40 patients)?\n\n")
cat("module-pair delta correlations, KNCC vs HHS: rho =", round(ct$estimate, 3),
    " p =", signif(ct$p.value, 3), " (", nrow(cmp), " pairs)\n")
cat("HHS-patient bootstrap 95% CI of structure consistency: [",
    round(ci[1], 3), ",", round(ci[2], 3), "]\n\n")
cat("key second-axis pairs:\n")
key <- cmp %>% filter((grepl("IL2_STAT5|REACTIVE_OXYGEN|EPITHELIAL|HYPOXIA", m1) |
                        grepl("IL2_STAT5|REACTIVE_OXYGEN|EPITHELIAL|HYPOXIA", m2))) %>%
  arrange(desc(abs(kncc)))
print(head(key, 10), row.names = FALSE)
sink()

# figure
p <- ggplot(cmp, aes(kncc, hhs)) +
  geom_point(color = "#7B4EA3", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  theme_bw(base_size = 12) +
  labs(x = "KNCC module-pair delta correlation (n = 105)",
       y = "HHS module-pair delta correlation (n = 40)",
       title = sprintf("Patient-level co-variation structure replicates longitudinally (rho = %.2f)",
                       ct$estimate),
       subtitle = sprintf("bootstrap 95%% CI [%.2f, %.2f]", ci[1], ci[2]))
ggsave("figures_kncc/hhs_structure_validation.pdf", p, width = 6, height = 5.5)
write.csv(cmp, "results/hhs/hhs_structure_pairs.csv", row.names = FALSE)
cat("Done.\n")
