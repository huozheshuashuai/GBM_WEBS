# =============================================================================
# 111. Task 9 (V4.0): CPTAC matched cross-sectional validation
#      Validate ONLY the KNCC-discovered structures (module axes, co-evolution
#      blocks, class-B modules) in the independent CPTAC-GBM proteome (n=99)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
pm <- pm[!rownames(pm) %in% c("Mean","Median","StdDev"), ]

# frozen modules (same definitions as 109)
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(pm)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(pm)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(pm)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)
cat("module sizes in CPTAC:", paste(sapply(mods, length), collapse = ", "), "\n")

S <- t(sapply(mods, function(g) colMeans(pm[g, , drop = FALSE], na.rm = TRUE)))

# KNCC reference: co-evolution correlation matrix from script 104 (recompute)
b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta
mods_k <- lapply(mods, function(g) intersect(g, rownames(delta)))
Md <- t(sapply(mods_k, function(g) colMeans(delta[g, , drop = FALSE], na.rm = TRUE)))

kncc_cor <- cor(t(Md), method = "spearman")
cptac_cor <- cor(t(S), method = "spearman")

sink("results/kncc/cptac_matched_validation.txt", split = TRUE)
cat("=== CPTAC matched cross-sectional validation ===\n\n")
cat("note: CPTAC is cross-sectional (99 primary tumors); it can validate\n")
cat("inter-patient STRUCTURE, not longitudinal trajectories.\n\n")

# 1. does the module correlation structure match KNCC?
pairs <- t(combn(names(mods), 2))
cmp <- data.frame(m1 = pairs[, 1], m2 = pairs[, 2],
                  kncc_delta_rho = mapply(function(a, b) kncc_cor[a, b], pairs[, 1], pairs[, 2]),
                  cptac_rho = mapply(function(a, b) cptac_cor[a, b], pairs[, 1], pairs[, 2]))
ct <- cor.test(cmp$kncc_delta_rho, cmp$cptac_rho, method = "spearman")
cat("module-pair correlation structure: KNCC delta vs CPTAC cross-sectional\n")
cat("  rho =", round(ct$estimate, 3), " p =", signif(ct$p.value, 3),
    " (", nrow(cmp), " module pairs)\n\n")
print(cmp[order(-cmp$kncc_delta_rho), ], row.names = FALSE)

# 2. key axis checks
cat("\nkey axes:\n")
cat("  PROLIF vs NEURO: KNCC delta rho =", round(kncc_cor["PROLIF","NEURO"], 2),
    " | CPTAC rho =", round(cptac_cor["PROLIF","NEURO"], 2), "\n")
cat("  IL2STAT5 vs ROS: KNCC delta rho =", round(kncc_cor["IL2STAT5","ROS"], 2),
    " | CPTAC rho =", round(cptac_cor["IL2STAT5","ROS"], 2), "\n")
cat("  HYPOXIA vs EMT: KNCC delta rho =", round(kncc_cor["HYPOXIA","EMT"], 2),
    " | CPTAC rho =", round(cptac_cor["HYPOXIA","EMT"], 2), "\n")

# 3. inter-patient spread of class-B modules (CPTAC cross-sectional variance)
cat("\ninter-patient variance (CPTAC cross-sectional module scores):\n")
print(round(apply(S, 1, var, na.rm = TRUE), 4))
sink()

# figure: KNCC delta-structure vs CPTAC structure
p <- ggplot(cmp, aes(kncc_delta_rho, cptac_rho)) +
  geom_point(color = "#7B4EA3", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  theme_bw(base_size = 12) +
  labs(x = "KNCC longitudinal delta correlation (module pair)",
       y = "CPTAC cross-sectional correlation",
       title = sprintf("Module co-variation structure replicates (rho = %.2f, p = %s)",
                       ct$estimate, signif(ct$p.value, 2)))
ggsave("figures_kncc/cptac_matched_validation.pdf", p, width = 6, height = 5.5)
write.csv(cmp, "results/kncc/cptac_matched_validation.csv", row.names = FALSE)
cat("Done.\n")
