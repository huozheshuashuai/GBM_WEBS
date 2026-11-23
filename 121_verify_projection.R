# =============================================================================
# 121. Verify HHS frozen PC2 projection: raw contrast vs proper frozen
#      scaled PCA projection (center+scale learned in KNCC, applied to HHS)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr) })
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

# KNCC module-delta matrix and its PCA (correlation matrix, scaled)
Md <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
Md <- Md[, colSums(is.na(Md)) == 0]
pc <- prcomp(t(Md), scale. = TRUE)
pc2_load <- pc$rotation[, 2]
pc2_center <- pc$center    # frozen center (KNCC mean delta per module)
pc2_scale  <- pc$scale     # frozen scale (KNCC sd of delta per module)

mods_h <- lapply(names(mods_k), function(mn) intersect(mods_k[[mn]], rownames(hd)))
names(mods_h) <- names(mods_k)
mods_h <- mods_h[sapply(mods_h, length) >= 5]
mods_k <- mods_k[names(mods_h)]
Md_h <- t(sapply(mods_h, function(g) colMeans(hd[g, , drop = FALSE], na.rm = TRUE)))
Md_k <- t(sapply(mods_k, function(g) colMeans(kdelta[g, , drop = FALSE], na.rm = TRUE)))
load2 <- pc2_load[names(mods_h)]

sink("results/hhs/hhs_projection_verify.txt", split = TRUE)
cat("=== HHS frozen PC2 projection: mathematical variants ===\n\n")
cat("frozen PCA: correlation matrix (scale. = TRUE) on KNCC module-delta matrix\n")
cat("frozen parameters: center = KNCC per-module mean delta; scale = KNCC per-module sd\n\n")

# ---- Variant A (as submitted): loadings x raw HHS module deltas ----
A_k <- as.numeric(t(load2) %*% Md_k[names(mods_h), ])
A_h <- as.numeric(t(load2) %*% Md_h[names(mods_h), ])
cat("--- Variant A: frozen PC2-derived contrast (loadings x raw deltas, no centering) ---\n")
cat("    zero point = 'no change along the axis'\n")
cat("    KNCC: median", round(median(A_k), 3), " p =",
    signif(wilcox.test(A_k, mu = 0)$p.value, 3), "\n")
cat("    HHS:  median", round(median(A_h), 3), " p =",
    signif(wilcox.test(A_h, mu = 0)$p.value, 3), "\n")

# ---- Variant B: proper frozen scaled PCA projection (center+scale applied) ----
cen <- pc2_center[names(mods_h)]; scl <- pc2_scale[names(mods_h)]
Zh <- sweep(Md_h[names(mods_h), ], 1, cen, "-"); Zh <- sweep(Zh, 1, scl, "/")
Zk <- sweep(Md_k[names(mods_h), ], 1, cen, "-"); Zk <- sweep(Zk, 1, scl, "/")
B_k <- as.numeric(t(load2) %*% Zk)
B_h <- as.numeric(t(load2) %*% Zh)
cat("\n--- Variant B: strict frozen PCA projection (KNCC center+scale applied to HHS) ---\n")
cat("    zero point = 'KNCC patient-mean position along PC2'\n")
cat("    KNCC: median", round(median(B_k), 3), " p =",
    signif(wilcox.test(B_k, mu = 0)$p.value, 3), "\n")
cat("    HHS:  median", round(median(B_h), 3), " p =",
    signif(wilcox.test(B_h, mu = 0)$p.value, 3), "\n")
# bootstrap CI of median for B_h
bs <- replicate(10000, median(sample(B_h, length(B_h), TRUE)))
cat("    HHS median bootstrap 95% CI: [", round(quantile(bs, 0.025), 3), ",",
    round(quantile(bs, 0.975), 3), "]\n")

# ---- Variant B sensitivities ----
cat("\n--- Variant B sensitivities ---\n")
# (i) high-coverage modules
cov_frac <- sapply(mods_h, length) / sapply(mods_k, length)
hc <- names(cov_frac)[cov_frac >= 0.6]
Zh2 <- sweep(Md_h[hc, , drop = FALSE], 1, cen[hc], "-"); Zh2 <- sweep(Zh2, 1, scl[hc], "/")
B_h2 <- as.numeric(t(load2[hc]) %*% Zh2)
cat("    (i) >=60% coverage modules (n =", length(hc), "): p =",
    signif(wilcox.test(B_h2, mu = 0)$p.value, 3), "\n")
# (ii) exclude bottom-tercile abundance proteins (imputation suspects)
pm1 <- h$mat[, h$ptype == "Primary"]
abund <- rowMeans(pm1, na.rm = TRUE)
lowab <- names(abund)[abund <= quantile(abund, 1/3, na.rm = TRUE)]
mods_hc <- lapply(mods_h, function(g) setdiff(g, lowab))
mods_hc <- mods_hc[sapply(mods_hc, length) >= 5]
Md_h3 <- t(sapply(mods_hc, function(g) colMeans(hd[g, , drop = FALSE], na.rm = TRUE)))
Zh3 <- sweep(Md_h3[names(mods_hc), , drop = FALSE], 1, cen[names(mods_hc)], "-")
Zh3 <- sweep(Zh3, 1, scl[names(mods_hc)], "/")
B_h3 <- as.numeric(t(load2[names(mods_hc)]) %*% Zh3)
cat("    (ii) excluding imputation-suspect proteins: p =",
    signif(wilcox.test(B_h3, mu = 0)$p.value, 3), "\n")
# (iii) leave-one-patient-out
loo <- sapply(seq_len(ncol(Zh)), function(j) wilcox.test(as.numeric(t(load2) %*% Zh[, -j]), mu = 0)$p.value)
cat("    (iii) LOPO p range: [", signif(min(loo), 2), ",", signif(max(loo), 2),
    "] | all < 0.05:", all(loo < 0.05), "\n")

# ---- direction check ----
cat("\n--- direction bookkeeping ---\n")
cat("PC2 loadings: inflammatory modules negative (IL2STAT5 -0.318, XENO -0.313,\n")
cat("HYPOXIA -0.294); proliferative modules positive (G2M +0.201, E2F +0.175).\n")
cat("Module deltas in HHS are predominantly positive for inflammatory modules\n")
cat("-> positive x negative = negative projection. A MORE NEGATIVE PC2 score\n")
cat("indicates STRONGER inflammatory-metabolic change at recurrence.\n")
sink()

saveRDS(list(variant_A = list(kncc = A_k, hhs = A_h),
             variant_B = list(kncc = B_k, hhs = B_h,
                              hhs_hc = B_h2, hhs_noimpute = B_h3, loo_p = loo),
             loadings = load2, center = cen, scale = scl),
        "results/hhs/projection_variants.rds")
cat("Done.\n")
