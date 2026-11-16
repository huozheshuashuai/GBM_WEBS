# =============================================================================
# 119. Archive HHS effect sizes with bootstrap CIs (pre-submission archive)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(dplyr))
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

boot_med <- function(x, B = 10000) {
  bs <- replicate(B, median(sample(x, length(x), TRUE)))
  c(mean = mean(x), median = median(x),
    ci_lo = unname(quantile(bs, 0.025)), ci_hi = unname(quantile(bs, 0.975)),
    wilcox_p = suppressWarnings(wilcox.test(x, mu = 0)$p.value))
}
boot_mean <- boot_med
ci <- rbind(KNCC = boot_mean(pc2_k), HHS = boot_mean(pc2_h))
out <- data.frame(cohort = rownames(ci), n = c(length(pc2_k), length(pc2_h)), round(ci[, c("median","ci_lo","ci_hi","wilcox_p")], 5))
write.csv(out, "results/hhs/hhs_pc2_effect_sizes.csv", row.names = FALSE)
print(out)

# module-level archive with CIs
mod_out <- lapply(names(mods_h), function(mn) {
  ck <- boot_mean(Md_k[mn, ]); ch <- boot_mean(Md_h[mn, ])
  data.frame(module = mn, n_shared = length(mods_h[[mn]]),
             kncc_med = ck["median"], kncc_ci_lo = ck["ci_lo"], kncc_ci_hi = ck["ci_hi"],
             hhs_med = ch["median"], hhs_ci_lo = ch["ci_lo"], hhs_ci_hi = ch["ci_hi"],
             hhs_p = ch["wilcox_p"])
})
mod_out <- do.call(rbind, mod_out)
mod_out$hhs_q <- p.adjust(mod_out$hhs_p, "BH")
write.csv(mod_out, "results/hhs/hhs_module_effect_sizes.csv", row.names = FALSE)
cat("module archive written:", nrow(mod_out), "rows\n")
