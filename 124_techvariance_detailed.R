# =============================================================================
# 124. Technical-variance detailed audit (advisor: per-duplicate reporting)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(dplyr))
set.seed(42)

f <- "E:/KIMI/data/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log KNCC_GBM", hdr))
cc <- rep("NULL", length(hdr)); cc[1] <- "character"; cc[log_idx] <- "numeric"
pr <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc, check.names = FALSE)
pr <- pr[!pr$Gene %in% c("Mean","Median","StdDev") & !is.na(pr$Gene), ]
raw_mat <- as.matrix(pr[, -1]); rownames(raw_mat) <- pr$Gene
colnames(raw_mat) <- sub(":.*", "", sub("^Log ", "", colnames(raw_mat)))
dup_names <- unique(colnames(raw_mat)[duplicated(colnames(raw_mat))])

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

b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta

sink("results/kncc/techvariance_detailed.txt", split = TRUE)
cat("=== per-duplicate technical differences (module-score level) ===\n\n")
dup_report <- lapply(dup_names, function(s) {
  idx <- which(colnames(raw_mat) == s)
  m1 <- raw_mat[, idx[1]]; m2 <- raw_mat[, idx[2]]
  d <- sapply(mods, function(g) mean(m1[g] - m2[g], na.rm = TRUE))
  c(sample = s, sd_across_modules = sd(d), max_abs = max(abs(d)),
    rmse = sqrt(mean(d^2)))
})
dup_df <- do.call(rbind, dup_report) %>% as.data.frame()
print(dup_df, row.names = FALSE)

cat("\n=== per-module: technical variance vs inter-patient variance (with CIs) ===\n\n")
cat("estimator: var(dup1-dup2 of module score) = technical variance of a two-\n")
cat("measurement delta (incl. plex effects -> conservative upper bound)\n")
cat("NOTE: duplicates measure ANALYTICAL variance of the same tissue lysate;\n")
cat("they do NOT capture regional sampling differences between T1 and T2.\n\n")

res <- lapply(names(mods), function(mn) {
  g <- mods[[mn]]
  tech_diffs <- sapply(dup_names, function(s) {
    idx <- which(colnames(raw_mat) == s)
    mean(raw_mat[g, idx[1]] - raw_mat[g, idx[2]], na.rm = TRUE)
  })
  tv <- var(tech_diffs)
  pd <- colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
  pv <- var(pd, na.rm = TRUE)
  # bootstrap CI of ratio over duplicate pairs
  bs <- replicate(5000, pv / var(sample(tech_diffs, length(tech_diffs), TRUE)))
  data.frame(module = mn, n_genes = length(g),
             tech_var = round(tv, 5), patient_var = round(pv, 4),
             ratio = round(pv / tv, 1),
             ratio_ci_lo = round(quantile(bs, 0.025, na.rm = TRUE), 1),
             ratio_ci_hi = round(quantile(bs, 0.975, na.rm = TRUE), 1))
})
res <- do.call(rbind, res)
print(res, row.names = FALSE)
write.csv(res, "results/kncc/techvariance_detailed.csv", row.names = FALSE)
sink()
cat("Done.\n")
