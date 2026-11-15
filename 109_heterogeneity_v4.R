# =============================================================================
# 109. Task 4+5 (V4.0): proteome QC audit + formal heterogeneity test
#      with full sensitivity battery for the 7 frozen modules
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot; delta <- b$prot_delta; paired <- b$prot_paired
dd <- read.csv("02_PATIENT_DATA_DICTIONARY_V2.csv")

# ---------------- Task 4: QC audit ----------------
sink("results/kncc/qc_audit_v4.txt", split = TRUE)
cat("=== Task 4: proteome QC audit ===\n\n")

# outlier samples: median correlation with all other samples
cmat <- cor(pm, use = "pairwise.complete.obs")
med_cor <- apply(cmat, 1, median, na.rm = TRUE)
outliers <- names(med_cor)[med_cor < quantile(med_cor, 0.01)]
cat("sample median correlation: min", round(min(med_cor), 3),
    "| 1st pct:", round(quantile(med_cor, 0.01), 3), "\n")
cat("outlier samples (<1st pct):", length(outliers), "\n")
if (length(outliers)) {
  odf <- cli[cli$sample %in% outliers, c("sample", "patient", "timepoint", "plex")]
  print(odf)
}

# missingness vs institution/plex/timepoint
cli$miss <- colMeans(is.na(pm))[cli$sample]
cli$inst <- dd$Institution[match(cli$patient, paste0("KNCC_", dd$PDC_ID))]
cat("\nmissingness by timepoint: T1", round(median(cli$miss[cli$timepoint == "T1"]), 3),
    " T2", round(median(cli$miss[cli$timepoint == "T2"]), 3),
    " Wilcoxon p =", signif(wilcox.test(miss ~ timepoint, data = cli)$p.value, 3), "\n")
kw_inst <- kruskal.test(miss ~ inst, data = cli)
cat("missingness by institution: KW p =", signif(kw_inst$p.value, 3), "\n")
kw_plex <- kruskal.test(miss ~ plex, data = cli)
cat("missingness by plex: KW p =", signif(kw_plex$p.value, 3), "\n")

# institution x timepoint structure
cat("\npatients per institution (paired):\n")
print(table(dd$Institution[dd$prot_paired]))
sink()

# ---------------- Task 5: formal heterogeneity ----------------
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(delta)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(delta)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)

# technical variance from duplicate aliquots (recompute from raw file)
f <- "E:/KIMI/data/PDC000514/KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv"
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log KNCC_GBM", hdr))
cc_cols <- rep("NULL", length(hdr)); cc_cols[1] <- "character"; cc_cols[log_idx] <- "numeric"
pr <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc_cols, check.names = FALSE)
pr <- pr[!pr$Gene %in% c("Mean","Median","StdDev") & !is.na(pr$Gene), ]
cn_base <- sub(":.*", "", sub("^Log ", "", colnames(pr)[-1]))
dup_names <- unique(cn_base[duplicated(cn_base)])
tech_var <- mean(sapply(dup_names, function(s) {
  m <- as.matrix(pr[, -1][, which(cn_base == s), drop = FALSE]); mode(m) <- "numeric"
  var(m[, 1] - m[, 2], na.rm = TRUE)
}))
cat("technical variance per protein:", round(tech_var, 4), "\n")

# main analysis function
analyze_set <- function(dmat, label) {
  t(sapply(names(mods), function(mn) {
    g <- mods[[mn]]
    d <- colMeans(dmat[g, , drop = FALSE], na.rm = TRUE)
    d <- d[!is.na(d)]
    n <- length(d)
    # patient bootstrap CI for SD
    bs <- replicate(2000, sd(sample(d, n, replace = TRUE)))
    wt <- suppressWarnings(wilcox.test(d, mu = 0))
    c(mean_delta = mean(d), sd_delta = sd(d),
      sd_ci_lo = quantile(bs, 0.025), sd_ci_hi = quantile(bs, 0.975),
      dir_consistency = max(mean(d > 0), mean(d < 0)),
      paired_p = wt$p.value, n = n,
      tech_ratio = var(d) / (tech_var / length(g)))
  }))
}

# (a) all 105 pairs
res_all <- analyze_set(delta, "all")
# (b) same-plex pairs only
same_plex_pts <- dd$PDC_ID[dd$same_plex_pair %in% TRUE]
same_plex_pts <- paste0("KNCC_", same_plex_pts)
d_spx <- delta[, colnames(delta) %in% same_plex_pts]
res_spx <- analyze_set(d_spx, "same-plex")
# (c) exclude pediatric
ped <- paste0("KNCC_", dd$PDC_ID[dd$pediatric_flag %in% TRUE])
d_nopd <- delta[, !colnames(delta) %in% ped]
res_nopd <- analyze_set(d_nopd, "no-pediatric")
# (d) leave-one-institution-out
insts <- unique(na.omit(dd$Institution))
res_loio <- lapply(insts, function(inst) {
  excl <- paste0("KNCC_", dd$PDC_ID[dd$Institution == inst])
  analyze_set(delta[, !colnames(delta) %in% excl], paste0("no-", inst))
})
# (e) missingness threshold variant: proteins quantified in >=80% of samples
keep80 <- rowSums(!is.na(pm)) >= ncol(pm) * 0.8
delta80 <- delta[keep80, ]
mods80 <- lapply(mods, function(g) intersect(g, rownames(delta80)))
mods80 <- mods80[sapply(mods80, length) >= 8]
mods_orig <- mods; mods <- mods80
res_m80 <- analyze_set(delta80, "miss80")
mods <- mods_orig

sink("results/kncc/heterogeneity_v4.txt", split = TRUE)
cat("=== Task 5: formal heterogeneity (7 frozen modules) ===\n\n")
print_one <- function(res, label) {
  cat("\n---", label, "---\n")
  df <- data.frame(module = rownames(res), round(res, 4))
  print(df, row.names = FALSE)
}
print_one(res_all, "(a) all 105 pairs")
print_one(res_spx, "(b) same-plex pairs (n=100)")
print_one(res_nopd, "(c) excluding 2 pediatric patients")
print_one(res_m80, "(e) proteins >=80% quantified")
cat("\n--- (d) leave-one-institution-out: sd_delta per module ---\n")
for (i in seq_along(insts)) {
  cat("excluding", as.character(insts[i]), ":\n")
  print(round(res_loio[[i]][, c("sd_delta", "tech_ratio")], 3))
}
sink()

# consolidate
write.csv(data.frame(module = rownames(res_all), res_all, row.names = NULL),
          "results/kncc/heterogeneity_all.csv", row.names = FALSE)

# figure: patient delta distributions per module with bootstrap SD CIs
plot_df <- do.call(rbind, lapply(names(mods), function(mn) {
  d <- colMeans(delta[mods[[mn]], , drop = FALSE], na.rm = TRUE)
  data.frame(module = mn, delta = d)
}))
sdd <- data.frame(module = rownames(res_all), sd = res_all[, "sd_delta"],
                  lo = res_all[, "sd_ci_lo.2.5%"], hi = res_all[, "sd_ci_hi.97.5%"])
p1 <- ggplot(plot_df, aes(module, delta, fill = module)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  coord_flip() + labs(y = "Per-patient module delta (T2 - T1)", x = NULL,
                      title = "Inter-patient heterogeneity of recurrence evolution")
sdd$module <- factor(sdd$module, levels = rownames(res_all))
p2 <- ggplot(sdd, aes(module, sd, ymin = lo, ymax = hi)) +
  geom_pointrange(color = "#7B4EA3") + coord_flip() +
  theme_bw(base_size = 11) +
  labs(y = "Inter-patient SD (bootstrap 95% CI)", x = NULL,
       title = "Heterogeneity magnitude with bootstrap CI")
library(patchwork)
ggsave("figures_kncc/heterogeneity_v4.pdf", p1 + p2, width = 12, height = 5)
cat("Done.\n")
