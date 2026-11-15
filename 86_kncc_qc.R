# =============================================================================
# 86. KNCC longitudinal proteome (PDC000514): QC + pairing + analysis-ready data
#     Module 1 of the KNCC longitudinal plan
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

dir.create("results/kncc", showWarnings = FALSE)
dir.create("figures_kncc", showWarnings = FALSE)

base <- "E:/KIMI/data/PDC000514/"
f <- paste0(base, "KNCC_Glioblastoma_Evolution_Proteome.tmt11.tsv")

# ---- read header, select Gene + Log columns (drop Unshared for main analysis) ----
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log ", hdr))
cc <- rep("NULL", length(hdr)); cc[1] <- "character"; cc[log_idx] <- "numeric"
prot <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc, check.names = FALSE)
cat("rows read:", nrow(prot), "\n")

# drop summary rows (Mean/Median/StdDev) and empty genes
prot <- prot[!prot$Gene %in% c("Mean", "Median", "StdDev") & !is.na(prot$Gene) & prot$Gene != "", ]
cat("protein rows after cleaning:", nrow(prot), "\n")

colnames(prot)[-1] <- sub(":.*", "", sub("^Log ", "", colnames(prot)[-1]))
pm <- as.matrix(prot[, -1]); rownames(pm) <- prot$Gene

# keep only KNCC GBM tumor aliquots; drop "Other" (non-GBM/QC) and "ref" channels
tumor_cols <- grepl("^KNCC_GBM", colnames(pm))
cat("dropping", sum(!tumor_cols), "non-tumor columns (Other/ref)\n")
pm <- pm[, tumor_cols]

# strip read.delim's ".1" de-duplication suffixes, then merge true duplicate aliquots
colnames(pm) <- sub("\\.[0-9]+$", "", colnames(pm))

# ---- sample annotation ----
samples <- colnames(pm)
dup_s <- samples[duplicated(samples)]
if (length(dup_s) > 0) {
  cat("duplicate sample columns:", length(dup_s), "->", paste(head(unique(dup_s), 5), collapse = ", "), "\n")
  # average duplicate aliquots (same sample assayed in >1 plex)
  pm_u <- matrix(NA_real_, nrow(pm), length(unique(samples)),
                 dimnames = list(rownames(pm), unique(samples)))
  for (s in unique(samples)) {
    idx <- which(samples == s)
    pm_u[, s] <- if (length(idx) == 1) pm[, idx] else rowMeans(pm[, idx, drop = FALSE], na.rm = TRUE)
  }
  pm <- pm_u
  samples <- colnames(pm)
  cat("matrix after merging duplicates:", ncol(pm), "samples\n")
}
patient <- sub("_(T[12])$", "", samples)          # KNCC_GBM0001
timepoint <- sub(".*_", "", samples)              # T1 / T2
annot <- data.frame(sample = samples, patient = patient, timepoint = timepoint)

# TMT plex from sample manifest (a sample may appear in >1 plex; keep all)
sman <- read.delim(paste0(base, "KNCC_Glioblastoma_Evolution_Proteome.sample.txt"),
                   check.names = FALSE)
sman <- sman[!duplicated(sman$AnalyticalSample), ]
plex_id <- sub(".*(S\\d+):.*", "\\1", sman$AnalyticalSample)
plex_map <- list()
for (ch in c("126C","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C")) {
  vals <- sman[[ch]]
  aliq <- sub(":.*", "", vals)
  for (i in seq_along(aliq)) if (aliq[i] != "" && !grepl("^ref", aliq[i]))
    plex_map[[aliq[i]]] <- unique(c(plex_map[[aliq[i]]], plex_id[i]))
}
annot$plex <- sapply(annot$sample, function(s)
  if (s %in% names(plex_map)) paste(plex_map[[s]], collapse = ";") else NA_character_)
cat("samples without plex mapping:", sum(is.na(annot$plex)),
    if (any(is.na(annot$plex))) paste("->", paste(annot$sample[is.na(annot$plex)], collapse = ", ")) else "", "\n")

# clinical join (Table S1: PDC_ID = GBM0001); dedupe clinical rows by PDC_ID
cli_tab <- read.csv(paste0(base, "KNCC_clinical_tableS1.csv"))
cli_tab <- cli_tab[!duplicated(cli_tab$PDC_ID) & grepl("^GBM", cli_tab$PDC_ID), ]
annot$PDC_ID <- sub("KNCC_", "", annot$patient)
cli <- annot %>% left_join(cli_tab, by = "PDC_ID")
cat("samples with clinical match:", sum(!is.na(cli$SEX)), "/", nrow(cli), "\n")

write.csv(cli, "results/kncc/sample_annotation.csv", row.names = FALSE)

# ---- pairing ----
ptab <- table(cli$patient, cli$timepoint)
paired_pts <- rownames(ptab)[ptab[, "T1"] > 0 & ptab[, "T2"] > 0]
cat("patients:", length(unique(cli$patient)),
    "| T1 samples:", sum(cli$timepoint == "T1"),
    "| T2 samples:", sum(cli$timepoint == "T2"),
    "| paired:", length(paired_pts), "\n")

# ---- QC metrics ----
miss_prot <- rowMeans(is.na(pm))
miss_samp <- colMeans(is.na(pm))
cat("protein missing rate: median", round(median(miss_prot), 3),
    "| >50% missing:", sum(miss_prot > 0.5), "proteins\n")
cat("sample missingness: median", round(median(miss_samp), 3),
    "range", round(range(miss_samp), 3), "\n")

# sample-sample correlation (complete pairwise)
set.seed(1)
cmat <- cor(pm, use = "pairwise.complete.obs")
diag(cmat) <- NA
within <- sapply(paired_pts, function(p) {
  s1 <- cli$sample[cli$patient == p & cli$timepoint == "T1"]
  s2 <- cli$sample[cli$patient == p & cli$timepoint == "T2"]
  if (length(s1) == 1 && length(s2) == 1) cmat[s1, s2] else NA
})
same_p <- outer(cli$patient, cli$patient, "==")
off_between <- cmat[upper.tri(cmat) & !same_p]
cat("within-patient T1-T2 correlation: median", round(median(within, na.rm = TRUE), 3), "\n")
cat("between-patient correlation: median", round(median(off_between, na.rm = TRUE), 3), "\n")

# ---- PCA (on complete-ish proteins) ----
keep <- rowSums(!is.na(pm)) >= ncol(pm) * 0.8
pmq <- pm[keep, ]
cat("proteins with >=80% quantification:", nrow(pmq), "\n")
pm_imp <- pmq
for (j in seq_len(ncol(pm_imp))) {
  nas <- is.na(pm_imp[, j])
  pm_imp[nas, j] <- rowMeans(pm_imp[nas, , drop = FALSE], na.rm = TRUE)
}
pca <- prcomp(t(pm_imp), scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pcdf <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   timepoint = cli$timepoint, plex = sub(";.*", "", cli$plex),
                   patient = cli$patient)

p1 <- ggplot(pcdf, aes(PC1, PC2, color = timepoint)) +
  geom_point(alpha = 0.8) + theme_bw(base_size = 12) +
  labs(title = "KNCC proteome PCA by timepoint",
       x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"))
ggsave("figures_kncc/qc_pca_timepoint.pdf", p1, width = 6, height = 5)

p2 <- ggplot(pcdf, aes(PC1, PC2, color = plex)) +
  geom_point(alpha = 0.8) + theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(title = "KNCC proteome PCA by TMT plex",
       x = paste0("PC1 (", round(var_exp[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2] * 100, 1), "%)"))
ggsave("figures_kncc/qc_pca_plex.pdf", p2, width = 6, height = 5)

# plex effect: ANOVA of PC1/PC2 by plex and timepoint
a_plex1 <- summary(aov(PC1 ~ plex, data = pcdf))[[1]][1, "Pr(>F)"]
a_plex2 <- summary(aov(PC2 ~ plex, data = pcdf))[[1]][1, "Pr(>F)"]
a_tp1 <- summary(aov(PC1 ~ timepoint, data = pcdf))[[1]][1, "Pr(>F)"]
a_tp2 <- summary(aov(PC2 ~ timepoint, data = pcdf))[[1]][1, "Pr(>F)"]

# ---- paired delta distribution (T2 - T1), all proteins ----
s1 <- cli$sample[match(paired_pts, cli$patient)]  # need T1 per patient
t1_map <- cli$sample[cli$timepoint == "T1"][match(paired_pts, cli$patient[cli$timepoint == "T1"])]
t2_map <- cli$sample[cli$timepoint == "T2"][match(paired_pts, cli$patient[cli$timepoint == "T2"])]
stopifnot(all(!is.na(t1_map)), all(!is.na(t2_map)))
delta <- pm[, t2_map] - pm[, t1_map]
colnames(delta) <- paired_pts
cat("delta matrix:", nrow(delta), "proteins x", ncol(delta), "paired patients\n")
cat("median protein delta (T2-T1):", round(median(apply(delta, 1, median, na.rm = TRUE), na.rm = TRUE), 4), "\n")

# ---- save analysis-ready objects ----
saveRDS(list(pm = pm, annot = cli, delta = delta, paired = paired_pts,
             pca_var = var_exp),
        "results/kncc/proteome_analysis_ready.rds")

sink("results/kncc/qc_report.txt", split = TRUE)
cat("=== KNCC PDC000514 QC report ===\n\n")
cat("matrix:", nrow(pm), "proteins x", ncol(pm), "samples (Log version)\n")
cat("patients:", length(unique(cli$patient)), "| paired T1+T2:", length(paired_pts), "\n")
cat("clinical matched:", sum(!is.na(cli$SEX)), "of", nrow(cli), "samples\n")
cat("plexes:", length(unique(cli$plex)), "\n")
cat("median sample missingness:", round(median(miss_samp), 3), "\n")
cat("proteins >=80% quantified:", nrow(pmq), "\n")
cat("within-patient T1-T2 cor: median", round(median(within, na.rm = TRUE), 3), "\n")
cat("between-patient cor: median", round(median(off_between, na.rm = TRUE), 3), "\n")
cat("PC1 variance:", round(var_exp[1] * 100, 1), "% | PC2:", round(var_exp[2] * 100, 1), "%\n")
cat("PC1 ~ plex p =", signif(a_plex1, 3), "| PC2 ~ plex p =", signif(a_plex2, 3), "\n")
cat("PC1 ~ timepoint p =", signif(a_tp1, 3), "| PC2 ~ timepoint p =", signif(a_tp2, 3), "\n")
cat("median T2-T1 delta across proteins:", round(median(apply(delta, 1, median, na.rm = TRUE), na.rm = TRUE), 4), "\n")
sink()
cat("Done.\n")
