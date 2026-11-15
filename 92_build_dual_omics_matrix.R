# =============================================================================
# 92. KNCC unified dual-omics T1/T2 analysis matrix
#     One bundle: proteome (per-sample + paired deltas) + phosphoproteome
#     (per-sample + paired deltas) + master patient table + clinical
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr) })
set.seed(42)

# ---------- proteome (from script 86 outputs) ----------
po <- readRDS("results/kncc/proteome_analysis_ready.rds")
pm <- po$pm; cli <- po$annot
prot_paired <- po$paired; prot_delta <- po$delta

# ---------- phosphoproteome: rebuild per-sample matrix (same rules as 90) ----------
f <- "E:/KIMI/data/PDC000515/KNCC_Glioblastoma_Evolution_Phosphoproteome.phosphosite.tmt11.tsv"
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
log_idx <- which(grepl("^Log KNCC_GBM", hdr))
gene_col <- which(hdr == "Gene")
cc <- rep("NULL", length(hdr)); cc[1] <- "character"; cc[gene_col] <- "character"; cc[log_idx] <- "numeric"
ph <- read.delim(f, sep = "\t", header = TRUE, colClasses = cc, check.names = FALSE)
colnames(ph)[-1] <- sub("\\.[0-9]+$", "", sub(":.*", "", sub("^Log ", "", colnames(ph)[-1])))
num_cols <- which(grepl("^KNCC_GBM", colnames(ph)))
m <- as.matrix(ph[, num_cols]); cn <- colnames(m)
mu <- matrix(NA_real_, nrow(m), length(unique(cn)), dimnames = list(NULL, unique(cn)))
for (s in unique(cn)) { idx <- which(cn == s); mu[, s] <- if (length(idx) == 1) m[, idx] else rowMeans(m[, idx, drop = FALSE], na.rm = TRUE) }
phm <- mu; rownames(phm) <- ph[[1]]
ph_gene <- ph$Gene

# phospho annotation + pairing
samples <- colnames(phm)
ph_annot <- data.frame(sample = samples,
                       patient = sub("_(T[12])$", "", samples),
                       timepoint = sub(".*_", "", samples))
ptab <- table(ph_annot$patient, ph_annot$timepoint)
ph_paired <- rownames(ptab)[ptab[, "T1"] > 0 & ptab[, "T2"] > 0]
t1 <- ph_annot$sample[ph_annot$timepoint == "T1"][match(ph_paired, ph_annot$patient[ph_annot$timepoint == "T1"])]
t2 <- ph_annot$sample[ph_annot$timepoint == "T2"][match(ph_paired, ph_annot$patient[ph_annot$timepoint == "T2"])]
ph_delta <- phm[, t2] - phm[, t1]
colnames(ph_delta) <- ph_paired

# ---------- master patient table ----------
cli_tab <- read.csv("E:/KIMI/data/PDC000514/KNCC_clinical_tableS1.csv")
cli_tab <- cli_tab[!duplicated(cli_tab$PDC_ID) & grepl("^GBM", cli_tab$PDC_ID), ]
all_pts <- sort(unique(c(cli$patient, ph_annot$patient)))
master <- data.frame(
  patient = all_pts,
  PDC_ID = sub("KNCC_", "", all_pts),
  prot_T1 = all_pts %in% cli$patient[cli$timepoint == "T1"],
  prot_T2 = all_pts %in% cli$patient[cli$timepoint == "T2"],
  phos_T1 = all_pts %in% ph_annot$patient[ph_annot$timepoint == "T1"],
  phos_T2 = all_pts %in% ph_annot$patient[ph_annot$timepoint == "T2"]
)
master$prot_paired <- master$prot_T1 & master$prot_T2
master$phos_paired <- master$phos_T1 & master$phos_T2
master$dual_paired <- master$prot_paired & master$phos_paired
master <- master %>% left_join(cli_tab %>% select(PDC_ID, Institution, SEX, AGE, PFS, OS,
                                                  `death.0.1.`), by = "PDC_ID")
write.csv(master, "results/kncc/master_patient_table.csv", row.names = FALSE)

sink("results/kncc/dual_omics_matrix_report.txt", split = TRUE)
cat("=== KNCC unified dual-omics matrix ===\n\n")
cat("proteome:   ", nrow(pm), "proteins x", ncol(pm), "samples |",
    sum(master$prot_paired), "paired patients\n")
cat("phospho:    ", nrow(phm), "sites x", ncol(phm), "samples |",
    sum(master$phos_paired), "paired patients\n")
cat("dual paired (all four layers T1/T2 x prot/phos):", sum(master$dual_paired), "\n")
cat("patients total:", nrow(master), "\n\n")
cat("clinical coverage in master: SEX", sum(!is.na(master$SEX)), "| PFS",
    sum(!is.na(master$PFS)), "| OS", sum(!is.na(master$OS)), "\n")
cat("\nproteome delta missingness (median per patient):",
    round(median(apply(prot_delta, 2, function(x) mean(is.na(x)))), 3), "\n")
cat("phospho delta missingness (median per patient):",
    round(median(apply(ph_delta, 2, function(x) mean(is.na(x)))), 3), "\n")
sink()

# save unified bundle
saveRDS(list(
  # per-sample matrices (log2 ratio vs pooled reference)
  prot_samples = pm, prot_annot = cli,
  phos_samples = phm, phos_annot = ph_annot, phos_gene = ph_gene,
  # paired patient deltas
  prot_delta = prot_delta, phos_delta = ph_delta,
  prot_paired = prot_paired, phos_paired = ph_paired,
  dual_paired = master$patient[master$dual_paired],
  master = master
), "results/kncc/dual_omics_bundle.rds")

cat("\nbundle saved. Dual-paired patients:", sum(master$dual_paired), "\n")
cat("Done.\n")
