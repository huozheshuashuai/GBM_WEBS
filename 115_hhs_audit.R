# =============================================================================
# 115. HHS (Tatari 2022) feasibility audit + analysis-ready matrix build
#      Table 2 = imputed log2 LFQ (6,977 proteins x 134 samples)
#      Table 1 = patient metadata (GBMxx-PM/RM matched pairs)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(readxl); library(dplyr) })
dir.create("results/hhs", showWarnings = FALSE)

f2 <- "E:/KIMI/data/HHS_Tatari2022/NIHMS1894660-supplement-Supplementary_Table_2.xlsx"
raw <- read_excel(f2, sheet = "Imputed values", col_names = FALSE, skip = 3)
colnames(raw) <- as.character(raw[1, ])
raw <- raw[-1, ]
cat("matrix:", nrow(raw), "proteins x", ncol(raw) - 3, "samples\n")
cat("first cols:", paste(colnames(raw)[1:8], collapse = " | "), "\n")

genes <- raw$Gene.names
acc <- raw$Accession
mat <- as.matrix(sapply(raw[, -(1:3)], as.numeric))
rownames(mat) <- make.unique(ifelse(is.na(genes) | genes == "", acc, genes))
cat("samples:", paste(head(colnames(mat), 10), collapse = ", "), "...\n")

# sample types from names
sn <- colnames(mat)
ptype <- ifelse(grepl("-N$", sn), "NAT",
         ifelse(grepl("-PM$", sn), "Primary",
         ifelse(grepl("-RM$", sn), "Recurrent",
         ifelse(grepl("-PNM$", sn), "PrimaryNonMatched", "RecurrentNonMatched"))))
cat("\nsample types:", paste(table(ptype), collapse = ", "), "\n")

# patient pairing: GBMxx-PM / GBMxx-RM
pid <- sub("-(PM|RM|PNM|RNM|N)$", "", sn)
tab <- table(pid[ptype %in% c("Primary", "Recurrent")], ptype[ptype %in% c("Primary", "Recurrent")])
paired <- rownames(tab)[tab[, "Primary"] > 0 & tab[, "Recurrent"] > 0]
cat("\npatients with PM+RM proteome:", length(paired), "\n")
cat("patients only PM:", sum(tab[, "Primary"] > 0 & tab[, "Recurrent"] == 0),
    " only RM:", sum(tab[, "Primary"] == 0 & tab[, "Recurrent"] > 0), "\n")

# multi-specimen rule: does any patient have >1 PM or >1 RM?
dup <- names(which(table(c(pid)) > 2))
cat("patients with >2 specimens:", paste(dup, collapse = ", "), "\n")

# metadata from Table 1
t1 <- read_excel("E:/KIMI/data/HHS_Tatari2022/NIHMS1894660-supplement-Supplementary_Table_1.xlsx",
                 col_names = FALSE, skip = 1)
meta <- data.frame(id = t1[[1]], sex = t1[[2]], age = t1[[3]],
                   ttr_months = t1[[4]], ttd_months = t1[[5]], treatment = t1[[6]])
meta <- meta[grepl("^GBM", meta$id), ]
meta$patient <- sub("-(PM|RM|PNM|RNM|N)$", "", meta$id)
meta$type <- sub(".*-", "", meta$id)
cat("\nTable 1 metadata rows:", nrow(meta), "| patients:", length(unique(meta$patient)), "\n")

# module coverage: frozen 7 + PC2 modules
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(mat)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(mat)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(mat)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)
cat("\nmodule coverage in HHS matrix:\n")
for (mn in names(mods)) cat(" ", mn, ":", length(mods[[mn]]), "\n")

# compartment proxies coverage
proxies <- list(
  Microglia = c("PTPRC","CD14","CSF1R","C1QA","AIF1","TMEM119","P2RY12","CX3CR1","CD68","CD163"),
  MES_state = c("CD44","VIM","S100A11","S100A10","ANXA1","LGALS1","MT2A","C1S","NAMPT","EMP3"),
  Endothelial = c("PECAM1","VWF","CLDN5","CDH5","KDR","ESAM","FLT1","ENG"),
  Neuron = c("RBFOX3","SNAP25","SYP","SYT1","STMN2","TUBB3","MAP2","NEFL")
)
for (pn in names(proxies))
  cat(" proxy", pn, ":", sum(proxies[[pn]] %in% rownames(mat)), "/", length(proxies[[pn]]), "\n")

# build delta matrix for paired patients
t1c <- mat[, ptype == "Primary"]; t2c <- mat[, ptype == "Recurrent"]
colnames(t1c) <- sub("-PM$", "", colnames(t1c))
colnames(t2c) <- sub("-RM$", "", colnames(t2c))
common_pt <- intersect(colnames(t1c), colnames(t2c))
delta <- t2c[, common_pt] - t1c[, common_pt]
cat("\ndelta matrix:", nrow(delta), "proteins x", ncol(delta), "paired patients\n")

saveRDS(list(mat = mat, delta = delta, paired = common_pt, meta = meta,
             ptype = ptype, pid = pid), "results/hhs/hhs_analysis_ready.rds")
write.csv(data.frame(patient = common_pt), "results/hhs/hhs_paired_patients.csv", row.names = FALSE)

sink("results/hhs/hhs_audit.txt", split = TRUE)
cat("=== HHS (Tatari 2022) feasibility audit ===\n\n")
cat("matrix: 6,977 proteins x 134 samples (imputed log2 LFQ-iBAQ)\n")
cat("value provenance: Gaussian downshift imputation (width 0.3, downshift 1.8)\n")
cat("CAUTION: imputation inflates apparent coverage; low-abundance module scores\n")
cat("may be compressed toward zero delta. Sensitivity: restrict to proteins\n")
cat("quantified in original (unimputed) matrix if available; otherwise treat\n")
cat("effect sizes as conservative.\n\n")
cat("sample types:", paste(names(table(ptype)), table(ptype), collapse = ", "), "\n")
cat("paired patients (PM+RM):", length(paired), "\n")
cat("module coverage:", paste(names(mods), sapply(mods, length), collapse = ", "), "\n")
sink()
cat("Done.\n")
