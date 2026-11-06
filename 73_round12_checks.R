# 73. Round-12 verification batch:
# (a) CGGA-693/325 patient overlap
# (b) phospho up/down breakdown
# (c) GENIE3 regulator-level paired test
# (d) protein WEBS z-score sensitivity
# (e) acetylome site-filter accounting + Jonckheere type
# (f) bootstrap pairing check (already paired by index reuse) - just document
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr) })

sink("results/round12_verifications.txt", split = TRUE)

# --- (a) CGGA overlap ---
c693 <- readRDS("data/raw/CGGA/cgga693_clinical.rds")
c325 <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
ov <- intersect(c693$sample, c325$sample)
cat("=== (a) CGGA-693 vs CGGA-325 ===\n")
cat("cgga693 samples:", length(c693$sample), " cgga325:", length(c325$sample),
    " overlap:", length(ov), "\n\n")

# --- (b) phospho up/down ---
ph <- read.csv("results/phospho_limma_LHS_vs_LLS.csv")
sig <- ph[ph$adj.P.Val < 0.05 & ph$is_spliceosome, ]
cat("=== (b) Spliceosome phosphosites (raw) ===\n")
cat("sig total:", nrow(sig), " up:", sum(sig$logFC > 0), " down:", sum(sig$logFC < 0), "\n")
phc <- read.csv("results/phosphoproteome_abundance_corrected.csv")
sigc <- phc[phc$adj.P.Val < 0.05 & phc$is_spliceosome, ]
cat("corrected sig total:", nrow(sigc), " up:", sum(sigc$logFC > 0), " down:", sum(sigc$logFC < 0), "\n")
srsf <- sig[grepl("^SRSF", sig$gene), ]
cat("SRSF sig sites:", nrow(srsf), " up:", sum(srsf$logFC > 0), "\n\n")

# --- (c) GENIE3 regulator-level ---
ko <- read.csv("results/virtual_ko_abs_summary.csv")
cat("=== (c) GENIE3 regulator-level test ===\n")
cat("14 regulator ratios (spliceosome/other):", paste(round(ko$ratio, 2), collapse = ", "), "\n")
wt <- wilcox.test(ko$ratio, mu = 1, alternative = "greater")
cat("One-sample Wilcoxon (ratio > 1, n =", nrow(ko), "regulators): p =", signif(wt$p.value, 3), "\n\n")

# --- (d) protein WEBS z-score sensitivity ---
prot <- read.delim("data/raw/CPTAC/HS_CPTAC_GBM_proteome_mssm_per_gene.cct", check.names = FALSE)
pm <- as.matrix(sapply(prot[, -1], as.numeric)); rownames(pm) <- prot[[1]]
webs <- read.csv("results/cptac_webs.csv")
writers <- c("HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT7","KAT8","CREBBP","EP300","KAT14")
erasers <- c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10",
             "HDAC11","SIRT1","SIRT2","SIRT3","SIRT5","SIRT6")
cs <- intersect(colnames(pm), webs$sample)
# z-score per protein across samples, then mean diff
z <- t(scale(t(pm[unique(c(writers, erasers)), cs])))   # proteins x samples, z per protein
w_z <- colMeans(z[writers, , drop = FALSE], na.rm = TRUE)
e_z <- colMeans(z[erasers, , drop = FALSE], na.rm = TRUE)
webs_z <- w_z - e_z
webs_orig <- webs$WEBS[match(cs, webs$sample)]
ct <- cor.test(webs_orig, webs_z, method = "spearman")
cat("=== (d) Protein WEBS z-score sensitivity ===\n")
cat("Samples with complete 10+16 proteins:", sum(complete.cases(t(z))), "/", length(cs), "\n")
cat("Spearman rho (original vs z-scored protein WEBS):", round(ct$estimate, 4),
    " p =", signif(ct$p.value, 3), "\n\n")

# --- (e) acetylome accounting ---
ac <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct", check.names = FALSE)
cat("=== (e) Acetylome site accounting ===\n")
cat("Total sites in file:", nrow(ac), "\n")
tt <- read.csv("results/acetyl_sites_limma.csv")
cat("Tested (limma):", nrow(tt), "\n")
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
tt$gene <- sub("-.*$", "", tt$site)
cat("spliceosome sites tested:", sum(tt$gene %in% splice), "\n")
jf <- read.csv("results/jonckheere_trend.csv")
cat("Jonckheere table rows:", nrow(jf), "\n")
pcol <- grep("^p|P\\.Value|pvalue", colnames(jf), value = TRUE, ignore.case = TRUE)[1]
if (!is.na(pcol)) cat("Jonckheere significant (p<0.05):", sum(jf[[pcol]] < 0.05), "/", nrow(jf), "\n")
sink()
cat(readLines("results/round12_verifications.txt"), sep = "\n")
