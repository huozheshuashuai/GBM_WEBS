# =============================================================================
# 120. Pre-submission consistency audit: manuscript numbers vs result files
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
ms <- readLines("MANUSCRIPT_V2_EVOLUTION.md", warn = FALSE)
mstext <- paste(ms, collapse = " ")

checks <- list()

# helper: verify a manuscript claim against a recomputed value
report <- function(label, in_ms, recomputed, tol = 0.05) {
  ok <- if (is.numeric(recomputed) && is.numeric(in_ms)) abs(in_ms - recomputed) <= tol else identical(in_ms, recomputed)
  cat(sprintf("%-55s ms=%s recomputed=%s %s\n", label, in_ms, round(recomputed, 3),
              ifelse(ok, "OK", "** MISMATCH **")))
  ok
}

cat("=== consistency audit ===\n")
b <- readRDS("results/kncc/dual_omics_bundle.rds")
h <- readRDS("results/hhs/hhs_analysis_ready.rds")
kdelta <- b$prot_delta; hd <- h$delta

# 1. spliceosome decline 105/110
src_sp <- read.csv("results/kncc/paired_differential_proteins.csv")
n_sig <- sum(src_sp$is_spliceosome & src_sp$q < 0.05, na.rm = TRUE)
n_sp <- sum(src_sp$is_spliceosome)
report("spliceosome proteins FDR<0.05 (105/110)", 105, n_sig, 0)

# 2. subcomplex KW p=0.29
kw <- readLines("results/kncc/final_adjudication.txt")
# 3. PC2 HHS p
pc2 <- read.csv("results/hhs/hhs_pc2_effect_sizes.csv")
report("HHS PC2 mean (-0.95)", -0.95, pc2$median[pc2$cohort == "HHS"], 0.05)
report("HHS PC2 p (2.5e-04)", 2.5e-04, pc2$wilcox_p[pc2$cohort == "HHS"], 1e-3)

# 4. module effect sizes from archive
me <- read.csv("results/hhs/hhs_module_effect_sizes.csv")
get <- function(m, col) me[me$module == m, col]
report("IL2STAT5 KNCC mean (0.045)", 0.045, get("HALLMARK_IL2_STAT5_SIGNALING", "kncc_med"), 0.01)
report("IL2STAT5 HHS mean (0.197)", 0.197, get("HALLMARK_IL2_STAT5_SIGNALING", "hhs_med"), 0.02)
report("ROS KNCC mean (0.120)", 0.120, get("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "kncc_med"), 0.01)
report("ROS HHS mean (0.185)", 0.185, get("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "hhs_med"), 0.02)
report("EMT KNCC mean (-0.049)", -0.049, get("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "kncc_med"), 0.01)
report("EMT HHS mean (0.528)", 0.528, get("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "hhs_med"), 0.02)
report("NEURONAL KNCC (0.45)", 0.45, get("NEURONAL_CURATED", "kncc_med"), 0.02)
report("NEURONAL HHS (-0.69)", -0.69, get("NEURONAL_CURATED", "hhs_med"), 0.03)

# 5. GLASS result
gl <- read.csv("results/kncc/glass_second_axis.csv")
cat("GLASS module p-values (all n.s. expected):", paste(round(gl$paired_p, 3), collapse = ", "), "\n")

# 6. GeoMx EMT stain contrast
gm <- readLines("results/geomx/spatial_second_axis.txt")
emt_line <- gm[grepl("EMT.*CD64-GFAP", gm)]
cat("GeoMx EMT contrast in manuscript: -0.13, p=0.020 | in file:", emt_line, "\n")

# 7. tech ratio range claimed 7.5-260
tv <- read.delim("results/kncc/techvariance_direct.txt", sep = "\n", header = FALSE)
ratios <- as.numeric(gsub(".* ([0-9.]+) .*", "\\1", grep("^[A-Z]", tv$V1[20:40], value = TRUE)))
cat("tech ratios in file vs claimed range 7.5-260: check techvariance_direct.txt manually\n")

# 8. check figure references consistency (Fig. S1 defined?)
cat("\nFig references in ms:", paste(unique(unlist(regmatches(mstext, gregexpr("Fig\\. [0-9]+(\\. S?[0-9]+)?", mstext)))), collapse = ", "), "\n")
cat("S1 defined in legends:", grepl("Fig. S1", mstext), "\n")
