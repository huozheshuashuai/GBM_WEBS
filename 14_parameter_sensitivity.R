# =============================================================================
# 14. WEBS Parameter Sensitivity Analysis
# =============================================================================
# (a) UCell maxRank sweep: 0.05 -> 0.90, rank correlation with default WEBS
# (b) Leave-one-out and leave-two-out writer/eraser gene removal
# (c) Cross-platform consistency: 10x vs Smart-seq2 (GSE131928)
# Input: data/processed/GSE138794_webs_annotated.rds (from run_02_scrna.R)
# Output: results/sensitivity_*.csv, figures/figS_sensitivity.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat)
  library(UCell)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
source("tools/WEBS/webs.R")

udiff_from_out <- function(w, e, ncells) {
  wv <- if (nrow(w) == ncells) w[, 1] else w[1, ]
  ev <- if (nrow(e) == ncells) e[, 1] else e[1, ]
  wv - ev
}

WRITER_GENES <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASER_GENES <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

obj <- readRDS("data/processed/GSE138794_webs_annotated.rds")
# subsample for runtime (robustness metrics converge well at 8k cells)
set.seed(42)
if (ncol(obj) > 8000) obj <- subset(obj, cells = sample(colnames(obj), 8000))
expr <- GetAssayData(obj, layer = "data")
webs_ref <- obj$WEBS

webs_from_scores <- function(w, e) w - e

# ---- (a) maxRank sweep ----
message("=== (a) UCell maxRank sweep ===")
maxranks <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50, 0.75, 0.90)
sweep_res <- lapply(maxranks, function(mr) {
  w <- ScoreSignatures_UCell(expr, features = list(W = WRITER_GENES), maxRank = mr * 20000)
  e <- ScoreSignatures_UCell(expr, features = list(E = ERASER_GENES), maxRank = mr * 20000)
  webs <- udiff_from_out(w, e, ncol(expr))
  data.frame(
    maxRank = mr,
    spearman_vs_default = cor(webs, webs_ref, method = "spearman"),
    pearson_vs_default = cor(webs, webs_ref, method = "pearson")
  )
})
sweep_df <- do.call(rbind, sweep_res)
write.csv(sweep_df, "results/sensitivity_maxrank.csv", row.names = FALSE)
print(sweep_df)

# ---- (b) leave-one-out / leave-two-out ----
message("=== (b) Leave-one-out / leave-two-out ===")
all_genes <- c(setNames(WRITER_GENES, "writer"), rep(NA, 0))
combos <- list()
for (g in WRITER_GENES) combos[[paste0("minus_", g)]] <- list(w = setdiff(WRITER_GENES, g), e = ERASER_GENES)
for (g in ERASER_GENES) combos[[paste0("minus_", g)]] <- list(w = WRITER_GENES, e = setdiff(ERASER_GENES, g))
# leave-two-out (all writer pairs and all eraser pairs)
wpairs <- combn(WRITER_GENES, 2, simplify = FALSE)
epairs <- combn(ERASER_GENES, 2, simplify = FALSE)
for (pr in wpairs) combos[[paste0("minus_", paste(pr, collapse = "_"))]] <- list(w = setdiff(WRITER_GENES, pr), e = ERASER_GENES)
for (pr in epairs) combos[[paste0("minus_", paste(pr, collapse = "_"))]] <- list(w = WRITER_GENES, e = setdiff(ERASER_GENES, pr))

loo_res <- lapply(names(combos), function(nm) {
  cfg <- combos[[nm]]
  w <- ScoreSignatures_UCell(expr, features = list(W = cfg$w))
  e <- ScoreSignatures_UCell(expr, features = list(E = cfg$e))
  webs <- udiff_from_out(w, e, ncol(expr))
  data.frame(
    perturbation = nm,
    n_writers = length(cfg$w),
    n_erasers = length(cfg$e),
    spearman_vs_default = cor(webs, webs_ref, method = "spearman")
  )
})
loo_df <- do.call(rbind, loo_res)
loo_df$type <- ifelse(grepl("_", gsub("^minus_", "", loo_df$perturbation)), "leave-two-out", "leave-one-out")
write.csv(loo_df, "results/sensitivity_leaveout.csv", row.names = FALSE)
print(summary(loo_df$spearman_vs_default))

# ---- (c) Cross-platform consistency (GSE131928: Smart-seq2 vs 10x) ----
message("=== (c) Platform consistency (GSE131928) ===")
platform_res <- NULL
ss2_file <- "data/raw/GSE131928/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz"
tenx_file <- "data/raw/GSE131928/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz"

if (file.exists(ss2_file) && file.exists(tenx_file)) {
  read_tpm <- function(f) {
    m <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
    as.matrix(m)
  }
  ss2 <- read_tpm(ss2_file)
  tx <- read_tpm(tenx_file)

  # patient = column prefix (columns are like "MGH42_A01_..." cell IDs)
  extract_patient <- function(cn) sub("-.*$", "", cn)

  patient_webs <- function(mat) {
    wg <- intersect(WRITER_GENES, rownames(mat))
    eg <- intersect(ERASER_GENES, rownames(mat))
    w <- ScoreSignatures_UCell(mat, features = list(W = wg))
    e <- ScoreSignatures_UCell(mat, features = list(E = eg))
    webs <- udiff_from_out(w, e, ncol(expr))
    pts <- extract_patient(colnames(mat))
    tapply(webs, pts, median, na.rm = TRUE)
  }

  pw_ss2 <- patient_webs(ss2)
  pw_10x <- patient_webs(tx)
  common <- intersect(names(pw_ss2), names(pw_10x))
  platform_res <- data.frame(
    patient = common,
    WEBS_smartseq2 = pw_ss2[common],
    WEBS_10x = pw_10x[common]
  )
  ct <- cor.test(platform_res$WEBS_smartseq2, platform_res$WEBS_10x, method = "spearman")
  message("Patient-level WEBS correlation (Smart-seq2 vs 10x): rho = ",
          round(ct$estimate, 3), ", p = ", signif(ct$p.value, 3))
  write.csv(platform_res, "results/sensitivity_platform.csv", row.names = FALSE)
} else {
  message("GSE131928 TPM files not found; skipping platform comparison.")
}

# ---- Figures ----
dir.create("figures", showWarnings = FALSE)

p1 <- ggplot(sweep_df, aes(x = factor(maxRank), y = spearman_vs_default)) +
  geom_point(size = 3) + geom_line(group = 1) +
  ylim(0.9, 1.0) + theme_bw(base_size = 12) +
  labs(x = "UCell maxRank (fraction of genes)", y = "Spearman rho vs default WEBS",
       title = "(a) WEBS is robust to UCell maxRank")

p2 <- ggplot(loo_df, aes(x = spearman_vs_default, fill = type)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.6) +
  theme_bw(base_size = 12) +
  labs(x = "Spearman rho vs default WEBS", y = "Count",
       title = "(b) Leave-one/two-out robustness")

p3 <- if (!is.null(platform_res)) {
  ggplot(platform_res, aes(x = WEBS_smartseq2, y = WEBS_10x)) +
    geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE, color = "red") +
    theme_bw(base_size = 12) +
    labs(x = "Patient WEBS (Smart-seq2)", y = "Patient WEBS (10x)",
         title = "(c) Cross-platform consistency")
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "GSE131928 not downloaded") + theme_void()
}

pdf("figures/figS_sensitivity.pdf", width = 15, height = 5)
print(p1 + p2 + p3)
dev.off()

message("Done. See results/sensitivity_*.csv and figures/figS_sensitivity.pdf")
