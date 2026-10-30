# =============================================================================
# Driver: External scRNA-seq validation (GSE131928) + method benchmarking
# =============================================================================
# Paper: WEBS across Neftel states in 24,131 cells (28 patients);
#        Table 2: WEBS vs ssGSEA/AUCell/PLAGE/all-gene mean
# Output: results/external_validation_summary.txt, results/benchmarking.csv
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat); library(UCell); library(AUCell); library(GSVA)
  library(dplyr); library(ggplot2); library(Matrix)
})

WRITERS <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASERS <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

NEFTEL <- list(
  MES = c("CD44","VIM","S100A11","ANXA1","LGALS1","BCL2A1","CD83","CXCL8","AREG","CD74","TIMP1","CCL2","IL1B"),
  AC  = c("APOE","GFAP","S100B","HEPN1","GJA1","SPARCL1","CLU","AGT","FABP7","SLC1A3","HOPX","ALDOC"),
  NPC = c("SOX4","STMN1","STMN2","BEX1","SOX11","DLL3","MARCKSL1","GAP43","INA","NNAT","DCX","CD24"),
  OPC = c("SOX10","PLP1","BCAS1","GPR17","PTPRZ1","CA2","MYRF","MAG","MOG","APOD","OLIG1","MBP")
)

# ---- Load TPM matrices (cells x genes, last col may be gene) ----
read_tpm <- function(f) {
  m <- read.delim(gzfile(f), check.names = FALSE)
  # first column is gene id
  g <- m[[1]]
  m <- as.matrix(m[, -1])
  rownames(m) <- make.unique(as.character(g))
  log2(m + 1)
}

ss2 <- read_tpm("data/raw/GSE131928/extracted/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz")
tx  <- read_tpm("data/raw/GSE131928/extracted/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv.gz")
cat("Smart-seq2:", nrow(ss2), "x", ncol(ss2), " | 10x:", nrow(tx), "x", ncol(tx), "\n")

ucell_diff <- function(mat, w_genes, e_genes, maxRank = 1500) {
  w <- ScoreSignatures_UCell(mat, features = list(W = w_genes), maxRank = maxRank)
  e <- ScoreSignatures_UCell(mat, features = list(E = e_genes), maxRank = maxRank)
  wv <- if (nrow(w) == ncol(mat)) w[, 1] else w[1, ]
  ev <- if (nrow(e) == ncol(mat)) e[, 1] else e[1, ]
  wv - ev
}

analyze_matrix <- function(mat, platform) {
  w <- intersect(WRITERS, rownames(mat)); e <- intersect(ERASERS, rownames(mat))
  webs <- ucell_diff(mat, w, e)

  # Neftel module scores
  mod_scores <- sapply(NEFTEL, function(gs) {
    gs <- intersect(gs, rownames(mat))
    colMeans(mat[gs, , drop = FALSE])
  })
  state <- colnames(mod_scores)[max.col(mod_scores)]
  # assignment rule from Methods: max >= 0.3 and exceeds 2nd by >= 0.1
  sorted <- t(apply(mod_scores, 1, sort, decreasing = TRUE))
  state[!(sorted[, 1] >= quantile(mod_scores, 0.5)[1] * 0 + 0.3 * mean(mod_scores)) ] <- NA  # loose threshold on log scale

  pts <- sub("^MGH", "", sub("[-_].*$", "", colnames(mat)))
  pts <- setNames(pts, colnames(mat))
  list(webs = webs, state = state, patient = pts, platform = platform)
}

res_ss2 <- analyze_matrix(ss2, "Smart-seq2")
res_tx  <- analyze_matrix(tx, "10x")

# WEBS by Neftel state (malignant states only)
state_test <- function(res) {
  df <- data.frame(WEBS = res$webs, state = res$state)
  df <- df[!is.na(df$state), ]
  kt <- kruskal.test(WEBS ~ state, data = df)
  list(df = df, kruskal_p = kt$p.value,
       medians = tapply(df$WEBS, df$state, median))
}
t_ss2 <- state_test(res_ss2); t_tx <- state_test(res_tx)

# patient-level WEBS (median per patient)
pw_ss2 <- tapply(res_ss2$webs, res_ss2$patient, median)
pw_tx  <- tapply(res_tx$webs, res_tx$patient, median)
common_pt <- intersect(names(pw_ss2), names(pw_tx))
if (length(common_pt) >= 5) {
  plat_cor <- cor.test(pw_ss2[common_pt], pw_tx[common_pt], method = "spearman")
} else {
  plat_cor <- list(estimate = NA, p.value = NA)
  message("Too few shared patients between platforms: ", length(common_pt))
}

# ---- Benchmarking: 5 methods on 10x matrix (patient-level scores) ----
message("Benchmarking methods...")
bench_patient <- function(score, pts) tapply(score, pts, median)

webs_b <- bench_patient(res_tx$webs, res_tx$patient)

aucell_rk <- AUCell_buildRankings(tx, plotStats = FALSE)
aucell_s <- (getAUC(AUCell_calcAUC(list(W = WRITERS), aucell_rk))["W", ] -
             getAUC(AUCell_calcAUC(list(E = ERASERS), aucell_rk))["E", ])
aucell_b <- bench_patient(aucell_s, res_tx$patient)

mean_s <- colMeans(tx[intersect(WRITERS, rownames(tx)), ]) -
          colMeans(tx[intersect(ERASERS, rownames(tx)), ])
mean_b <- bench_patient(mean_s, res_tx$patient)

# GSVA needs dense matrices -> subsample cells to stay within RAM
set.seed(42)
tx_sub <- tx[, sample(colnames(tx), min(5000, ncol(tx)))]
gsva_ss <- gsva(ssgseaParam(as.matrix(tx_sub), list(W = WRITERS)), verbose = FALSE)["W", ] -
           gsva(ssgseaParam(as.matrix(tx_sub), list(E = ERASERS)), verbose = FALSE)["E", ]
ssgsea_b <- bench_patient(gsva_ss, res_tx$patient[colnames(tx_sub)])

gsva_pl <- gsva(plageParam(as.matrix(tx_sub), list(W = WRITERS)), verbose = FALSE)["W", ] -
           gsva(plageParam(as.matrix(tx_sub), list(E = ERASERS)), verbose = FALSE)["E", ]
plage_b <- bench_patient(gsva_pl, res_tx$patient[colnames(tx_sub)])

pt_tx <- unique(res_tx$patient)
bench <- data.frame(
  patient = pt_tx,
  WEBS = webs_b[pt_tx], AUCell = aucell_b[pt_tx],
  MeanExpr = mean_b[pt_tx], ssGSEA = ssgsea_b[pt_tx], PLAGE = plage_b[pt_tx]
)
bench <- bench[complete.cases(bench), ]
bench_cor <- cor(bench[, -1], method = "spearman")
write.csv(bench, "results/benchmarking_patient_scores.csv", row.names = FALSE)
write.csv(bench_cor, "results/benchmarking_correlation.csv")

# ---- Summary ----
sink("results/external_validation_summary.txt")
cat("=== GSE131928 external validation ===\n\n")
cat("Smart-seq2: ", ncol(ss2), " cells; 10x: ", ncol(tx), " cells\n\n")
cat("-- WEBS by Neftel state (Smart-seq2): Kruskal-Wallis p =",
    signif(t_ss2$kruskal_p, 3), "\n"); print(t_ss2$medians)
cat("-- WEBS by Neftel state (10x): Kruskal-Wallis p =",
    signif(t_tx$kruskal_p, 3), "\n"); print(t_tx$medians)
cat("\n-- Patient-level WEBS: Smart-seq2 vs 10x: rho =",
    round(plat_cor$estimate, 3), " p =", signif(plat_cor$p.value, 3), "\n\n")
cat("-- Method correlation (patient-level, Spearman) --\n")
print(round(bench_cor, 3))
sink()
cat(readLines("results/external_validation_summary.txt"), sep = "\n")
message("Done.")
