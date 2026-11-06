# =============================================================================
# 46. Extended benchmark: WEBS vs AUCell/ssGSEA/PLAGE/mean/Z-score on
#     malignant-cell biology + survival C-index
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(UCell); library(AUCell); library(GSVA)
                   library(survival); library(dplyr) })
set.seed(42)

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")

# ---- 1. scRNA benchmark (GSE138794 malignant cells, subsample for speed) ----
mal <- readRDS("data/processed/GSE138794_malignant.rds")
set.seed(42)
mal <- subset(mal, cells = sample(colnames(mal), 5000))
expr <- GetAssayData(mal, layer = "data")
w <- intersect(WRITERS, rownames(expr)); e <- intersect(ERASERS, rownames(expr))

ucell_diff <- function(mat, maxRank = 1500) {
  wu <- ScoreSignatures_UCell(mat, features = list(W = w), maxRank = maxRank)
  eu <- ScoreSignatures_UCell(mat, features = list(E = e), maxRank = maxRank)
  wv <- if (nrow(wu) == ncol(mat)) wu[, 1] else wu[1, ]
  ev <- if (nrow(eu) == ncol(mat)) eu[, 1] else eu[1, ]
  wv - ev
}

cat("Computing 6 methods on", ncol(expr), "malignant cells...\n")
scores <- list()
scores$WEBS <- ucell_diff(expr)
scores$AUCell <- {
  rk <- AUCell_buildRankings(expr, plotStats = FALSE)
  getAUC(AUCell_calcAUC(list(W = w), rk))["W", ] - getAUC(AUCell_calcAUC(list(E = e), rk))["E", ]
}
m_dense <- as.matrix(expr)
scores$ssGSEA <- gsva(ssgseaParam(m_dense, list(W = w)), verbose = FALSE)["W", ] -
                 gsva(ssgseaParam(m_dense, list(E = e)), verbose = FALSE)["E", ]
scores$PLAGE <- gsva(plageParam(m_dense, list(W = w)), verbose = FALSE)["W", ] -
                gsva(plageParam(m_dense, list(E = e)), verbose = FALSE)["E", ]
scores$MeanExpr <- colMeans(m_dense[w, ]) - colMeans(m_dense[e, ])
scores$Zscore <- {
  zw <- scale(t(m_dense[w, ])); ze <- scale(t(m_dense[e, ]))
  rowMeans(zw) - rowMeans(ze)
}

# evaluation: correlation with malignant markers + with pseudotime
markers <- intersect(c("MKI67", "TOP2A", "NES", "SOX2"), rownames(expr))
mal$pt <- readRDS("data/processed/GSE138794_malignant_pt.rds")$pseudotime[match(colnames(expr), colnames(readRDS("data/processed/GSE138794_malignant_pt.rds")))]

eval_df <- do.call(rbind, lapply(names(scores), function(m) {
  s <- scores[[m]]
  data.frame(
    method = m,
    cor_MKI67 = cor(s, as.numeric(expr["MKI67", ]), method = "spearman"),
    cor_pseudotime = cor(s, mal$pt, method = "spearman", use = "complete.obs")
  )
}))
write.csv(eval_df, "results/benchmark_scRNA_extended.csv", row.names = FALSE)
print(eval_df)

# ---- 2. bulk survival C-index per method (TCGA + CGGA) ----
expr_t <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin_t <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")

bulk_methods <- list(
  WEBS = colMeans(log2(expr_t[intersect(w, rownames(expr_t)), ] + 1)) -
         colMeans(log2(expr_t[intersect(e, rownames(expr_t)), ] + 1)),
  MeanExpr = colMeans(log2(expr_t[intersect(w, rownames(expr_t)), ] + 1)) -
             colMeans(log2(expr_t[intersect(e, rownames(expr_t)), ] + 1)),
  ssGSEA = gsva(ssgseaParam(as.matrix(expr_t), list(W = w)), verbose = FALSE)["W", ] -
           gsva(ssgseaParam(as.matrix(expr_t), list(E = e)), verbose = FALSE)["E", ],
  PLAGE = gsva(plageParam(as.matrix(expr_t), list(W = w)), verbose = FALSE)["W", ] -
          gsva(plageParam(as.matrix(expr_t), list(E = e)), verbose = FALSE)["E", ],
  Zscore = { zw <- scale(t(expr_t[intersect(w, rownames(expr_t)), ]))
             ze <- scale(t(expr_t[intersect(e, rownames(expr_t)), ]))
             rowMeans(zw) - rowMeans(ze) }
)

ci_df <- do.call(rbind, lapply(names(bulk_methods), function(m) {
  clin_t$score <- bulk_methods[[m]][clin_t$sample]
  fit <- coxph(Surv(OS.time, OS) ~ score, data = clin_t[!is.na(clin_t$score), ])
  data.frame(method = m, C_index = summary(fit)$concordance[1])
}))
write.csv(ci_df, "results/benchmark_bulk_cindex.csv", row.names = FALSE)
print(ci_df)

sink("results/benchmark_extended_summary.txt", split = TRUE)
cat("=== Extended benchmark ===\n\nscRNA (malignant cells, n=5000):\n")
print(eval_df)
cat("\nBulk TCGA-GBM C-index:\n"); print(ci_df)
sink()
cat("Done.\n")
