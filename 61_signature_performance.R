# =============================================================================
# 61. Signature performance: HR+CI, C-index, time-dependent AUC (all cohorts)
# 62. 107/111 threshold sensitivity (different WEBS cutoffs)
# 63. GLASS recurrence: independent patients
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(timeROC); library(dplyr) })

PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
score_cohort <- function(expr_m) {
  g <- intersect(names(PAPER_COEFS), rownames(expr_m))
  as.vector(t(expr_m[g, , drop = FALSE]) %*% PAPER_COEFS[g])
}

sink("results/signature_performance_full.txt", split = TRUE)
cat("=== Six-gene signature: full performance report ===\n\n")
cat(sprintf("%-24s %4s %4s  %-18s %-9s %-7s %s\n",
            "cohort", "n", "ev", "HR (95% CI)", "Cox p", "C-index", "AUC(12/36/60mo)"))

report <- function(name, risk, time, event) {
  ok <- !is.na(time) & !is.na(event) & time > 0 & !is.na(risk)
  risk <- as.numeric(scale(risk[ok]))   # per-SD HR
  time <- as.numeric(time[ok]); event <- as.numeric(event[ok])
  fit <- coxph(Surv(time, event) ~ risk)
  co <- summary(fit)$coefficients; ci <- summary(fit)$conf.int
  cidx <- summary(fit)$concordance[1]
  auc <- tryCatch({
    tr <- timeROC(T = time, delta = event, marker = risk, cause = 1,
                  times = c(12, 36, 60), iid = FALSE)
    paste(round(tr$AUC, 2), collapse = "/")
  }, error = function(e) "NA")
  cat(sprintf("%-24s %4d %4d  %5.2f (%4.2f-%5.2f)   %.2e  %.3f   %s\n",
              name, length(risk), sum(event), ci[1, "exp(coef)"], ci[1, 3], ci[1, 4],
              co[1, "Pr(>|z|)"], cidx, auc))
}

# TCGA-GBM training (LinkedOmics cache)
e <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
c0 <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
report("TCGA-GBM (train)", score_cohort(as.matrix(e)), c0$OS.time, c0$OS)

# CGGA
e693 <- readRDS("data/raw/CGGA/cgga693_expr.rds"); c693 <- readRDS("data/raw/CGGA/cgga693_clinical.rds")
report("CGGA-693", score_cohort(as.matrix(e693)), c693$OS.time, c693$OS)
e325 <- readRDS("data/raw/CGGA/cgga325_expr.rds"); c325 <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
report("CGGA-325", score_cohort(as.matrix(e325)), c325$OS.time, c325$OS)

# GLASS
g <- read.csv("results/glass_primary_scores.csv")
report("GLASS primary", g$risk, g$OS.time, g$OS)

# GSE83300 (recompute quickly)
suppressMessages(library(GEOquery))
g8 <- getGEO(filename = "data/raw/GSE83300/GSE83300_series_matrix.txt.gz")
ex8 <- exprs(g8)
gpl <- getGEO("GPL6480", destdir = "data/raw/GPL")
fd <- Table(gpl); fd <- fd[!is.na(fd$ID) & fd$ID != "", ]; rownames(fd) <- fd$ID
sym <- fd$GENE_SYMBOL
cp <- intersect(rownames(ex8), rownames(fd))
ex8 <- ex8[cp, ]; sym <- sym[match(cp, rownames(fd))]
okp <- !is.na(sym) & sym != "" & sym != "---"
ex8 <- ex8[okp, ]; sym <- sym[okp]
o <- order(rowMeans(ex8, na.rm = TRUE), decreasing = TRUE)
ex8 <- ex8[o, ]; sym <- sym[o]
ex8 <- ex8[!duplicated(sym), ]; rownames(ex8) <- sym[!duplicated(sym)]
ex8 <- log2(ex8 + 1)
pd <- pData(g8)
cc <- grep("characteristics", colnames(pd), value = TRUE)
get_field <- function(pat) {
  for (c0 in cc) {
    v <- as.character(pd[[c0]])
    if (any(grepl(pat, v, ignore.case = TRUE))) return(sub("^.*?:\\s*", "", v))
  }
  NULL
}
tim <- as.numeric(get_field("overall survival"))
ev  <- as.integer(get_field("events"))
report("GSE83300", score_cohort(ex8), tim, ev)

cat("\n=== GLASS recurrence: independent patients ===\n")
lp <- read.csv("results/glass_recurrence_interval.csv")
cat("Pairs:", nrow(lp), " Independent patients:", length(unique(lp$case_barcode)), "\n")
dups <- sum(duplicated(lp$case_barcode))
cat("Patients contributing >1 pair:", dups, "\n")
sink()
cat(readLines("results/signature_performance_full.txt"), sep = "\n")
