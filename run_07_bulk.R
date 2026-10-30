# =============================================================================
# Driver: Bulk survival analysis across cohorts
# =============================================================================
# Cohorts: TCGA-GBM (LinkedOmics), CGGA-693, CGGA-325, GSE7696, GSE13041,
#          REMBRANDT (GSE108474)
# Output: results/bulk_survival_summary.csv + KM plots (figures/fig_survival_*)
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(survival)
  library(survminer)
  library(timeROC)
  library(GEOquery)
  library(dplyr)
})
source("tools/WEBS/webs.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

WRITERS <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASERS <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

analyze_bulk <- function(expr, clin, cohort) {
  # expr: genes x samples (log-scale); clin: sample, OS.time (months), OS
  w <- intersect(WRITERS, rownames(expr)); e <- intersect(ERASERS, rownames(expr))
  webs <- colMeans(expr[w, , drop = FALSE], na.rm = TRUE) -
          colMeans(expr[e, , drop = FALSE], na.rm = TRUE)
  clin$WEBS <- webs[match(clin$sample, colnames(expr))]
  clin <- clin[!is.na(clin$WEBS), ]
  q <- quantile(clin$WEBS, c(0.25, 0.75), na.rm = TRUE)
  clin$WEBS_group <- factor(ifelse(clin$WEBS <= q[1], "LLS",
                            ifelse(clin$WEBS >= q[2], "LHS", "LDTS")),
                            levels = c("LLS", "LDTS", "LHS"))

  cox <- summary(coxph(Surv(OS.time, OS) ~ WEBS, data = clin))
  fit <- survfit(Surv(OS.time, OS) ~ WEBS_group, data = clin)
  km_p <- surv_pvalue(fit, data = clin)$pval

  # KM plot
  p <- ggsurvplot(fit, data = clin, pval = TRUE, title = cohort,
                  palette = c("#1f77b4", "#ff7f0e", "#d62728"),
                  xlab = "Time (months)", ylab = "Overall survival",
                  legend.title = "WEBS", legend.labs = c("LLS", "LDTS", "LHS"))
  pdf(paste0("figures/fig_survival_", gsub("[^A-Za-z0-9]", "_", cohort), ".pdf"),
      width = 6, height = 5)
  print(p)
  dev.off()

  data.frame(cohort = cohort, n = nrow(clin), events = sum(clin$OS),
             HR = cox$conf.int["WEBS", "exp(coef)"],
             HR_lower = cox$conf.int["WEBS", "lower .95"],
             HR_upper = cox$conf.int["WEBS", "upper .95"],
             cox_p = cox$coefficients["WEBS", "Pr(>|z|)"],
             km_p = km_p)
}

results <- list()

# ---- TCGA-GBM (LinkedOmics) ----
expr <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_expr_log2rsem.rds")
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
results[["TCGA-GBM"]] <- analyze_bulk(expr, clin, "TCGA-GBM")

# ---- CGGA-693 / CGGA-325 (all glioma + IDH-wt GBM subset) ----
for (cg in c("cgga693", "cgga325")) {
  expr <- readRDS(paste0("data/raw/CGGA/", cg, "_expr.rds"))
  clin <- readRDS(paste0("data/raw/CGGA/", cg, "_clinical.rds"))
  results[[cg]] <- analyze_bulk(expr, clin, cg)
  # IDH-wildtype GBM subset
  sub <- clin$Histology %in% c("GBM", "rGBM") & clin$IDH_mutation == "Wildtype"
  sub_ids <- intersect(clin$sample[sub], colnames(expr))
  if (length(sub_ids) >= 40) {
    results[[paste0(cg, "_IDHwtGBM")]] <- analyze_bulk(
      expr[, sub_ids], clin[match(sub_ids, clin$sample), ], paste0(cg, " IDHwt-GBM"))
  }
}

# ---- GEO cohorts: GSE7696, GSE13041, REMBRANDT (GSE108474) ----
# Local series matrix files (downloaded via pget) + Bioc annotation packages
load_geo_local <- function(acc, gpl_pkg) {
  f <- list.files(file.path("data/raw", acc), pattern = "series_matrix.txt.gz$",
                  full.names = TRUE)[1]
  g <- getGEO(filename = f)
  exprs_g <- exprs(g)
  if (!require(gpl_pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(gpl_pkg, update = FALSE, ask = FALSE)
  }
  lib <- getExportedValue(gpl_pkg, paste0(sub("\\.db$", "", gpl_pkg), "SYMBOL"))
  sym <- unlist(as.list(lib))
  probe2sym <- sym[intersect(rownames(exprs_g), names(sym))]
  exprs_g <- exprs_g[names(probe2sym), ]
  rownames(exprs_g) <- unname(probe2sym)
  exprs_g <- exprs_g[!is.na(rownames(exprs_g)) & rownames(exprs_g) != "NA", ]
  # collapse duplicates
  o <- order(rowMeans(exprs_g, na.rm = TRUE), decreasing = TRUE)
  exprs_g <- exprs_g[o, ]
  exprs_g <- exprs_g[!duplicated(rownames(exprs_g)), ]
  list(expr = exprs_g, pheno = pData(g))
}

extract_survival_geo <- function(pd) {
  cn <- colnames(pd)
  time_col <- cn[grepl("overall.?survival|os\\b|survival.*day|days", cn, ignore.case = TRUE)][1]
  event_col <- cn[grepl("vital|death|event|status", cn, ignore.case = TRUE)][1]
  if (is.na(time_col) || is.na(event_col)) return(NULL)
  tim <- suppressWarnings(as.numeric(sub("[^0-9\\.].*$", "", as.character(pd[[time_col]]))))
  ev <- pd[[event_col]]
  if (is.numeric(ev)) {
    ev <- as.integer(ev)
  } else {
    ev <- as.integer(grepl("dead|deceased|1", as.character(ev), ignore.case = TRUE))
  }
  data.frame(sample = rownames(pd), OS.time = tim / 30.44, OS = ev)
}

extract_char <- function(pd, pattern) {
  cols <- grep("characteristics", colnames(pd), value = TRUE)
  for (cc in cols) {
    vals <- as.character(pd[[cc]])
    if (any(grepl(pattern, vals, ignore.case = TRUE), na.rm = TRUE)) {
      v <- sub(sprintf("^.*%s:? *", pattern), "", vals, ignore.case = TRUE)
      return(v)
    }
  }
  NULL
}

geo_cohorts <- list(
  list(acc = "GSE7696", gpl = "hgu133plus2.db",
       time_pat = "survival time in months", event_pat = "survival status",
       time_unit = "months", event_dead = "1"),
  list(acc = "GSE13041", gpl = "hgu133plus2.db",
       time_pat = "tts\\(days\\)", event_pat = "vital status",
       time_unit = "days", event_dead = "DECEASED")
)

for (cfg in geo_cohorts) {
  nm <- cfg$acc
  message("=== ", nm, " ===")
  dat <- tryCatch(load_geo_local(cfg$acc, cfg$gpl), error = function(e) { message(e$message); NULL })
  if (is.null(dat)) next
  tim <- suppressWarnings(as.numeric(extract_char(dat$pheno, cfg$time_pat)))
  ev <- as.character(extract_char(dat$pheno, cfg$event_pat))
  if (all(is.na(tim)) || all(is.na(ev))) { message(nm, ": survival fields missing"); next }
  clin <- data.frame(
    sample = rownames(dat$pheno),
    OS.time = ifelse(cfg$time_unit == "days", tim / 30.44, tim),
    OS = as.integer(ev == cfg$event_dead)
  )
  clin <- clin[!is.na(clin$OS.time) & !is.na(clin$OS) & clin$OS.time > 0, ]
  clin <- clin[clin$sample %in% colnames(dat$expr), ]
  if (nrow(clin) < 30) { message(nm, ": too few samples with survival"); next }
  results[[nm]] <- analyze_bulk(dat$expr[, clin$sample], clin, nm)
}

sumdf <- do.call(rbind, results)
write.csv(sumdf, "results/bulk_survival_summary.csv", row.names = FALSE)
print(sumdf)
message("Done.")
