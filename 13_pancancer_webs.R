# =============================================================================
# 13. Pan-Cancer WEBS Validation (TCGA)
# =============================================================================
# Purpose: Test whether WEBS stratifies survival beyond GBM (generalizability)
# Cohorts: TCGA-LUAD, TCGA-BRCA, TCGA-SKCM (+ TCGA-GBM as reference)
# Output: results/pancancer_webs_summary.txt, figures/figS_pancancer.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})
source("tools/WEBS/webs.R")

WRITER_GENES <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASER_GENES <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

CANCERS <- c("TCGA-GBM", "TCGA-LUAD", "TCGA-BRCA", "TCGA-SKCM")

# ---- 1. Download & prepare each cohort ----
prepare_cohort <- function(project) {
  cache <- file.path("data/raw", project, paste0(project, "_prepared.rds"))
  if (file.exists(cache)) return(readRDS(cache))

  dir.create(file.path("data/raw", project), recursive = TRUE, showWarnings = FALSE)

  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, directory = file.path("data/raw", project))
  se <- GDCprepare(query, directory = file.path("data/raw", project))

  # keep primary tumor (01); for SKCM also metastatic (06)
  keep_types <- if (project == "TCGA-SKCM") c("01", "06") else "01"
  st <- substr(colnames(se), 14, 15)
  se <- se[, st %in% keep_types]

  # TPM matrix, gene symbols as rows
  tpm <- assay(se, "tpm_unstrand")
  symbols <- rowData(se)$gene_name
  keep <- !is.na(symbols) & symbols != ""
  tpm <- tpm[keep, ]
  symbols <- symbols[keep]
  # collapse duplicate symbols by max mean
  ord <- order(rowMeans(tpm), decreasing = TRUE)
  tpm <- tpm[ord, ]; symbols <- symbols[ord]
  tpm <- tpm[!duplicated(symbols), ]
  rownames(tpm) <- symbols[!duplicated(symbols)]

  # Clinical
  clin <- as.data.frame(colData(se))
  os_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death, clin$days_to_last_follow_up)
  os_event <- ifelse(clin$vital_status == "Dead", 1, 0)
  clinical <- data.frame(
    sample = colnames(tpm),
    OS.time = as.numeric(os_time) / 30.44,  # months
    OS = os_event
  )
  clinical <- clinical[!is.na(clinical$OS.time) & clinical$OS.time > 0, ]
  tpm <- tpm[, clinical$sample]

  saveRDS(list(tpm = tpm, clinical = clinical), cache)
  list(tpm = tpm, clinical = clinical)
}

# ---- 2. Per-cohort WEBS + survival ----
analyze_cohort <- function(project, dat) {
  webs_df <- webs_compute_bulk(dat$tpm, WRITER_GENES, ERASER_GENES)
  clin <- dat$clinical
  clin$WEBS <- webs_df$WEBS[match(clin$sample, webs_df$sample)]
  q <- quantile(clin$WEBS, c(0.25, 0.75), na.rm = TRUE)
  clin$WEBS_group <- ifelse(clin$WEBS <= q[1], "LLS",
                            ifelse(clin$WEBS >= q[2], "LHS", "LDTS"))
  clin$WEBS_group <- factor(clin$WEBS_group, levels = c("LLS", "LDTS", "LHS"))

  # Cox: continuous WEBS
  cox <- coxph(Surv(OS.time, OS) ~ WEBS, data = clin)
  s <- summary(cox)

  # KM: LHS vs rest
  clin$binary <- ifelse(clin$WEBS_group == "LHS", "LHS (top 25%)", "Other")
  km_fit <- survfit(Surv(OS.time, OS) ~ binary, data = clin)
  km_p <- surv_pvalue(km_fit, data = clin)$pval

  list(clinical = clin,
       HR = s$conf.int["WEBS", "exp(coef)"],
       HR_lower = s$conf.int["WEBS", "lower .95"],
       HR_upper = s$conf.int["WEBS", "upper .95"],
       cox_p = s$coefficients["WEBS", "Pr(>|z|)"],
       km_p = km_p)
}

results <- list()
for (proj in CANCERS) {
  message("=== ", proj, " ===")
  dat <- tryCatch(prepare_cohort(proj), error = function(e) { message(e$message); NULL })
  if (is.null(dat)) next
  results[[proj]] <- analyze_cohort(proj, dat)
  message(proj, ": n=", nrow(results[[proj]]$clinical),
          " HR=", round(results[[proj]]$HR, 3),
          " p=", signif(results[[proj]]$cox_p, 3))
}

# ---- 3. Summary table + forest plot + KM panels ----
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

sumdf <- do.call(rbind, lapply(names(results), function(p) {
  r <- results[[p]]
  data.frame(cohort = p, n = nrow(r$clinical),
             HR = r$HR, lower = r$HR_lower, upper = r$HR_upper,
             cox_p = r$cox_p, km_p = r$km_p)
}))
write.csv(sumdf, "results/pancancer_webs_summary.csv", row.names = FALSE)

sink("results/pancancer_webs_summary.txt")
cat("=== Pan-cancer WEBS survival analysis (Cox PH, continuous WEBS) ===\n\n")
print(sumdf)
sink()

# Forest plot
p_forest <- ggplot(sumdf, aes(x = cohort, y = HR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(title = "WEBS and overall survival across TCGA cohorts",
       y = "Hazard Ratio (per unit WEBS)", x = NULL)

# KM curves per cohort
km_plots <- lapply(names(results), function(p) {
  clin <- results[[p]]$clinical
  clin$binary <- ifelse(clin$WEBS_group == "LHS", "LHS (top 25%)", "Other")
  fit <- survfit(Surv(OS.time, OS) ~ binary, data = clin)
  ggsurvplot(fit, data = clin, pval = TRUE, title = p,
             palette = c("#d62728", "#7f7f7f"),
             xlab = "Time (months)", ylab = "Overall survival",
             legend.labs = c("LHS (top 25%)", "Other"),
             ggtheme = theme_bw(base_size = 11))
})

pdf("figures/figS_pancancer.pdf", width = 12, height = 9)
print(p_forest)
for (kp in km_plots) print(kp)
dev.off()

message("Done. See results/pancancer_webs_summary.csv and figures/figS_pancancer.pdf")
