# =============================================================================
# 38. Confounder adjustment + clinical subgroups + treatment interaction
#     (six-gene signature, CGGA-693 & CGGA-325)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(dplyr); library(ggplot2) })

PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

load_cohort <- function(tag) {
  expr <- readRDS(paste0("data/raw/CGGA/", tag, "_expr.rds"))
  clin <- readRDS(paste0("data/raw/CGGA/", tag, "_clinical.rds"))
  # merge in MGMT / treatment fields from the original clinical file
  clin_file <- if (tag == "cgga693") "data/raw/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt"
               else "data/raw/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt"
  raw <- read.delim(clin_file, check.names = FALSE)
  colnames(raw)[1] <- "sample"
  extra <- raw[, c("sample", grep("MGMT|Radio|Chemo", colnames(raw), value = TRUE))]
  clin <- merge(clin, extra, by = "sample", all.x = TRUE)
  g <- intersect(names(PAPER_COEFS), rownames(expr))
  rs <- as.vector(t(expr[g, , drop = FALSE]) %*% PAPER_COEFS[g])
  clin$risk <- rs[match(clin$sample, colnames(expr))]
  clin <- clin[!is.na(clin$risk), ]
  clin$risk_group <- factor(ifelse(clin$risk > median(clin$risk), "High", "Low"),
                            levels = c("Low", "High"))
  clin
}

dir.create("results", showWarnings = FALSE)
sink("results/confounder_subgroup_summary.txt", split = TRUE)

for (tag in c("cgga693", "cgga325")) {
  clin <- load_cohort(tag)
  cat("\n##########", tag, "(n =", nrow(clin), ") ##########\n")

  # ---- 1. Multivariable Cox (confounder adjustment) ----
  clin$Age_s <- scale(clin$Age)
  clin$IDHwt <- as.integer(clin$IDH_mutation == "Wildtype")
  clin$MGMTm <- as.integer(tolower(clin$MGMTp_methylation_status) == "methylated")
  clin$GradeIV <- as.integer(clin$Grade == "WHO IV")
  d <- clin[complete.cases(clin[, c("Age_s", "IDHwt", "MGMTm", "GradeIV")]), ]
  m_uni <- coxph(Surv(OS.time, OS) ~ risk, data = d)
  m_multi <- coxph(Surv(OS.time, OS) ~ risk + Age_s + Gender + IDHwt + MGMTm + GradeIV, data = d)
  cat("\n-- Confounder adjustment --\n")
  cat("Univariate:      HR =", round(exp(coef(m_uni)["risk"]), 3),
      " p =", signif(summary(m_uni)$coefficients["risk", 5], 3), "\n")
  cat("Multivariable:   HR =", round(exp(coef(m_multi)["risk"]), 3),
      " p =", signif(summary(m_multi)$coefficients["risk", 5], 3),
      " (adjusted for age, gender, IDH, MGMT, grade)\n")

  # ---- 2. Subgroups ----
  cat("\n-- Subgroup HRs (per-unit risk) --\n")
  subgroups <- list(
    "IDH-wildtype" = d$IDHwt == 1,
    "IDH-mutant" = d$IDHwt == 0,
    "MGMT methylated" = d$MGMTm == 1,
    "MGMT unmethylated" = d$MGMTm == 0,
    "Primary" = d$PRS_type == "Primary",
    "Recurrent" = d$PRS_type == "Recurrent",
    "WHO grade IV" = d$GradeIV == 1
  )
  for (nm in names(subgroups)) {
    dd <- d[subgroups[[nm]], ]
    if (nrow(dd) >= 40 && sum(dd$OS) >= 10) {
      m <- coxph(Surv(OS.time, OS) ~ risk, data = dd)
      cat(sprintf("  %-20s n=%4d  HR=%.3f  p=%.3g\n", nm, nrow(dd),
                  exp(coef(m)["risk"]), summary(m)$coefficients["risk", 5]))
    }
  }

  # ---- 3. Treatment interaction ----
  d2 <- clin[!is.na(clin$`Radio_status`) & !is.na(clin$`Chemo_status`), ]
  if ("Chemo_status" %in% colnames(d2)) {}
  chemo_col <- grep("Chemo", colnames(d2), value = TRUE)[1]
  radio_col <- grep("Radio", colnames(d2), value = TRUE)[1]
  d2$TMZ <- as.integer(d2[[chemo_col]] == 1)
  d2$RT <- as.integer(d2[[radio_col]] == 1)
  cat("\n-- Treatment interaction (risk × treatment) --\n")
  for (tx in c("TMZ", "RT")) {
    dd <- d2[!is.na(d2[[tx]]), ]
    if (length(unique(dd[[tx]])) == 2) {
      f <- as.formula(paste("Surv(OS.time, OS) ~ risk *", tx))
      m <- tryCatch(coxph(f, data = dd), error = function(e) NULL)
      if (!is.null(m)) {
        int_term <- grep("risk:", names(coef(m)), value = TRUE)
        cat(sprintf("  risk × %s interaction: HR = %.3f, p = %.3g\n",
                    tx, exp(coef(m)[int_term]), summary(m)$coefficients[int_term, 5]))
      }
    }
  }
}
sink()
cat(readLines("results/confounder_subgroup_summary.txt"), sep = "\n")
