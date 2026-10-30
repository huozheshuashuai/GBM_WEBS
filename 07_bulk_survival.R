# =============================================================================
# Bulk RNA-seq Survival Analysis (8 Cohorts)
# =============================================================================
# Cohorts: TCGA-GBM, CGGA-693, CGGA-156, GSE7696, GSE13041, REMBRANDT, GLASS
# Methods: K-M curves, Cox PH, Harrell's C-index, AUC(t) from timeROC
# =============================================================================

library(survival)
library(survminer)
library(timeROC)
library(rms)

# ---- 1. Compute WEBS for Bulk Cohorts ----
compute_bulk_webs <- function(expr_matrix, clinical_data) {
  source("tools/WEBS/webs.R")
  
  webs_df <- webs_compute_bulk(
    expr_matrix,
    writer_genes = c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP"),
    eraser_genes = c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")
  )
  
  # Merge with clinical
  clinical_data$WEBS <- webs_df$WEBS[match(clinical_data$sample, webs_df$sample)]
  clinical_data$WEBS_group <- ifelse(
    clinical_data$WEBS <= quantile(clinical_data$WEBS, 0.25, na.rm = TRUE), "LLS",
    ifelse(clinical_data$WEBS >= quantile(clinical_data$WEBS, 0.75, na.rm = TRUE), "LHS", "LDTS")
  )
  
  return(clinical_data)
}

# ---- 2. Kaplan-Meier Survival Curves ----
plot_km <- function(clinical_data, cohort_name, time_col = "OS.time", event_col = "OS") {
  
  fit <- survfit(as.formula(paste("Surv(", time_col, ",", event_col, ") ~ WEBS_group")),
                 data = clinical_data)
  
  p <- ggsurvplot(fit, data = clinical_data,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, ncensor.plot = FALSE,
                  conf.int = FALSE,
                  xlab = "Time (months)", ylab = "Overall Survival",
                  title = paste("K-M Survival:", cohort_name),
                  palette = c("#E41A1C", "#377EB8", "#4DAF4A"),
                  legend.title = "WEBS Group",
                  legend.labs = c("LDTS", "LHS", "LLS"),
                  ggtheme = theme_bw(base_size = 12))
  
  return(p)
}

# ---- 3. Cox Proportional Hazards ----
run_cox <- function(clinical_data, time_col = "OS.time", event_col = "OS", 
                    covariates = "WEBS") {
  
  formula <- as.formula(paste("Surv(", time_col, ",", event_col, ") ~", covariates))
  cox_fit <- coxph(formula, data = clinical_data)
  
  summary_cox <- summary(cox_fit)
  
  results <- data.frame(
    variable = names(cox_fit$coefficients),
    HR = exp(cox_fit$coefficients),
    lower_95 = summary_cox$conf.int[, "lower .95"],
    upper_95 = summary_cox$conf.int[, "upper .95"],
    p_value = summary_cox$coefficients[, "Pr(>|z|)"]
  )
  
  return(list(fit = cox_fit, table = results))
}

# ---- 4. Harrell's C-index ----
compute_cindex <- function(clinical_data, time_col = "OS.time", event_col = "OS",
                           predictor = "WEBS") {
  
  cindex <- survcomp::concordance.index(
    x = clinical_data[[predictor]],
    surv.time = clinical_data[[time_col]],
    surv.event = clinical_data[[event_col]],
    method = "noether"
  )
  
  return(data.frame(
    C_index = cindex$c.index,
    lower_95 = cindex$lower,
    upper_95 = cindex$upper,
    p_value = cindex$p.value
  ))
}

# ---- 5. Time-Dependent AUC (timeROC) ----
compute_time_auc <- function(clinical_data, times = c(12, 24, 36), 
                             time_col = "OS.time", event_col = "OS",
                             marker = "WEBS") {
  
  roc_results <- lapply(times, function(t) {
    roc <- timeROC(
      T = clinical_data[[time_col]],
      delta = clinical_data[[event_col]],
      marker = clinical_data[[marker]],
      cause = 1,
      times = t,
      iid = FALSE
    )
    
    data.frame(
      time = t,
      AUC = roc$AUC[2],
      marker = marker
    )
  })
  
  do.call(rbind, roc_results)
}

# ---- 6. Calibrate Cox Model ----
calibrate_cox <- function(cox_fit, clinical_data, time_col = "OS.time", event_col = "OS") {
  
  # Predict risk scores
  risk_score <- predict(cox_fit, type = "risk")
  clinical_data$risk_score <- risk_score
  
  # Calibration plot
  cal <- calibrate(cox_fit, u = 24, method = "boot", B = 200)
  
  return(cal)
}

# ---- 7. Compare Across Cohorts ----
meta_analysis_cohorts <- function(cohort_results_list) {
  # cohort_results_list: list of data frames with HR and 95% CI per cohort
  
  library(meta)
  
  # Combine all cohort results
  all_results <- do.call(rbind, cohort_results_list)
  
  # Random-effects meta-analysis
  meta_res <- metagen(
    TE = log(all_results$HR),
    seTE = (log(all_results$upper_95) - log(all_results$lower_95)) / (2 * 1.96),
    studlab = all_results$cohort,
    sm = "HR",
    method.tau = "REML"
  )
  
  return(meta_res)
}
