# =============================================================================
# Six-Gene Risk Score (LASSO-Cox)
# =============================================================================
# Training: 693-sample CGGA cohort
# Genes: SIRT2, HDAC1, MMP2, EP300, KAT2A, HDAC2 (LASSO-selected)
# Validation: 7 external cohorts + IMvigor210
# =============================================================================

library(glmnet)
library(survival)
library(survminer)

# ---- 1. Prepare Training Data ----
prepare_training_data <- function(expr_matrix, clinical_data, time_col = "OS.time", event_col = "OS") {
  
  # Ensure samples match
  common_samples <- intersect(colnames(expr_matrix), clinical_data$sample)
  expr <- expr_matrix[, common_samples]
  
  clinical <- clinical_data[match(common_samples, clinical_data$sample), ]
  
  # Time in months, event = 1 (death), 0 (censored)
  y <- Surv(clinical[[time_col]], clinical[[event_col]])
  
  # Transpose: samples x genes
  x <- t(expr)
  
  return(list(x = x, y = y, clinical = clinical))
}

# ---- 2. LASSO-Cox with Cross-Validation ----
run_lasso_cox <- function(x, y, nfolds = 10, alpha = 1) {
  
  # Cross-validation to find optimal lambda
  cv_fit <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = alpha,  # 1 = LASSO, 0 = Ridge
    nfolds = nfolds,
    type.measure = "C"
  )
  
  # Best lambda (1-SE rule for parsimony)
  best_lambda <- cv_fit$lambda.1se
  
  # Fit final model
  final_fit <- glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = alpha,
    lambda = best_lambda
  )
  
  # Extract coefficients
  coefs <- as.matrix(coef(final_fit))
  coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]
  
  print("Selected genes and coefficients:")
  print(coefs)
  
  return(list(
    cv_fit = cv_fit,
    final_fit = final_fit,
    coefficients = coefs,
    lambda = best_lambda
  ))
}

# ---- 3. Compute Risk Score ----
compute_risk_score <- function(expr_matrix, coefs) {
  
  # Risk score = sum(beta_i * expression_i)
  genes <- rownames(coefs)
  genes_present <- intersect(genes, rownames(expr_matrix))
  
  if (length(genes_present) < length(genes)) {
    warning("Missing genes: ", paste(setdiff(genes, genes_present), collapse = ", "))
  }
  
  # Subset to present genes
  coefs_sub <- coefs[genes_present, , drop = FALSE]
  expr_sub <- expr_matrix[genes_present, ]
  
  # Linear predictor
  risk_score <- as.vector(t(expr_sub) %*% coefs_sub)
  names(risk_score) <- colnames(expr_matrix)
  
  return(risk_score)
}

# ---- 4. Stratify by Risk Score ----
stratify_risk <- function(risk_score, method = "median") {
  
  if (method == "median") {
    cutoff <- median(risk_score, na.rm = TRUE)
    group <- ifelse(risk_score > cutoff, "High", "Low")
  } else if (method == "tertile") {
    q <- quantile(risk_score, c(1/3, 2/3))
    group <- ifelse(risk_score <= q[1], "Low",
                    ifelse(risk_score >= q[2], "High", "Intermediate"))
  }
  
  return(factor(group, levels = c("Low", "Intermediate", "High")))
}

# ---- 5. Validation Cohort Analysis ----
validate_signature <- function(expr_matrix, clinical_data, coefs, 
                                cohort_name, time_col = "OS.time", event_col = "OS") {
  
  # Compute risk score
  risk_score <- compute_risk_score(expr_matrix, coefs)
  
  # Add to clinical
  clinical_data$risk_score <- risk_score[match(clinical_data$sample, names(risk_score))]
  clinical_data$risk_group <- stratify_risk(clinical_data$risk_score, method = "median")
  
  # Cox model
  cox_res <- run_cox(clinical_data, time_col = time_col, event_col = event_col, 
                     covariates = "risk_score")
  
  # C-index
  cindex <- compute_cindex(clinical_data, time_col = time_col, event_col = event_col,
                           predictor = "risk_score")
  
  # Time-dependent AUC
  auc <- compute_time_auc(clinical_data, times = c(12, 24, 36),
                          time_col = time_col, event_col = event_col,
                          marker = "risk_score")
  
  return(list(
    clinical = clinical_data,
    cox = cox_res,
    cindex = cindex,
    auc = auc
  ))
}

# ---- 6. IMvigor210 Validation (Immunotherapy) ----
validate_imvigor210 <- function(imvigor_expr, imvigor_clinical, coefs) {
  
  risk_score <- compute_risk_score(imvigor_expr, coefs)
  imvigor_clinical$risk_score <- risk_score[match(imvigor_clinical$sample, names(risk_score))]
  imvigor_clinical$risk_group <- stratify_risk(imvigor_clinical$risk_score)
  
  # Response prediction (CR/PR vs SD/PD)
  library(pROC)
  roc_obj <- roc(imvigor_clinical$response %in% c("CR", "PR"), 
                 imvigor_clinical$risk_score)
  
  return(list(
    clinical = imvigor_clinical,
    AUC = auc(roc_obj),
    CI = ci.auc(roc_obj)
  ))
}

# ---- 7. Calibration Plot ----
plot_calibration <- function(clinical_data, time_point = 24, 
                              time_col = "OS.time", event_col = "OS") {
  
  library(rms)
  
  # Fit Cox model
  f <- cph(as.formula(paste("Surv(", time_col, ",", event_col, ") ~ risk_score")),
           data = clinical_data, x = TRUE, y = TRUE, surv = TRUE)
  
  # Calibrate at specified time point
  cal <- calibrate(f, u = time_point, cmethod = "KM", m = 50, B = 200)
  
  plot(cal, xlab = "Predicted Probability", ylab = "Observed Probability",
       main = paste("Calibration at", time_point, "months"))
  abline(0, 1, lty = 2, col = "red")
  
  return(cal)
}

# ---- 8. Nomogram ----
build_nomogram <- function(clinical_data, time_col = "OS.time", event_col = "OS") {
  
  library(rms)
  
  dd <- datadist(clinical_data)
  options(datadist = "dd")
  
  f <- cph(as.formula(paste("Surv(", time_col, ",", event_col, ") ~ risk_score + age + gender")),
           data = clinical_data, x = TRUE, y = TRUE, surv = TRUE)
  
  nom <- nomogram(f, fun = list(function(x) surv(12, x), function(x) surv(24, x)),
                  funlabel = c("1-Year Survival", "2-Year Survival"))
  
  plot(nom)
  return(nom)
}
