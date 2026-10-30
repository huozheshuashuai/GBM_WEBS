# =============================================================================
# 19. SHAP Interpretability for the Six-Gene Risk Model (optional)
# =============================================================================
# SHAP values for the Cox linear predictor (risk score) using fastshap,
# showing each gene's contribution across risk strata.
# Input:  training cohort expression + the fitted glmnet model
#         (defaults reuse the CONFIG of 15_bootstrap_stability.R)
# Output: results/shap_values.csv, figures/figS_shap.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
if (!require("kernelshap", quietly = TRUE)) {
  install.packages("kernelshap", repos = "https://cloud.r-project.org/")
}
suppressMessages({
  library(kernelshap)
  library(glmnet)
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
set.seed(42)

EXPR_FILE <- "data/raw/CGGA/cgga693_expr.rds"
CLIN_FILE <- "data/raw/CGGA/cgga693_clinical.rds"
# Paper's final six-gene model coefficients (Methods):
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

if (!file.exists(EXPR_FILE)) stop("CGGA training data not found. See 15_bootstrap_stability.R CONFIG.")

expr <- readRDS(EXPR_FILE)
clin <- readRDS(CLIN_FILE)
common <- intersect(colnames(expr), clin$sample)
X <- t(expr[names(PAPER_COEFS), common])   # samples x 6 genes
clin <- clin[match(common, clin$sample), ]

# Refit six-gene Cox (fixed gene set, coefficients re-estimated on this cohort)
y <- Surv(clin$OS.time, clin$OS)
fit <- glmnet(X, y, family = "cox", alpha = 0, lambda = 1e-4)  # near-unpenalized refit

pred_fun <- function(object, newdata) {
  as.vector(predict(object, newx = as.matrix(newdata), type = "link"))
}

message("Computing SHAP values (kernelshap)...")
ks <- kernelshap(fit, X = as.data.frame(X), pred_fun = pred_fun, verbose = FALSE)
shap <- ks$S

shap_df <- as.data.frame(shap)
shap_df$sample <- rownames(X)
risk <- as.vector(predict(fit, newx = X, type = "link"))
shap_df$risk_score <- risk
shap_df$risk_group <- ifelse(risk > median(risk), "High", "Low")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write.csv(shap_df, "results/shap_values.csv", row.names = FALSE)

# Summary: mean |SHAP| per gene, overall and by risk group
long <- shap_df %>%
  pivot_longer(cols = all_of(names(PAPER_COEFS)), names_to = "gene", values_to = "shap")

summary_df <- long %>%
  group_by(gene, risk_group) %>%
  summarise(mean_abs_shap = mean(abs(shap)), mean_shap = mean(shap), .groups = "drop")
write.csv(summary_df, "results/shap_summary.csv", row.names = FALSE)

p1 <- ggplot(summary_df, aes(x = reorder(gene, mean_abs_shap), y = mean_abs_shap, fill = risk_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("High" = "#d62728", "Low" = "#1f77b4")) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Mean |SHAP|", fill = "Risk group",
       title = "Feature importance by risk group (SHAP)")

p2 <- ggplot(long, aes(x = shap, y = gene, color = risk_group)) +
  geom_violin(scale = "width") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("High" = "#d62728", "Low" = "#1f77b4")) +
  theme_bw(base_size = 12) +
  labs(x = "SHAP value (impact on risk score)", y = NULL, color = "Risk group",
       title = "SHAP distribution per gene")

pdf("figures/figS_shap.pdf", width = 12, height = 5)
print(p1 + p2)
dev.off()

message("Done. See results/shap_summary.csv and figures/figS_shap.pdf")
