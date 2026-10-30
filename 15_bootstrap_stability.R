# =============================================================================
# 15. Bootstrap Stability of LASSO-Cox Feature Selection (six-gene signature)
# =============================================================================
# 100 bootstrap resamples of the training cohort; record how often each gene
# is selected (non-zero coefficient). Demonstrates the signature is not an
# overfitting artifact.
# Input:  data/raw/CGGA/ expression matrix + clinical (see CONFIG below)
# Output: results/bootstrap_selection_frequency.csv, figures/figS_bootstrap.pdf
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(glmnet)
  library(survival)
  library(ggplot2)
  library(dplyr)
})
set.seed(42)

# ---- CONFIG: adjust to the actual CGGA files after download ----
EXPR_FILE <- "data/raw/CGGA/cgga693_expr.rds"      # genes x samples matrix
CLIN_FILE <- "data/raw/CGGA/cgga693_clinical.rds"  # data.frame: sample, OS.time (months), OS (0/1)
# Candidate gene pool: lactylation / hypoxia / invasion biology (as in Methods)
CANDIDATE_GENES <- c(
  "EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP",
  "HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7",
  "HIF1A", "CA9", "VEGFA", "LDHA", "LDHB", "PKM", "MCT4", "ENO1", "ENO2", "GAPDH",
  "MMP2", "ITGA5", "SPP1", "FN1"
)
N_BOOT <- 100

if (!file.exists(EXPR_FILE) || !file.exists(CLIN_FILE)) {
  stop("CGGA training data not found. Download from http://www.cgga.org.cn/ ",
       "and place/prepare as specified in CONFIG at the top of this script.")
}

expr <- readRDS(EXPR_FILE)
clin <- readRDS(CLIN_FILE)

common <- intersect(colnames(expr), clin$sample)
expr <- expr[, common]
clin <- clin[match(common, clin$sample), ]

genes <- intersect(CANDIDATE_GENES, rownames(expr))
x_full <- t(expr[genes, ])          # samples x genes
y_full <- Surv(clin$OS.time, clin$OS)

# ---- Reference fit on full data ----
cv_ref <- cv.glmnet(x_full, y_full, family = "cox", alpha = 1, nfolds = 10, type.measure = "C")
fit_ref <- glmnet(x_full, y_full, family = "cox", alpha = 1, lambda = cv_ref$lambda.1se)
ref_genes <- rownames(coef(fit_ref))[as.vector(coef(fit_ref) != 0)]
cat("Reference model genes:", paste(ref_genes, collapse = ", "), "\n")

# ---- Bootstrap ----
selection <- matrix(0, nrow = N_BOOT, ncol = length(genes),
                    dimnames = list(NULL, genes))

for (b in seq_len(N_BOOT)) {
  idx <- sample(seq_len(nrow(x_full)), replace = TRUE)
  xb <- x_full[idx, , drop = FALSE]
  yb <- y_full[idx]
  fit <- tryCatch({
    cv <- cv.glmnet(xb, yb, family = "cox", alpha = 1, nfolds = 10, type.measure = "C")
    f <- glmnet(xb, yb, family = "cox", alpha = 1, lambda = cv$lambda.1se)
    as.vector(coef(f))
  }, error = function(e) NULL)
  if (!is.null(fit)) selection[b, ] <- as.integer(fit != 0)
  if (b %% 20 == 0) message("bootstrap ", b, "/", N_BOOT)
}

freq <- colMeans(selection)
freq_df <- data.frame(gene = names(freq), selection_frequency = freq,
                      in_reference = names(freq) %in% ref_genes) %>%
  arrange(desc(selection_frequency))

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
write.csv(freq_df, "results/bootstrap_selection_frequency.csv", row.names = FALSE)

sink("results/bootstrap_stability_summary.txt")
cat("=== Bootstrap stability of LASSO-Cox feature selection ===\n")
cat("N bootstrap:", N_BOOT, "\n")
cat("Reference genes:", paste(ref_genes, collapse = ", "), "\n\n")
cat("Selection frequency of reference genes:\n")
print(freq_df[freq_df$in_reference, ])
cat("\nTop 15 genes overall:\n")
print(head(freq_df, 15))
sink()

p <- ggplot(freq_df, aes(x = reorder(gene, selection_frequency),
                         y = selection_frequency, fill = in_reference)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#7f7f7f"),
                    labels = c("Not in final model", "In final six-gene model")) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Selection frequency (100 bootstraps)",
       fill = NULL,
       title = "Bootstrap stability of LASSO-Cox feature selection")

pdf("figures/figS_bootstrap.pdf", width = 7, height = 6)
print(p)
dev.off()

message("Done. See results/bootstrap_selection_frequency.csv and figures/figS_bootstrap.pdf")
