# =============================================================================
# 106. Baseline prediction: can primary-tumor (T1) molecular state predict
#      recurrence (T2) endpoints? Baseline-adjusted models with nested CV.
#      M0 = clinical+technical | M1 = M0 + target at T1 | M2 = M1 + candidate T1 modules
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
pm <- b$prot_samples; cli <- b$prot_annot; master <- b$master

# ---- per-sample module scores at T1 and T2 ----
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(pm)))

mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(pm)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(pm)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)
mods <- mods[sapply(mods, length) >= 8]

score_mat <- t(sapply(mods, function(g) colMeans(pm[g, , drop = FALSE], na.rm = TRUE)))
colnames(score_mat) <- colnames(pm)

paired <- b$prot_paired
t1s <- cli$sample[cli$timepoint == "T1"][match(paired, cli$patient[cli$timepoint == "T1"])]
t2s <- cli$sample[cli$timepoint == "T2"][match(paired, cli$patient[cli$timepoint == "T2"])]
S1 <- score_mat[, t1s]; S2 <- score_mat[, t2s]
colnames(S1) <- colnames(S2) <- paired

# clinical covariates
clin <- master %>% filter(patient %in% paired) %>%
  mutate(AGE = as.numeric(AGE), SEX = as.factor(SEX))
clin <- clin[!is.na(clin$AGE) & !is.na(clin$SEX), ]
clin <- clin[match(clin$patient, clin$patient), ]
# restrict all matrices to patients with complete covariates
keep_pts <- clin$patient
S1 <- S1[, keep_pts]; S2 <- S2[, keep_pts]
paired <- keep_pts
C0 <- model.matrix(~ AGE + SEX, data = clin)
rownames(C0) <- clin$patient

targets <- c("NEURO", "ROS", "IL2STAT5", "PROLIF")
candidates <- setdiff(rownames(score_mat), NULL)

# ---- nested CV: outer 10-fold patient-level x 10 repeats; inner not needed
# (no feature selection; fixed small predictor sets) ----
cv_eval <- function(target) {
  y1 <- S1[target, ]; y2 <- S2[target, ]
  X1 <- t(S1[setdiff(candidates, target), , drop = FALSE])
  n <- length(paired)
  mse <- function(e) mean(e^2)
  err <- c(M0 = c(), M1 = c(), M2 = c())
  for (rep in 1:10) {
    folds <- sample(rep(1:10, length.out = n))
    for (k in 1:10) {
      te <- folds == k; tr <- !te
      # M0: clinical only
      f0 <- lm(y2[tr] ~ C0[tr, -1])
      e0 <- y2[te] - predict(f0, as.data.frame(C0[te, -1, drop = FALSE]))
      # M1: + target baseline
      d1 <- data.frame(y = y2[tr], y1 = y1[tr], C0[tr, -1, drop = FALSE])
      f1 <- lm(y ~ ., data = d1)
      d1t <- data.frame(y1 = y1[te], C0[te, -1, drop = FALSE])
      e1 <- y2[te] - predict(f1, d1t)
      # M2: + other module baselines
      d2 <- data.frame(y = y2[tr], y1 = y1[tr], X1[tr, , drop = FALSE], C0[tr, -1, drop = FALSE])
      f2 <- lm(y ~ ., data = d2)
      d2t <- data.frame(y1 = y1[te], X1[te, , drop = FALSE], C0[te, -1, drop = FALSE])
      e2 <- y2[te] - predict(f2, d2t)
      err$M0 <- c(err$M0, mse(e0)); err$M1 <- c(err$M1, mse(e1)); err$M2 <- c(err$M2, mse(e2))
    }
  }
  data.frame(target = target,
             M0 = mean(err$M0), M1 = mean(err$M1), M2 = mean(err$M2),
             gain_M2_vs_M1 = (mean(err$M1) - mean(err$M2)) / mean(err$M1))
}
res <- do.call(rbind, lapply(targets, cv_eval))

# permutation test for M2 gain (shuffle T2 endpoint across patients), 500x
perm_test <- function(target, n_perm = 500) {
  obs <- res$gain_M2_vs_M1[res$target == target]
  y1 <- S1[target, ]; y2 <- S2[target, ]
  X1 <- t(S1[setdiff(candidates, target), , drop = FALSE])
  n <- length(paired)
  perm_gain <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    y2p <- sample(y2)
    folds <- sample(rep(1:10, length.out = n))
    e1 <- e2 <- c()
    for (k in 1:10) {
      te <- folds == k; tr <- !te
      d1 <- data.frame(y = y2p[tr], y1 = y1[tr], C0[tr, -1, drop = FALSE])
      f1 <- lm(y ~ ., data = d1)
      e1 <- c(e1, (y2p[te] - predict(f1, data.frame(y1 = y1[te], C0[te, -1, drop = FALSE])))^2)
      d2 <- data.frame(y = y2p[tr], y1 = y1[tr], X1[tr, , drop = FALSE], C0[tr, -1, drop = FALSE])
      f2 <- lm(y ~ ., data = d2)
      e2 <- c(e2, (y2p[te] - predict(f2, data.frame(y1 = y1[te], X1[te, , drop = FALSE], C0[te, -1, drop = FALSE])))^2)
    }
    perm_gain[p] <- (mean(e1) - mean(e2)) / mean(e1)
  }
  (sum(perm_gain >= obs) + 1) / (n_perm + 1)
}
cat("running permutation tests...\n")
res$perm_p <- sapply(targets, perm_test)

sink("results/kncc/baseline_prediction.txt", split = TRUE)
cat("=== Baseline prediction (nested repeated 10-fold CV, 10 repeats) ===\n")
cat("M0 = age+sex | M1 = M0 + target T1 | M2 = M1 + other module T1 scores\n\n")
print(res, row.names = FALSE)
cat("\nInterpretation: gain_M2_vs_M1 > 0 with perm_p < 0.05 means other T1 module\n")
cat("states carry incremental information about the T2 endpoint beyond its own baseline.\n")
sink()

res$target <- factor(res$target, levels = targets)
p <- ggplot(res, aes(target, gain_M2_vs_M1, fill = perm_p < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = paste0("p=", signif(perm_p, 2))), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("TRUE" = "#7B4EA3", "FALSE" = "grey75"),
                    labels = c("n.s.", "perm p<0.05")) +
  theme_bw(base_size = 12) +
  labs(x = "T2 endpoint module", y = "CV MSE reduction (M2 vs M1)",
       fill = NULL, title = "Does other primary-state information predict recurrence endpoints?")
ggsave("figures_kncc/baseline_prediction.pdf", p, width = 6.5, height = 5)
write.csv(res, "results/kncc/baseline_prediction.csv", row.names = FALSE)
cat("Done.\n")
