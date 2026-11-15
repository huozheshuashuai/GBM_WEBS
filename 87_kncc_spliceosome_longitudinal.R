# =============================================================================
# 87. KNCC module 2: spliceosome proteome remodeling in paired T1->T2
#     paired Wilcoxon per protein; complex-level score; specificity controls
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })
set.seed(42)

obj <- readRDS("results/kncc/proteome_analysis_ready.rds")
pm <- obj$pm; cli <- obj$annot; delta <- obj$delta; paired_pts <- obj$paired

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
sp <- intersect(splice, rownames(pm))
cat("spliceosome proteins in KNCC matrix:", length(sp), "\n")

# ---- per-protein paired test on delta (T2 - T1) ----
# Wilcoxon signed-rank on delta != 0; also report median delta and direction
res <- t(sapply(rownames(delta), function(g) {
  d <- delta[g, ]
  d <- d[!is.na(d)]
  if (length(d) < 30) return(c(n = length(d), med = NA, p = NA))
  wt <- suppressWarnings(wilcox.test(d, mu = 0))
  c(n = length(d), med = median(d), p = wt$p.value)
}))
res <- data.frame(gene = rownames(delta), res, row.names = NULL)
res$q <- p.adjust(res$p, "BH")
res$is_spliceosome <- res$gene %in% sp
write.csv(res, "results/kncc/paired_differential_proteins.csv", row.names = FALSE)

sig_sp <- res %>% filter(is_spliceosome, q < 0.05)
cat("\n=== spliceosome paired T2-T1 results ===\n")
cat("significant (FDR<0.05):", nrow(sig_sp), "/", sum(res$is_spliceosome), "\n")
cat("  up in T2:", sum(sig_sp$med > 0), "| down in T2:", sum(sig_sp$med < 0), "\n")
cat("direction among ALL spliceosome: up", sum(res$med[res$is_spliceosome] > 0, na.rm = TRUE),
    "/ down", sum(res$med[res$is_spliceosome] < 0, na.rm = TRUE), "\n")
cat("median delta spliceosome:", round(median(res$med[res$is_spliceosome], na.rm = TRUE), 3),
    " vs genome:", round(median(res$med[!res$is_spliceosome], na.rm = TRUE), 3), "\n")
wt_dir <- wilcox.test(res$med[res$is_spliceosome], res$med[!res$is_spliceosome])
cat("spliceosome vs genome delta distribution: p =", signif(wt_dir$p.value, 3), "\n")

# ---- complex-level score per sample (mean of spliceosome log-ratios) ----
score <- colMeans(pm[sp, , drop = FALSE], na.rm = TRUE)
cli$sp_score <- score[cli$sample]
t1 <- cli$sp_score[cli$timepoint == "T1"][match(paired_pts, cli$patient[cli$timepoint == "T1"])]
t2 <- cli$sp_score[cli$timepoint == "T2"][match(paired_pts, cli$patient[cli$timepoint == "T2"])]
wt_score <- wilcox.test(t2, t1, paired = TRUE)
cat("\nspliceosome complex score T1 vs T2 (paired, n =", length(paired_pts), "):\n")
cat("  median T1:", round(median(t1), 3), " T2:", round(median(t2), 3),
    " paired Wilcoxon p =", signif(wt_score$p.value, 3), "\n")

# ---- specificity: strict matched permutation (same scheme as fixed script 80) ----
# delta-score per patient for random abundance+variance-matched non-spliceosome sets
expr_mean <- rowMeans(pm, na.rm = TRUE)
expr_sd <- apply(pm, 1, sd, na.rm = TRUE)
dec <- cut(expr_mean, breaks = quantile(expr_mean, probs = 0:10/10, na.rm = TRUE),
           include.lowest = TRUE, labels = FALSE)
names(dec) <- rownames(pm)
bin <- integer(nrow(pm)); names(bin) <- rownames(pm)
for (d in 1:10) {
  idx <- which(dec == d)
  vt <- cut(expr_sd[idx], breaks = unique(quantile(expr_sd[idx], probs = 0:3/3, na.rm = TRUE)),
            include.lowest = TRUE, labels = FALSE)
  bin[idx] <- (d - 1) * 3 + as.integer(vt)
}
pool <- setdiff(rownames(pm), sp)
sp_bins <- bin[sp]
bin_levels <- sort(unique(sp_bins))
bin_pool <- lapply(bin_levels, function(b) which(bin == b & names(bin) %in% pool))
names(bin_pool) <- bin_levels
dec_pool <- lapply(1:10, function(d) which(dec == d & names(dec) %in% pool))
bin_dec <- ceiling(bin_levels / 3)
n_need <- sapply(bin_levels, function(b) sum(sp_bins == b))

sp_delta_score <- colMeans(delta[sp, , drop = FALSE], na.rm = TRUE)
obs_stat <- median(sp_delta_score, na.rm = TRUE)

n_iter <- 10000
rand_stat <- numeric(n_iter)
for (i in seq_len(n_iter)) {
  used <- logical(nrow(pm)); gs <- integer(0)
  for (k in seq_along(bin_levels)) {
    cand <- bin_pool[[k]]; avail <- cand[!used[cand]]; need <- n_need[k]
    if (length(avail) >= need) {
      pick <- sample(avail, need)
    } else {
      pick <- avail
      cand2 <- dec_pool[[bin_dec[k]]]
      avail2 <- cand2[!used[cand2] & !(cand2 %in% pick)]
      pick <- c(pick, sample(avail2, need - length(avail)))
    }
    used[pick] <- TRUE; gs <- c(gs, pick)
  }
  rand_stat[i] <- median(colMeans(delta[gs, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
}
p_emp <- (sum(abs(rand_stat) >= abs(obs_stat)) + 1) / (n_iter + 1)
cat("\nstrict matched permutation (10,000 sets): observed median delta-score =",
    round(obs_stat, 4), " empirical p =", signif(p_emp, 3), "\n")

# ---- figures ----
res$cat <- ifelse(res$is_spliceosome & res$q < 0.05, "Spliceosome FDR<0.05",
           ifelse(res$is_spliceosome, "Spliceosome n.s.", "Other"))
res <- res[order(res$cat != "Spliceosome FDR<0.05"), ]
res$sz <- ifelse(res$cat == "Spliceosome FDR<0.05", 1.8, 0.4)
res$shp <- ifelse(res$is_spliceosome, 17, 16)
p_vol <- ggplot(res, aes(med, -log10(p))) +
  geom_point(aes(color = cat, size = sz, shape = shp), alpha = 0.65) +
  scale_size_identity() + scale_shape_identity() +
  scale_color_manual(values = c("Spliceosome FDR<0.05" = "#7B4EA3",
                                "Spliceosome n.s." = "#3D9B8C", "Other" = "grey85")) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Median paired delta (T2 − T1, log2)", y = "-log10(p)",
       title = sprintf("Paired T2 vs T1: spliceosome %d/%d FDR<0.05",
                       nrow(sig_sp), sum(res$is_spliceosome)))
ggsave("figures_kncc/m2_paired_volcano.pdf", p_vol, width = 6, height = 5)

score_df <- data.frame(patient = rep(paired_pts, 2),
                       timepoint = rep(c("T1", "T2"), each = length(paired_pts)),
                       score = c(t1, t2))
p_sc <- ggplot(score_df, aes(timepoint, score)) +
  geom_line(aes(group = patient), color = "grey80") +
  geom_boxplot(aes(fill = timepoint), alpha = 0.85, outlier.shape = NA) +
  scale_fill_manual(values = c("T1" = "#1f77b4", "T2" = "#d62728")) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "Spliceosome protein score",
       title = sprintf("Complex score T1 vs T2 (paired p = %s)", signif(wt_score$p.value, 2)))
ggsave("figures_kncc/m2_complex_score_paired.pdf", p_sc, width = 4.5, height = 5)

null_df <- data.frame(x = rand_stat)
p_null <- ggplot(null_df, aes(x)) +
  geom_histogram(bins = 50, fill = "grey75") +
  geom_vline(xintercept = obs_stat, color = "#7B4EA3", linewidth = 1) +
  theme_bw(base_size = 11) +
  labs(x = "Median delta-score of matched random sets", y = "Count",
       title = sprintf("Spliceosome observed = %.3f (empirical p = %s)",
                       obs_stat, signif(p_emp, 2)))
ggsave("figures_kncc/m2_permutation_null.pdf", p_null, width = 6, height = 4.5)

sink("results/kncc/module2_results.txt", split = TRUE)
cat("=== Module 2: spliceosome longitudinal remodeling (KNCC paired) ===\n\n")
cat("spliceosome proteins tested:", sum(res$is_spliceosome), "\n")
cat("FDR<0.05:", nrow(sig_sp), "( up:", sum(sig_sp$med > 0), " down:", sum(sig_sp$med < 0), ")\n")
cat("direction all spliceosome: up", sum(res$med[res$is_spliceosome] > 0, na.rm = TRUE),
    "/ down", sum(res$med[res$is_spliceosome] < 0, na.rm = TRUE), "\n")
cat("median delta: spliceosome", round(median(res$med[res$is_spliceosome], na.rm = TRUE), 3),
    "vs genome", round(median(res$med[!res$is_spliceosome], na.rm = TRUE), 3),
    " (Wilcoxon p =", signif(wt_dir$p.value, 3), ")\n")
cat("complex score paired p =", signif(wt_score$p.value, 3), "\n")
cat("strict permutation empirical p =", signif(p_emp, 3), "\n\n")
cat("top 20 spliceosome proteins by |delta|:\n")
top <- res %>% filter(is_spliceosome) %>% arrange(desc(abs(med))) %>% head(20)
print(top[, c("gene", "med", "q")], row.names = FALSE)
sink()
cat("\nDone.\n")
