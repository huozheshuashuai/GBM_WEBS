# Fast MIF permutation test: LLS vs LHS interaction score difference
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(Matrix) })
set.seed(42)

obj <- readRDS("data/processed/GSE138794_final.rds")
data_mat <- GetAssayData(obj, layer = "data")

is_mac <- (!is.na(obj$original_type) & obj$original_type == "Myeloid") |
          (is.na(obj$original_type) & obj$cell_type == "Macrophage_Microglia")
is_mal <- (!is.na(obj$original_type) & obj$original_type == "Neoplastic_cell") |
          (is.na(obj$original_type) & obj$cell_type == "Malignant")

v_mif <- as.numeric(data_mat["MIF", ])
v_rec <- as.numeric(colMeans(data_mat[c("CD74", "CXCR4"), ]))
grp <- obj$WEBS_group

mac_idx <- which(is_mac & grp %in% c("LLS", "LHS"))
mal_idx <- which(is_mal & grp %in% c("LLS", "LHS"))
grp_mac <- factor(as.character(grp[mac_idx]), levels = c("LLS", "LHS"))
grp_mal <- factor(as.character(grp[mal_idx]), levels = c("LLS", "LHS"))

score_diff <- function(gm, gv) {
  s <- tapply(v_mif[mac_idx], gm, mean) * tapply(v_rec[mal_idx], gv, mean)
  s["LHS"] - s["LLS"]
}

obs <- score_diff(grp_mac, grp_mal)
cat("Observed LHS - LLS interaction score:", round(unname(obs), 4), "\n")

n_perm <- 10000
nulls <- replicate(n_perm, score_diff(sample(grp_mac), sample(grp_mal)))
p_two <- mean(abs(nulls) >= abs(as.numeric(obs)))
cat("Permutation p (|null| >= |obs|, n=10000):", signif(p_two, 4), "\n")
cat("Null mean:", round(mean(nulls), 4), " sd:", round(sd(nulls), 4), "\n")

res <- data.frame(test = "MIF->CD74/CXCR4 LHS vs LLS",
                  obs_diff = as.numeric(obs),
                  perm_p = p_two,
                  null_mean = mean(nulls), null_sd = sd(nulls))
write.csv(res, "results/mif_permutation.csv", row.names = FALSE)
cat("Saved results/mif_permutation.csv\n")
