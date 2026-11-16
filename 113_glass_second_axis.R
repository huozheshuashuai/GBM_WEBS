# =============================================================================
# 113. Task 3 (V6): GLASS independent longitudinal validation of the
#      second axis (ROS / IL2-STAT5 / EMT) at RNA level
#      GLASS paired primary-recurrent RNA-seq (independent of KNCC)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

tpm <- read.delim("E:/KIMI/data/GLASS_extracted/gene_tpm_matrix_all_samples.tsv",
                  check.names = FALSE)
pairs <- read.csv("E:/KIMI/data/GLASS_extracted/analysis_rnaseq_pairs.csv")
cat("GLASS pairs table cols:", paste(colnames(pairs), collapse = ", "), "\n")
print(head(pairs, 3))

# gene matrix
g <- tpm[[1]]
m <- as.matrix(sapply(tpm[, -1], as.numeric)); rownames(m) <- g
m <- m[!is.na(rownames(m)) & rownames(m) != "", ]
cat("GLASS matrix:", nrow(m), "genes x", ncol(m), "samples\n")

# frozen second-axis modules
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f2 <- strsplit(l, "\t")[[1]]; f2[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(x) intersect(x, rownames(m)))
mods <- list(
  PROLIF  = intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                        "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(m)),
  NEURO   = intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                        "GAP43","TUBB3","MAP2"), rownames(m)),
  ROS     = hm$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
  IL2STAT5 = hm$HALLMARK_IL2_STAT5_SIGNALING,
  HYPOXIA = hm$HALLMARK_HYPOXIA,
  EMT     = hm$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  OXPHOS  = hm$HALLMARK_OXIDATIVE_PHOSPHORYLATION
)
cat("module sizes in GLASS:", paste(sapply(mods, length), collapse = ", "), "\n")

# identify primary/recurrent pair columns from pairs table
cat("\nGLASS pairs total:", nrow(pairs), "\n")

# match barcodes to matrix columns (matrix colnames may contain the barcode as prefix)
match_col <- function(bc) {
  bc2 <- gsub("-", ".", bc, fixed = TRUE)
  hit <- grep(paste0("^", bc2), colnames(m), value = TRUE)
  if (!length(hit)) {
    bc3 <- gsub("-", ".", substr(bc, 1, 19), fixed = TRUE)
    hit <- grep(paste0("^", bc3), colnames(m), value = TRUE)
  }
  if (length(hit) >= 1) hit[1] else NA
}
pairs$col_a <- sapply(pairs$tumor_barcode_a, match_col)
pairs$col_b <- sapply(pairs$tumor_barcode_b, match_col)
okp <- pairs[!is.na(pairs$col_a) & !is.na(pairs$col_b), ]
cat("pairs matched to matrix:", nrow(okp), "\n")

# log-TPM module scores
lmat <- log2(m + 0.5)
score <- t(sapply(mods, function(g) colMeans(lmat[g, , drop = FALSE], na.rm = TRUE)))

res <- lapply(names(mods), function(mn) {
  s1 <- score[mn, okp$col_a]; s2 <- score[mn, okp$col_b]
  d <- s2 - s1
  wt <- suppressWarnings(wilcox.test(d, mu = 0))
  data.frame(module = mn, n_pairs = sum(!is.na(d)),
             mean_delta = mean(d, na.rm = TRUE), sd_delta = sd(d, na.rm = TRUE),
             dir_consistency = max(mean(d > 0, na.rm = TRUE), mean(d < 0, na.rm = TRUE)),
             paired_p = wt$p.value)
})
res <- do.call(rbind, res)
res$q <- p.adjust(res$paired_p, "BH")

sink("results/kncc/glass_second_axis.txt", split = TRUE)
cat("=== GLASS independent longitudinal validation (RNA level) ===\n")
cat("independent of KNCC; pairs:", nrow(okp), "\n\n")
print(res, row.names = FALSE)
cat("\nKNCC protein-level reference: ROS +0.12 (77% up), IL2STAT5 +0.04 (63% up),\n")
cat("EMT -0.04 (bidirectional 56%), HYPOXIA +0.02 (bidirectional 55%),\n")
cat("PROLIF -0.22 (71% down), NEURO +0.52 (74% up)\n")
sink()

# figure
res$module <- factor(res$module, levels = res$module[order(res$mean_delta)])
p <- ggplot(res, aes(module, mean_delta)) +
  geom_col(aes(fill = q < 0.05)) +
  geom_errorbar(aes(ymin = mean_delta - sd_delta / sqrt(n_pairs),
                    ymax = mean_delta + sd_delta / sqrt(n_pairs)), width = 0.3) +
  geom_hline(yintercept = 0) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#7B4EA3", "FALSE" = "grey70")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(x = NULL, y = "Mean paired delta (recurrent - primary, log2 TPM)",
       title = sprintf("GLASS independent longitudinal validation (n = %d pairs, RNA level)",
                       nrow(okp)))
ggsave("figures_kncc/glass_second_axis.pdf", p, width = 7, height = 4.5)
write.csv(res, "results/kncc/glass_second_axis.csv", row.names = FALSE)
cat("Done.\n")
