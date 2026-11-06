# =============================================================================
# 60. Cross-cohort replication of WEBS-associated splicing events (TCGA-LGG)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(org.Hs.eg.db); library(AnnotationDbi) })

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
meta_cols <- c("symbol","as_id","splice_type","exons","from_exon","to_exon",
               "novel_splice","pct_with_values","psi_range","std_psi")

# --- LGG expression (Xena GDC hub) -> WEBS ---
tpm <- read.delim(gzfile("E:/KIMI/data/Xena_GDC/TCGA-LGG.star_tpm.tsv.gz"), check.names = FALSE)
ens <- sub("\\..*$", "", tpm[[1]])
ens_map <- mapIds(org.Hs.eg.db, keys = unique(ens), column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")
sym <- unname(ens_map[ens])
keep <- !is.na(sym)
m <- as.matrix(tpm[keep, -1]); rownames(m) <- sym[keep]
o <- order(rowMeans(m), decreasing = TRUE)
m <- m[o, ]; m <- m[!duplicated(rownames(m)), ]
lgg_webs <- colMeans(m[intersect(WRITERS, rownames(m)), , drop = FALSE]) -
            colMeans(m[intersect(ERASERS, rownames(m)), , drop = FALSE])

# --- LGG PSI (7 event types) ---
psi_list <- list()
for (st in c("ES","RI","AP","AT","AA","AD","ME")) {
  f <- file.path("data/raw/TCGASpliceSeq", paste0("LGG_", st), "PSI_download_LGG.txt")
  d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$splice_type %in% st, ]
  ki <- which(!colnames(d) %in% meta_cols & colnames(d) != "")
  mm <- as.matrix(sapply(d[, ki], as.numeric))
  colnames(mm) <- gsub("_", "-", colnames(d)[ki])
  rownames(mm) <- paste(d$symbol, d$as_id, d$splice_type, sep = "|")
  psi_list[[st]] <- mm
}
# align events
all_ev <- unique(unlist(lapply(psi_list, rownames)))

# --- GBM significant events ---
gbm <- read.csv("results/splicing_webs_corr.csv")
gbm_sig <- gbm[gbm$q < 0.05 & abs(gbm$rho) > 0.2, ]

# --- replicate in LGG ---
rep_one <- function(ev) {
  for (st in names(psi_list)) {
    mm <- psi_list[[st]]
    if (ev %in% rownames(mm)) {
      x <- mm[ev, ]
      common <- intersect(names(x), substr(names(lgg_webs), 1, 12))
      common <- intersect(common, colnames(mm))
      xx <- x[common]; ww <- lgg_webs[match(common, substr(names(lgg_webs), 1, 12))]
      ok <- !is.na(xx) & !is.na(ww)
      if (sum(ok) < 30 || sd(xx[ok]) == 0 || sd(ww[ok]) == 0) return(c(NA, NA, 0))
      ct <- suppressWarnings(cor.test(xx[ok], ww[ok], method = "spearman"))
      return(c(unname(ct$estimate), ct$p.value, sum(ok)))
    }
  }
  c(NA, NA, 0)
}

cat("Testing", nrow(gbm_sig), "GBM events in LGG...\n")
rep <- t(sapply(gbm_sig$event, rep_one))
out <- data.frame(event = gbm_sig$event, rho_gbm = gbm_sig$rho,
                  rho_lgg = rep[, 1], p_lgg = rep[, 2], n_lgg = rep[, 3])
out <- out[!is.na(out$rho_lgg), ]
out$q_lgg <- p.adjust(out$p_lgg, method = "BH")
out$concordant <- sign(out$rho_gbm) == sign(out$rho_lgg)
write.csv(out, "results/splicing_lgg_replication.csv", row.names = FALSE)

sink("results/splicing_lgg_summary.txt", split = TRUE)
cat("=== Cross-cohort splicing replication (GBM -> LGG) ===\n\n")
cat("GBM WEBS-associated events:", nrow(gbm_sig), "\n")
cat("Tested in LGG:", nrow(out), "\n")
cat("Sign concordant:", sum(out$concordant), "/", nrow(out),
    sprintf("(%.0f%%)\n", 100 * mean(out$concordant)))
bt <- binom.test(sum(out$concordant), nrow(out), p = 0.5)
cat("Binomial test vs 50%: p =", signif(bt$p.value, 3), "\n")
sig_lgg <- out[out$q_lgg < 0.05, ]
cat("Also significant in LGG (q<0.05):", nrow(sig_lgg),
    sprintf("(%.0f%%), of which concordant: %d/%d\n",
            100 * nrow(sig_lgg) / nrow(out), sum(sig_lgg$concordant), nrow(sig_lgg)))
ct <- cor.test(out$rho_gbm, out$rho_lgg, method = "spearman")
cat(sprintf("Effect-size correlation GBM vs LGG: rho = %.3f, p = %.3g\n",
            ct$estimate, ct$p.value))
sink()
cat(readLines("results/splicing_lgg_summary.txt"), sep = "\n")

p <- ggplot(out, aes(rho_gbm, rho_lgg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.3, size = 0.8, color = "#7B4EA3") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(title = "Cross-cohort replication of WEBS-associated splicing events",
       subtitle = paste0(sum(out$concordant), "/", nrow(out), " sign-concordant (",
                         round(100 * mean(out$concordant)), "%), effect rho = ",
                         round(ct$estimate, 2)),
       x = "rho (PSI vs WEBS) in TCGA-GBM", y = "rho (PSI vs WEBS) in TCGA-LGG")
ggsave("figures_final/fig_splicing_lgg_replication.pdf", p, width = 7, height = 6)
cat("Done.\n")
