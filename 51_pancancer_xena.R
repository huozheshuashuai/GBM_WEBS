# =============================================================================
# 51. Pan-cancer analysis (33 TCGA cohorts, Xena GDC hub TPM)
#     WEBS x spliceosome correlation + six-gene signature survival per cancer
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(survival); library(ggplot2); library(dplyr)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

XD <- "E:/KIMI/data/Xena_GDC"
WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)
source("08_cptac_proteomics.R")
SPLICE <- get_spliceosome_genes()

# Ensembl (no version) -> symbol map
ens_all <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
sym_map <- mapIds(org.Hs.eg.db, keys = ens_all, column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")

files <- list.files(XD, pattern = "star_tpm.tsv.gz$")
cancers <- sub("TCGA-(.*)\\.star_tpm\\.tsv\\.gz", "\\1", files)

results <- list()
for (ca in cancers) {
  cat("Processing", ca, "...\n")
  tryCatch({
    tpm <- read.delim(gzfile(file.path(XD, paste0("TCGA-", ca, ".star_tpm.tsv.gz"))),
                      check.names = FALSE)
    ens <- sub("\\..*$", "", tpm[[1]])
    sym <- unname(sym_map[ens])
    keep <- !is.na(sym)
    m <- as.matrix(tpm[keep, -1])
    rownames(m) <- sym[keep]
    # dedup by max mean
    o <- order(rowMeans(m), decreasing = TRUE)
    m <- m[o, ]
    m <- m[!duplicated(rownames(m)), ]
    # primary tumors only (sample type 01); LAML/hematologic keep all
    st <- substr(colnames(m), 14, 15)
    if (!(ca %in% c("LAML", "DLBC", "THYM"))) m <- m[, st == "01", drop = FALSE]
    if (ncol(m) < 30) { cat("  too few samples, skip\n"); next }

    w_in <- intersect(WRITERS, rownames(m)); e_in <- intersect(ERASERS, rownames(m))
    webs <- colMeans(m[w_in, , drop = FALSE]) - colMeans(m[e_in, , drop = FALSE])
    sp_in <- intersect(SPLICE, rownames(m))
    splice_score <- colMeans(m[sp_in, , drop = FALSE])
    g6 <- intersect(names(PAPER_COEFS), rownames(m))
    risk <- as.vector(t(m[g6, , drop = FALSE]) %*% PAPER_COEFS[g6])

    ct <- suppressWarnings(cor.test(webs, splice_score, method = "spearman"))

    # survival
    svf <- file.path(XD, paste0("TCGA-", ca, ".survival.tsv.gz"))
    HR <- NA; p_sig <- NA
    if (file.exists(svf)) {
      sv <- read.delim(gzfile(svf))
      sv$sample <- substr(sv$sample, 1, 16)
      df <- data.frame(sample = substr(colnames(m), 1, 16), risk = as.numeric(risk))
      df <- merge(df, sv[, c("sample", "OS.time", "OS")], by = "sample")
      df <- df[!is.na(df$OS.time) & df$OS.time > 0 & !is.na(df$OS), ]
      if (nrow(df) >= 30 && sum(df$OS) >= 5) {
        fit <- coxph(Surv(OS.time, OS) ~ risk, data = df)
        co <- summary(fit)$coefficients
        HR <- unname(co[1, "exp(coef)"]); p_sig <- unname(co[1, "Pr(>|z|)"])
      }
    }
    results[[ca]] <- data.frame(cancer = ca, n = ncol(m),
                                rho_webs_splice = unname(ct$estimate),
                                p_webs_splice = ct$p.value,
                                n_splice_genes = length(sp_in),
                                sig_HR = HR, sig_p = p_sig)
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
res <- do.call(rbind, results)
res$q_webs_splice <- p.adjust(res$p_webs_splice, method = "BH")
write.csv(res, "results/pancancer_33_summary.csv", row.names = FALSE)

sink("results/pancancer_33_summary.txt", split = TRUE)
cat("=== Pan-cancer WEBS x spliceosome correlation (33 TCGA cohorts) ===\n\n")
print(res %>% arrange(desc(rho_webs_splice)))
cat("\nCohorts with q<0.05:", sum(res$q_webs_splice < 0.05, na.rm = TRUE), "/",
    nrow(res), "\n")
cat("GBM rank by rho:", which(res$cancer[order(-res$rho_webs_splice)] == "GBM"), "\n")
sink()

# forest-style plot: WEBS-spliceosome correlation
res$cancer <- factor(res$cancer, levels = res$cancer[order(res$rho_webs_splice)])
res$sig_flag <- ifelse(res$q_webs_splice < 0.05, "q<0.05", "n.s.")
p1 <- ggplot(res, aes(rho_webs_splice, cancer, color = sig_flag)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("q<0.05" = "#7B4EA3", "n.s." = "grey60")) +
  theme_bw(base_size = 11) +
  labs(title = "Pan-cancer: WEBS vs spliceosome pathway correlation",
       x = "Spearman rho", y = NULL, color = NULL)
ggsave("figures_final/fig_pancancer_webs_splice.pdf", p1, width = 6, height = 8)

# signature HR plot
res2 <- res[!is.na(res$sig_HR), ]
res2$cancer <- factor(res2$cancer, levels = res2$cancer[order(res2$sig_HR)])
res2$sig_flag <- ifelse(res2$sig_p < 0.05, "p<0.05", "n.s.")
p2 <- ggplot(res2, aes(sig_HR, cancer, color = sig_flag)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("p<0.05" = "#d62728", "n.s." = "grey60")) +
  theme_bw(base_size = 11) +
  labs(title = "Pan-cancer: six-gene signature Cox HR",
       x = "Hazard ratio", y = NULL, color = NULL)
ggsave("figures_final/fig_pancancer_signature_hr.pdf", p2, width = 6, height = 8)
cat("Done.\n")
