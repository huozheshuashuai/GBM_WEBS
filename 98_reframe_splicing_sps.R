# =============================================================================
# 98. WEBS-free splicing analysis: PSI ~ spliceosome mRNA score
#     TCGA-GBM (Firehose RSEM) + TCGA-LGG (Xena TPM) cross-cohort replication
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

splice_score_rsem <- function(mat_genes, mat) {
  g <- intersect(splice, mat_genes)
  colMeans(mat[mat_genes[mat_genes %in% g] |> (\(x) x)() %in% g, , drop = FALSE])
}

# ---------- TCGA-GBM ----------
rsem <- read.delim(gzfile("E:/KIMI/data/Human__TCGA_GBM__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"),
                   check.names = FALSE)
g <- rsem[[1]]; m <- as.matrix(sapply(rsem[, -1], as.numeric)); rownames(m) <- g
m <- m[!is.na(rownames(m)) & rownames(m) != "?", ]
sp_g <- intersect(splice, rownames(m))
sps <- colMeans(m[sp_g, , drop = FALSE])
gbm_score <- data.frame(sample = gsub("\\.", "-", substr(names(sps), 1, 12)), SPS = as.numeric(sps))

psi_raw <- read.delim("data/raw/TCGASpliceSeq/PSI_download_GBM.txt",
                      check.names = FALSE, stringsAsFactors = FALSE)
ev_types <- c("ES","AA","AD","RI","AP","AT","ME")
psi_ev <- psi_raw[psi_raw$splice_type %in% ev_types, ]
event_id <- paste(psi_ev$symbol, psi_ev$as_id, psi_ev$splice_type, sep = "|")
meta_cols <- c("symbol","as_id","splice_type","exons","from_exon","to_exon",
               "novel_splice","pct_with_values","psi_range","std_psi")
keep_idx <- which(!colnames(psi_ev) %in% meta_cols & colnames(psi_ev) != "")
psi_m <- as.matrix(sapply(psi_ev[, keep_idx], as.numeric))
rownames(psi_m) <- event_id
colnames(psi_m) <- gsub("_", "-", colnames(psi_ev)[keep_idx])

common <- intersect(colnames(psi_m), gbm_score$sample)
cat("GBM PSI events:", nrow(psi_m), "| matched samples:", length(common), "\n")
psi_g <- psi_m[, common]
x <- gbm_score$SPS[match(common, gbm_score$sample)]

cor_one <- function(v, x) {
  ok <- !is.na(v) & !is.na(x)
  if (sum(ok) < 30 || sd(v[ok]) == 0) return(c(NA, NA))
  ct <- suppressWarnings(cor.test(v[ok], x[ok], method = "spearman"))
  c(unname(ct$estimate), ct$p.value)
}
res_g <- t(apply(psi_g, 1, cor_one, x))
gbm_res <- data.frame(event = rownames(psi_g), rho = res_g[, 1], p = res_g[, 2])
gbm_res$q <- p.adjust(gbm_res$p, "BH")
gbm_res$symbol <- sub("\\|.*", "", gbm_res$event)
gbm_res$type <- sub(".*\\|", "", gbm_res$event)
sig <- !is.na(gbm_res$q) & gbm_res$q < 0.05 & abs(gbm_res$rho) > 0.2
cat("GBM: significant events (q<0.05 & |rho|>0.2):", sum(sig), "/", sum(!is.na(gbm_res$q)), "\n")
write.csv(gbm_res, "results/reframe/splicing_sps_gbm.csv", row.names = FALSE)

# ---------- TCGA-LGG ----------
tpm <- read.delim(gzfile("E:/KIMI/data/Xena_GDC/TCGA-LGG.star_tpm.tsv.gz"), check.names = FALSE)
lg <- tpm[[1]]
lm_ <- as.matrix(sapply(tpm[, -1], as.numeric)); rownames(lm_) <- lg
lm_ <- lm_[!is.na(rownames(lm_)), ]
if (grepl("^ENSG", rownames(lm_)[1])) {
  suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
  gm2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
          keys = sub("\\..*", "", rownames(lm_)), keytype = "ENSEMBL", columns = "SYMBOL"))
  sym <- gm2$SYMBOL[match(sub("\\..*", "", rownames(lm_)), gm2$ENSEMBL)]
  keep <- !is.na(sym)
  lm_ <- lm_[keep, ]; rownames(lm_) <- sym[keep]
  lm_ <- rowsum(lm_, rownames(lm_))
}
sp_gl <- intersect(splice, rownames(lm_))
lcpm <- log2(lm_ + 1)
sps_l <- colMeans(lcpm[sp_gl, , drop = FALSE])
lgg_score <- data.frame(sample = gsub("\\.", "-", substr(names(sps_l), 1, 12)), SPS = as.numeric(sps_l))

psi_l_list <- list()
for (st in ev_types) {
  f <- file.path("data/raw/TCGASpliceSeq", paste0("LGG_", st), "PSI_download_LGG.txt")
  if (file.exists(f)) psi_l_list[[st]] <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
}
psi_l <- do.call(rbind, psi_l_list)
ev_id_l <- paste(psi_l$symbol, psi_l$as_id, psi_l$splice_type, sep = "|")
keep_l <- which(!colnames(psi_l) %in% meta_cols & colnames(psi_l) != "")
psi_lm <- as.matrix(sapply(psi_l[, keep_l], as.numeric))
rownames(psi_lm) <- ev_id_l
colnames(psi_lm) <- gsub("_", "-", colnames(psi_l)[keep_l])

common_l <- intersect(colnames(psi_lm), lgg_score$sample)
cat("LGG PSI events:", nrow(psi_lm), "| matched samples:", length(common_l), "\n")
psi_lm <- psi_lm[, common_l]
xl <- lgg_score$SPS[match(common_l, lgg_score$sample)]
res_l <- t(apply(psi_lm, 1, cor_one, xl))
lgg_res <- data.frame(event = rownames(psi_lm), rho = res_l[, 1], p = res_l[, 2])
lgg_res$q <- p.adjust(lgg_res$p, "BH")
write.csv(lgg_res, "results/reframe/splicing_sps_lgg.csv", row.names = FALSE)

# ---------- cross-cohort replication ----------
mg <- merge(gbm_res, lgg_res, by = "event", suffixes = c("_gbm", "_lgg"))
mg <- mg[!is.na(mg$rho_gbm) & !is.na(mg$rho_lgg), ]
sig_gbm <- !is.na(mg$q_gbm) & mg$q_gbm < 0.05 & abs(mg$rho_gbm) > 0.2
rep <- mg[sig_gbm, ]
concord <- mean(sign(rep$rho_gbm) == sign(rep$rho_lgg))
ct <- cor.test(rep$rho_gbm, rep$rho_lgg, method = "spearman")
cat("\nreplication: testable GBM-sig events in LGG:", nrow(rep),
    "| sign-concordant:", round(concord * 100, 1), "%",
    "| effect rho =", round(ct$estimate, 2), " p =", signif(ct$p.value, 3), "\n")

sink("results/reframe/splicing_sps_summary.txt", split = TRUE)
cat("=== PSI ~ spliceosome mRNA score (WEBS-free) ===\n")
cat("GBM: events", nrow(gbm_res), "samples", length(common),
    " significant:", sum(sig), "\n")
cat("LGG: events", nrow(lgg_res), "samples", length(common_l), "\n")
cat("replication:", nrow(rep), "events,", round(concord * 100, 1), "% concordant,",
    "rho =", round(ct$estimate, 2), "\n")
cat("\ntop events by |rho| in GBM:\n")
print(head(gbm_res[order(-abs(gbm_res$rho)), c("event", "rho", "q")], 15), row.names = FALSE)
sink()

# figure
p <- ggplot(rep, aes(rho_gbm, rho_lgg)) +
  geom_point(alpha = 0.35, size = 0.8, color = "#7B4EA3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  theme_bw(base_size = 12) +
  labs(x = "rho (PSI vs spliceosome score) TCGA-GBM",
       y = "rho (PSI vs spliceosome score) TCGA-LGG",
       title = sprintf("%d GBM-identified events: %.0f%% sign-concordant in LGG",
                       nrow(rep), concord * 100))
ggsave("figures_reframe/splicing_sps_replication.pdf", p, width = 5.5, height = 5.5)
cat("Done.\n")
