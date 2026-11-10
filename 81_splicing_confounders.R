# =============================================================================
# 81. Splicing confounder control: PSI ~ WEBS with multi-covariate adjustment
#     (immune score + proliferation + MES + hypoxia programs)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ppcor) })

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")

rsem <- read.delim(gzfile("E:/KIMI/data/Human__TCGA_GBM__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"),
                   check.names = FALSE)
m <- as.matrix(sapply(rsem[, -1], as.numeric)); rownames(m) <- rsem[[1]]
m <- m[!is.na(rownames(m)) & rownames(m) != "?", ]
colnames(m) <- gsub("\\.", "-", substr(colnames(m), 1, 12))
webs <- colMeans(m[intersect(WRITERS, rownames(m)), , drop = FALSE]) -
        colMeans(m[intersect(ERASERS, rownames(m)), , drop = FALSE])
score_of <- function(gs) colMeans(m[intersect(gs, rownames(m)), , drop = FALSE])
covs <- data.frame(
  immune = score_of(c("PTPRC","CD68","CD163","CD14","CD3D","CD3E","CD8A","MS4A1","NCAM1","FOXP3")),
  prolif = score_of(c("MKI67","TOP2A","PCNA","MCM2","MCM6")),
  mes    = score_of(c("VIM","CD44","FN1","SERPINE1","TGFBI","COL1A2")),
  proneural = score_of(c("OLIG2","SOX10","PDGFRA","DLL3","ASCL1")),
  hypox  = score_of(c("CA9","VEGFA","LDHA","SLC2A1","ADM"))
)

psi_raw <- read.delim("data/raw/TCGASpliceSeq/PSI_download_GBM.txt",
                      check.names = FALSE, stringsAsFactors = FALSE)
ev_types <- c("ES","RI","AP","AT","AA","AD","ME")
psi_ev <- psi_raw[psi_raw$splice_type %in% ev_types, ]
meta_cols <- c("symbol","as_id","splice_type","exons","from_exon","to_exon",
               "novel_splice","pct_with_values","psi_range","std_psi")
ki <- which(!colnames(psi_ev) %in% meta_cols & colnames(psi_ev) != "")
psi_m <- as.matrix(sapply(psi_ev[, ki], as.numeric))
colnames(psi_m) <- gsub("_", "-", colnames(psi_ev)[ki])

res_old <- read.csv("results/splicing_webs_corr.csv")
sig_ev <- res_old$event[res_old$q < 0.05 & abs(res_old$rho) > 0.2]
event_id <- paste(psi_ev$symbol, psi_ev$as_id, psi_ev$splice_type, sep = "|")
idx <- match(sig_ev, event_id)
idx <- idx[!is.na(idx)]

common <- Reduce(intersect, list(colnames(psi_m), names(webs), rownames(covs)))
psi_m <- psi_m[, common]
w <- webs[common]; C <- as.matrix(covs[common, ])

sink("results/splicing_confounders.txt", split = TRUE)
cat("=== Splicing confounder control (n =", length(common), "samples) ===\n\n")
cat("Covariates: immune, proliferation, MES, proneural, hypoxia\n\n")

# stepwise: original rho -> +immune -> +all covariates
retain_step <- function(events_idx) {
  n_all <- length(events_idx)
  r_immune <- r_full <- numeric(n_all)
  for (j in seq_along(events_idx)) {
    x <- psi_m[events_idx[j], ]
    ok <- !is.na(x)
    if (sum(ok) < 30) { r_immune[j] <- r_full[j] <- NA; next }
    pc1 <- suppressWarnings(pcor.test(x[ok], w[ok], C[ok, "immune"], method = "spearman"))
    pc2 <- suppressWarnings(pcor.test(x[ok], w[ok], C[ok, ], method = "spearman"))
    r_immune[j] <- pc1$p.value; r_full[j] <- pc2$p.value
  }
  data.frame(q_immune = p.adjust(r_immune, "BH"), q_full = p.adjust(r_full, "BH"))
}
q <- retain_step(idx)
n_ret_immune <- sum(q$q_immune < 0.05, na.rm = TRUE)
n_ret_full <- sum(q$q_full < 0.05, na.rm = TRUE)
cat("WEBS-associated events:", length(idx), "\n")
cat("Retained after immune adjustment:", n_ret_immune,
    sprintf("(%.0f%%)\n", 100 * n_ret_immune / length(idx)))
cat("Retained after ALL five covariates:", n_ret_full,
    sprintf("(%.0f%%)\n", 100 * n_ret_full / length(idx)))
sink()
cat(readLines("results/splicing_confounders.txt"), sep = "\n")
write.csv(data.frame(event = sig_ev[seq_along(idx)], q_immune = q$q_immune, q_full = q$q_full),
          "results/splicing_confounders.csv", row.names = FALSE)
cat("Done.\n")
