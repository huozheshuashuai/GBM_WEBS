# =============================================================================
# 56. Bulk splicing: immune-fraction correction (confounder control)
#     Partial correlation PSI ~ WEBS | TME score
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ppcor) })

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
TME_MARKERS <- c("PTPRC","CD68","CD163","CD14","CD3D","CD3E","CD8A","MS4A1","NCAM1","FOXP3")

# TCGA-GBM expression (same source as 48)
rsem <- read.delim(gzfile("E:/KIMI/data/Human__TCGA_GBM__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"),
                   check.names = FALSE)
g <- rsem[[1]]
m <- as.matrix(sapply(rsem[, -1], as.numeric)); rownames(m) <- g
m <- m[!is.na(rownames(m)) & rownames(m) != "?", ]
colnames(m) <- gsub("\\.", "-", substr(colnames(m), 1, 12))
webs <- colMeans(m[intersect(WRITERS, rownames(m)), , drop = FALSE]) -
        colMeans(m[intersect(ERASERS, rownames(m)), , drop = FALSE])
tme <- colMeans(m[intersect(TME_MARKERS, rownames(m)), , drop = FALSE])

# PSI matrix
psi_raw <- read.delim("data/raw/TCGASpliceSeq/PSI_download_GBM.txt",
                      check.names = FALSE, stringsAsFactors = FALSE)
ev_types <- c("ES","AA","AD","RI","AP","AT","ME")
psi_ev <- psi_raw[psi_raw$splice_type %in% ev_types, ]
meta_cols <- c("symbol","as_id","splice_type","exons","from_exon","to_exon",
               "novel_splice","pct_with_values","psi_range","std_psi")
keep_idx <- which(!colnames(psi_ev) %in% meta_cols & colnames(psi_ev) != "")
psi_m <- as.matrix(sapply(psi_ev[, keep_idx], as.numeric))
colnames(psi_m) <- gsub("_", "-", colnames(psi_ev)[keep_idx])
event_id <- paste(psi_ev$symbol, psi_ev$as_id, psi_ev$splice_type, sep = "|")

common <- Reduce(intersect, list(colnames(psi_m), names(webs), names(tme)))
cat("Matched samples:", length(common), "\n")
psi_m <- psi_m[, common]
w <- webs[common]; tm <- tme[common]

# correlate TME with WEBS first
ct0 <- cor.test(w, tm, method = "spearman")
cat("WEBS vs TME score: rho =", round(ct0$estimate, 3), " p =", signif(ct0$p.value, 3), "\n")

# partial correlation per event (subset to significant events from 48 for speed)
res_old <- read.csv("results/splicing_webs_corr.csv")
sig_events <- res_old$event[res_old$q < 0.05 & abs(res_old$rho) > 0.2]
sig_events <- intersect(sig_events, paste(psi_ev$symbol, psi_ev$as_id, psi_ev$splice_type, sep = "|"))
idx <- match(sig_events, event_id)
cat("Re-testing", length(idx), "WEBS-associated events with TME correction\n")

partial_one <- function(i) {
  x <- psi_m[i, ]
  ok <- !is.na(x)
  if (sum(ok) < 30) return(c(NA, NA))
  pc <- suppressWarnings(pcor.test(x[ok], w[ok], tm[ok], method = "spearman"))
  c(unname(pc$estimate), pc$p.value)
}
pres <- t(sapply(idx, partial_one))
out <- data.frame(event = sig_events, rho_partial = pres[, 1], p_partial = pres[, 2])
out$q_partial <- p.adjust(out$p_partial, method = "BH")
out <- merge(out, res_old[, c("event","rho","q")], by = "event")
out$retained <- out$q_partial < 0.05 & abs(out$rho_partial) > 0.2
write.csv(out, "results/splicing_tme_corrected.csv", row.names = FALSE)

sink("results/splicing_tme_summary.txt", split = TRUE)
cat("=== Bulk splicing: immune-fraction correction ===\n")
cat("WEBS vs TME score: rho =", round(ct0$estimate, 3), " p =", signif(ct0$p.value, 3), "\n")
cat("WEBS-associated events:", nrow(out), "\n")
cat("Retained after TME correction:", sum(out$retained, na.rm = TRUE),
    sprintf("(%.0f%%)\n", 100 * mean(out$retained, na.rm = TRUE)))
cat("Median |rho| before:", round(median(abs(out$rho), na.rm = TRUE), 3),
    " after:", round(median(abs(out$rho_partial), na.rm = TRUE), 3), "\n")
sink()
cat(readLines("results/splicing_tme_summary.txt"), sep = "\n")

