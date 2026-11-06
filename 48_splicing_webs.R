# =============================================================================
# 48. Functional consequence of spliceosome accumulation: alternative splicing
#     TCGASpliceSeq PSI (TCGA-GBM) correlated with bulk WEBS
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db) })

WRITERS <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
ERASERS <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")

# --- 1. TCGA-GBM bulk WEBS from Firehose RSEM ---
rsem <- read.delim(gzfile("E:/KIMI/data/Human__TCGA_GBM__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"),
                   check.names = FALSE)
g <- rsem[[1]]
m <- as.matrix(sapply(rsem[, -1], as.numeric))
rownames(m) <- g
m <- m[!is.na(rownames(m)) & rownames(m) != "?", ]
webs <- colMeans(m[intersect(WRITERS, rownames(m)), , drop = FALSE]) -
        colMeans(m[intersect(ERASERS, rownames(m)), , drop = FALSE])
tcga_webs <- data.frame(sample = names(webs), WEBS = as.numeric(webs))
tcga_webs$sample <- substr(tcga_webs$sample, 1, 12)   # TCGA-02-0047
tcga_webs$sample <- gsub("\\.", "-", tcga_webs$sample)

# --- 2. PSI matrix ---
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

common <- intersect(colnames(psi_m), tcga_webs$sample)
cat("PSI events:", nrow(psi_m), " PSI samples:", ncol(psi_m),
    " matched with WEBS:", length(common), "\n")
psi_m <- psi_m[, common]
w <- tcga_webs$WEBS[match(common, tcga_webs$sample)]

# --- 3. Spearman correlation PSI ~ WEBS ---
cor_one <- function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 30 || sd(x[ok]) == 0) return(c(NA, NA))
  ct <- suppressWarnings(cor.test(x[ok], w[ok], method = "spearman"))
  c(unname(ct$estimate), ct$p.value)
}
res <- t(apply(psi_m, 1, cor_one))
res <- data.frame(event = rownames(psi_m), rho = res[, 1], p = res[, 2])
res$q <- p.adjust(res$p, method = "BH")
res$symbol <- psi_ev$symbol
res$splice_type <- psi_ev$splice_type
res <- res[!is.na(res$p), ]
write.csv(res, "results/splicing_webs_corr.csv", row.names = FALSE)

sig <- res[res$q < 0.05 & abs(res$rho) > 0.2, ]
cat("Significant WEBS-associated AS events (q<0.05, |rho|>0.2):", nrow(sig), "\n")
print(table(sig$splice_type, sign(sig$rho)))

# --- 4. summary ---
sink("results/splicing_summary.txt", split = TRUE)
cat("=== Alternative splicing vs WEBS (TCGA-GBM, TCGASpliceSeq) ===\n")
cat("Events tested:", nrow(res), " Samples:", length(common), "\n")
cat("Significant (q<0.05 & |rho|>0.2):", nrow(sig),
    sprintf("(%.1f%%)\n", 100 * nrow(sig) / nrow(res)))
cat("Median |rho| among significant:", round(median(abs(sig$rho)), 3), "\n\n")
cat("Event-type breakdown among significant:\n")
print(table(sig$splice_type))
cat("\nTop 20 events by |rho|:\n")
print(sig %>% arrange(desc(abs(rho))) %>% head(20) %>%
        dplyr::select(event, splice_type, rho, q))
sink()

# --- 5. enrichment of genes with significant events ---
sig_genes <- unique(sig$symbol)
bg_genes <- unique(res$symbol)
ego <- enrichGO(sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                ont = "BP", universe = bg_genes, pAdjustMethod = "BH",
                qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  write.csv(as.data.frame(ego), "results/splicing_enrichment.csv", row.names = FALSE)
  print(head(as.data.frame(ego)[, c("Description","p.adjust","Count")], 15))
}

# --- 6. figures ---
res$sig_flag <- ifelse(res$q < 0.05 & abs(res$rho) > 0.2, "Significant", "NS")
p_vol <- ggplot(res, aes(rho, -log10(p), color = sig_flag)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("Significant" = "#7B4EA3", "NS" = "grey75")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 12) +
  labs(title = paste0("AS events vs WEBS (n=", nrow(res), " events, ",
                      length(common), " samples)"),
       subtitle = paste0(nrow(sig), " significant events (q<0.05, |rho|>0.2)"),
       x = "Spearman rho (PSI vs WEBS)", y = "-log10(p)", color = NULL)
ggsave("figures_final/fig_splicing_volcano.pdf", p_vol, width = 7, height = 5)

tt <- as.data.frame(table(sig$splice_type))
colnames(tt) <- c("type", "count")
p_bar <- ggplot(tt, aes(reorder(type, -count), count, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "WEBS-associated AS events by type", x = "Splice event type",
       y = "Number of significant events")
ggsave("figures_final/fig_splicing_eventtypes.pdf", p_bar, width = 5.5, height = 4.5)

cat("Done.\n")
