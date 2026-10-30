# =============================================================================
# 33. Fig 8: snATAC-seq promoter accessibility at writer/eraser loci
# =============================================================================
# Fragment files not deposited; use peak x cell matrices directly:
# per-sample accessibility at writer/eraser promoter peaks (peaks overlapping
# TSS +/- 2kb), aggregated, compared across samples' scRNA WEBS.
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Matrix); library(dplyr); library(ggplot2); library(patchwork)
})
set.seed(42)
dir.create("figures_v2", showWarnings = FALSE)

# writer/eraser TSS coordinates (GSE138794 ATAC is hg19/GRCh37 -> EnsDb v75)
suppressMessages(library(EnsDb.Hsapiens.v75))
GENES <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP",
           "HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")
gn <- genes(EnsDb.Hsapiens.v75, filter = GeneNameFilter(GENES))
tss <- data.frame(gene = gn$gene_name, seqnames = as.character(seqnames(gn)),
                  tss = ifelse(strand(gn) == "+", start(gn), end(gn)))

raw <- "data/raw/GSE138794/extracted"
atac_samples <- list.files(raw, pattern = "snATAC_matrix.mtx.gz")
atac_samples <- unique(sub("_snATAC_matrix.mtx.gz$", "", atac_samples))

results <- list()
for (s in atac_samples) {
  message("Processing ", s)
  peaks <- read.delim(gzfile(file.path(raw, paste0(s, "_snATAC_peaks.bed.gz"))),
                      header = FALSE, sep = "", stringsAsFactors = FALSE)
  if (ncol(peaks) == 1) {
    # format: "chr1:180949-181688" single string
    parts <- do.call(rbind, strsplit(peaks[[1]], "[:-]"))
    peaks <- data.frame(chr = parts[, 1], start = as.numeric(parts[, 2]),
                        end = as.numeric(parts[, 3]))
  } else {
    colnames(peaks)[1:3] <- c("chr", "start", "end")
  }
  # match TSS chromosomes (peaks may lack "chr" prefix)
  has_chr <- grepl("^chr", peaks$chr[1])
  tss_c <- tss; if (has_chr) tss_c$seqnames <- paste0("chr", tss_c$seqnames)

  # promoter peaks: overlap TSS +/- 2kb
  peak_gene <- list()
  for (i in seq_len(nrow(tss_c))) {
    g <- tss_c$gene[i]; chr <- tss_c$seqnames[i]; t0 <- tss_c$tss[i]
    hit <- which(peaks$chr == chr & peaks$start <= t0 + 2000 & peaks$end >= t0 - 2000)
    if (length(hit) > 0) peak_gene[[g]] <- hit
  }

  mat <- readMM(gzfile(file.path(raw, paste0(s, "_snATAC_matrix.mtx.gz"))))
  # orientation: peaks x cells (standard 10x ATAC)
  if (nrow(mat) != nrow(peaks)) mat <- t(mat)

  per_gene <- sapply(names(peak_gene), function(g) {
    colSums(mat[peak_gene[[g]], , drop = FALSE])  # per cell total accessibility
  })
  if (length(peak_gene) == 0) { message("  no promoter peaks, skipping"); next }
  df <- data.frame(sample = s,
                   gene = rep(names(peak_gene), each = ncol(mat)),
                   access = as.vector(per_gene))
  results[[s]] <- df
}

all_df <- do.call(rbind, results)
all_df$type <- ifelse(all_df$gene %in% c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP"),
                      "Writer", "Eraser")
write.csv(all_df, "results/snatac_promoter_accessibility.csv", row.names = FALSE)

p <- ggplot(all_df, aes(x = gene, y = log1p(access), fill = type)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c("Writer" = "#d62728", "Eraser" = "#1f77b4")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "snATAC-seq promoter accessibility (TSS +/- 2kb) at writer/eraser loci",
       x = NULL, y = "log1p(peak counts)", fill = NULL)

pdf("figures_v2/Fig8_snATAC_promoter.pdf", width = 10, height = 5.5)
print(p)
dev.off()
message("Fig8 done")

# ---- Panel B: per-sample mean accessibility vs patient-level scRNA WEBS ----
obj <- readRDS("data/processed/GSE138794_final.rds")
pt_webs <- tapply(obj$WEBS, obj$sample, median)
# snATAC sample names like GSM4119513_SF11964 -> patient SF11964
summ <- all_df %>% mutate(patient = sub("^GSM[0-9]+_", "", sample)) %>%
  group_by(patient, gene, type) %>%
  summarise(mean_access = mean(log1p(access)), .groups = "drop") %>%
  mutate(WEBS = pt_webs[patient])
summ <- summ[!is.na(summ$WEBS), ]

cors <- summ %>% group_by(gene) %>%
  summarise(r = tryCatch(cor(WEBS, mean_access, method = "spearman"), error = function(e) NA_real_),
            p = tryCatch(cor.test(WEBS, mean_access, method = "spearman")$p.value, error = function(e) NA_real_),
            type = dplyr::first(type), .groups = "drop")
write.csv(cors, "results/snatac_webs_correlation.csv", row.names = FALSE)

pB <- ggplot(summ, aes(x = WEBS, y = mean_access, color = type)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, aes(group = gene), alpha = 0.3, linewidth = 0.4) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Writer" = "#d62728", "Eraser" = "#1f77b4")) +
  theme_bw(base_size = 10) +
  labs(title = "B. Promoter accessibility vs patient WEBS (6 snATAC samples)",
       x = "Patient median WEBS (scRNA-seq)", y = "Mean promoter accessibility (log1p)",
       color = NULL)

pdf("figures_v2/Fig8_snATAC_promoter.pdf", width = 13, height = 8)
print(p / pB)
dev.off()
message("Fig8 v2 done")
