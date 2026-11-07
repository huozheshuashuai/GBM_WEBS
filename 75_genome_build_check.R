setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(EnsDb.Hsapiens.v75); library(EnsDb.Hsapiens.v86) })

lines <- readLines(gzfile("data/raw/GSE138794/extracted/GSM4119513_SF11964_snATAC_peaks.bed.gz"))
m <- regmatches(lines, regexec("^(chr[0-9XY]+):(\\d+)-(\\d+)", lines))
m <- m[lengths(m) > 0]
peaks <- data.frame(chr = sapply(m, `[`, 2),
                    start = as.numeric(sapply(m, `[`, 3)),
                    end = as.numeric(sapply(m, `[`, 4)))
peaks <- peaks[peaks$chr %in% paste0("chr", 1:22), ]
cat("peaks parsed:", nrow(peaks), "\n")

tss_overlap <- function(edb, label) {
  gn <- genes(edb)
  tss <- data.frame(chr = paste0("chr", seqnames(gn)), pos = start(gn))
  tss <- tss[tss$chr %in% paste0("chr", 1:22), ]
  tss$lo <- tss$pos - 2000; tss$hi <- tss$pos + 2000
  set.seed(1)
  pk <- peaks[sample(nrow(peaks), min(5000, nrow(peaks))), ]
  hit <- 0
  for (c in unique(pk$chr)) {
    pc <- pk[pk$chr == c, ]; tc <- tss[tss$chr == c, ]
    for (i in seq_len(nrow(pc))) {
      if (any(tc$lo <= pc$end[i] & tc$hi >= pc$start[i])) hit <- hit + 1
    }
  }
  cat(label, "TSS overlap:", hit, "/", nrow(pk), "=", round(100 * hit / nrow(pk), 1), "%\n")
}
tss_overlap(EnsDb.Hsapiens.v75, "hg19 (v75)")
tss_overlap(EnsDb.Hsapiens.v86, "hg38 (v86)")
