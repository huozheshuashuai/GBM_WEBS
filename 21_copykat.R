# =============================================================================
# 21. CopyKAT purity control - parallel shard version
# Usage: Rscript 21_copykat.R <shard_index> <n_shards>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
shard <- as.integer(args[1]); n_shards <- as.integer(args[2])

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(copykat); library(dplyr) })
set.seed(42)

obj <- readRDS("data/processed/GSE138794_final.rds")

samples <- unique(obj$sample)
samples <- samples[seq(shard, length(samples), by = n_shards)]

calls <- list()
for (smp in samples) {
  cells <- WhichCells(obj, expression = sample == smp)
  if (length(cells) < 100) next
  # subsample for speed (CopyKAT is per-sample robust at 1000+ cells)
  if (length(cells) > 1200) cells <- sample(cells, 1200)
  mat <- as.matrix(GetAssayData(obj, layer = "counts")[, cells])
  message("CopyKAT on ", smp, " (", length(cells), " cells)")
  ck <- tryCatch(
    copykat(rawmat = mat, id.type = "S", ngene.chr = 5,
            win.size = 25, KS.cut = 0.1, sam.name = smp,
            distance = "euclidean", n.cores = 1, output.seg = "FALSE",
            plot.genes = "FALSE", genome = "hg20"),
    error = function(e) { message("  failed: ", e$message); NULL })
  if (!is.null(ck)) {
    pr <- ck$prediction
    pr$sample <- smp
    calls[[smp]] <- pr
  }
}

if (length(calls) > 0) {
  pred <- do.call(rbind, calls)
  saveRDS(pred, paste0("data/processed/copykat_shard", shard, ".rds"))
}
message("Shard ", shard, " done: ", length(calls), " samples")
