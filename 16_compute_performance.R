# =============================================================================
# 16. Computational Performance Benchmark
# =============================================================================
# Runtime & memory of WEBS (UCell) vs ssGSEA (GSVA) vs AUCell vs all-gene mean
# at increasing cell counts (1k / 5k / 10k / 25k / 50k).
# Input:  data/processed/GSE138794_webs_annotated.rds
# Output: results/computational_performance.csv + console table for Methods
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(Seurat)
  library(UCell)
  library(AUCell)
  library(GSVA)
  library(dplyr)
})

WRITER_GENES <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASER_GENES <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")

obj <- readRDS("data/processed/GSE138794_webs_annotated.rds")
expr_full <- GetAssayData(obj, layer = "data")

bench_one <- function(mat, method) {
  gc()
  mem_before <- sum(gc()[, 2])  # Mb used
  t <- system.time({
    score <- switch(method,
      WEBS = {
        w <- ScoreSignatures_UCell(mat, features = list(W = WRITER_GENES))
        e <- ScoreSignatures_UCell(mat, features = list(E = ERASER_GENES))
        w[1, ] - e[1, ]
      },
      AUCell = {
        rk <- AUCell_buildRankings(mat, plotStats = FALSE)
        wa <- AUCell_calcAUC(list(W = WRITER_GENES), rk)
        ea <- AUCell_calcAUC(list(E = ERASER_GENES), rk)
        getAUC(wa)["W", ] - getAUC(ea)["E", ]
      },
      ssGSEA = {
        m <- as.matrix(mat)
        w <- gsva(ssgseaParam(m, list(W = WRITER_GENES)), verbose = FALSE)
        e <- gsva(ssgseaParam(m, list(E = ERASER_GENES)), verbose = FALSE)
        w[1, ] - e[1, ]
      },
      MeanExpr = {
        Matrix::colMeans(mat[intersect(WRITER_GENES, rownames(mat)), , drop = FALSE]) -
          Matrix::colMeans(mat[intersect(ERASER_GENES, rownames(mat)), , drop = FALSE])
      }
    )
  })
  mem_after <- sum(gc()[, 2])
  data.frame(elapsed_sec = t["elapsed"], delta_mem_mb = max(0, mem_after - mem_before))
}

sizes <- c(1000, 5000, 10000, 25000, min(50000, ncol(expr_full)))
methods <- c("WEBS", "AUCell", "ssGSEA", "MeanExpr")

results <- list()
for (n in sizes) {
  cells <- sample(colnames(expr_full), n)
  mat <- expr_full[, cells]
  for (m in methods) {
    message("n=", n, " method=", m)
    r <- tryCatch(bench_one(mat, m), error = function(e) {
      message("  FAILED: ", e$message)
      data.frame(elapsed_sec = NA, delta_mem_mb = NA)
    })
    results[[length(results) + 1]] <- data.frame(n_cells = n, method = m, r)
  }
}

perf <- do.call(rbind, results)
dir.create("results", showWarnings = FALSE)
write.csv(perf, "results/computational_performance.csv", row.names = FALSE)
print(perf)

# Methods-ready summary text
sink("results/computational_performance_summary.txt")
cat("Computational performance (WEBS vs alternatives)\n")
cat("Machine:", Sys.info()["machine"], "| RAM:", "16 GB", "\n\n")
print(perf %>% arrange(method, n_cells))
cat("\nWEBS at 10,000 cells: ",
    round(perf$elapsed_sec[perf$method == "WEBS" & perf$n_cells == 10000], 1), " s\n")
sink()

library(ggplot2)
p <- ggplot(perf, aes(x = n_cells, y = elapsed_sec, color = method)) +
  geom_line() + geom_point(size = 2) +
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
  theme_bw(base_size = 12) +
  labs(x = "Number of cells (log scale)", y = "Runtime, seconds (log scale)",
       color = "Method", title = "Computational scalability of WEBS and alternatives")
dir.create("figures", showWarnings = FALSE)
pdf("figures/figS_performance.pdf", width = 7, height = 5)
print(p)
dev.off()

message("Done. See results/computational_performance.csv and figures/figS_performance.pdf")
