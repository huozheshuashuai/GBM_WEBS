# Merge CopyKAT shard results and produce purity-control summary
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(Seurat); library(dplyr) })

shards <- list.files("data/processed", pattern = "copykat_shard[0-9]+\\.rds$", full.names = TRUE)
pred <- do.call(rbind, lapply(shards, readRDS))
# copykat renames cells: "<sam>.<orig cell name with - -> .>" -> restore
rn <- rownames(pred)
rn <- sub("^[^.]+\\.", "", rn)                 # strip "<sam>." prefix
rn <- sub("\\.([0-9]+)$", "-\\1", rn)          # trailing .1 -> -1
pred$cell <- rn
saveRDS(pred, "data/processed/copykat_calls.rds")

obj <- readRDS("data/processed/GSE138794_final.rds")
obj$copykat <- pred$copykat.pred[match(colnames(obj), pred$cell)]
obj$copykat[obj$copykat == "not.defined"] <- NA

aneu <- subset(obj, subset = copykat == "aneuploid")
dipl <- subset(obj, subset = copykat == "diploid")

sink("results/copykat_summary.txt")
cat("=== CopyKAT purity control ===\n\n")
cat("Called cells:", sum(!is.na(obj$copykat)), "/", ncol(obj),
    " across", length(unique(pred$sample)), "samples\n")
print(table(obj$copykat, useNA = "ifany"))
cat("\nAneuploid (tumor) cells:", ncol(aneu), "\n")
cat("WEBS groups among aneuploid cells:\n"); print(table(aneu$WEBS_group))
cat("\nConcordance: aneuploid vs Neoplastic_cell annotation:\n")
print(table(aneu$original_type, useNA = "ifany"))
kt <- kruskal.test(WEBS ~ WEBS_group, data = aneu@meta.data)
cat("\nWEBS group differences within aneuploid-only cells: KW p =",
    signif(kt$p.value, 3), "\n")
cat("(Acceptance: p < 0.001 means stratification is not a normal-cell artifact)\n")
sink()
cat(readLines("results/copykat_summary.txt"), sep = "\n")
