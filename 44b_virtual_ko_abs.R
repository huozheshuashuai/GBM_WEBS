# =============================================================================
# 44b. Virtual KO redo: absolute connectivity-based sensitivity
#      (GENIE3 weights are unsigned -> use |shift|, not signed shift)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
writers <- c("EP300","KAT2A","KAT2B","KAT5","KAT8","CREBBP")
erasers <- c("HDAC1","HDAC2","HDAC3","SIRT1","SIRT2","SIRT3","SIRT6","SIRT7")

wm <- readRDS("results/genie3_grn.rds")   # regulators x targets
adj_full <- t(as.matrix(wm))
adj_full <- adj_full / (rowSums(abs(adj_full)) + 1e-9)
regs <- colnames(adj_full)

predict_ko <- function(ko_gene) {
  v <- rep(0, length(regs)); names(v) <- regs
  v[ko_gene] <- -1
  step1 <- as.vector(adj_full %*% v)
  v2 <- rep(0, length(regs)); names(v2) <- regs
  shared <- intersect(rownames(adj_full), regs)
  v2[shared] <- step1[match(shared, rownames(adj_full))]
  step2 <- as.vector(adj_full %*% v2)
  out <- step1 + 0.5 * step2
  names(out) <- rownames(adj_full)
  out
}

ko_df <- do.call(rbind, lapply(c(writers, erasers), function(g) {
  if (!g %in% rownames(wm)) return(NULL)
  pred <- predict_ko(g)
  data.frame(ko = g, gene = names(pred), abs_shift = abs(pred))
}))
ko_df$is_spliceosome <- ko_df$gene %in% splice

summ <- ko_df %>%
  group_by(ko) %>%
  summarise(
    spliceosome_abs = mean(abs_shift[is_spliceosome], na.rm = TRUE),
    other_abs = mean(abs_shift[!is_spliceosome], na.rm = TRUE),
    n_splice = sum(is_spliceosome), .groups = "drop")
summ$ratio <- summ$spliceosome_abs / summ$other_abs
summ$type <- ifelse(summ$ko %in% writers, "Writer", "Eraser")
write.csv(summ, "results/virtual_ko_abs_summary.csv", row.names = FALSE)

sink("results/virtual_ko_abs.txt", split = TRUE)
cat("=== Virtual KO (absolute connectivity-based sensitivity) ===\n\n")
cat("S_r = mean |predicted shift| per KO; ratio = spliceosome / other genes\n\n")
print(summ %>% arrange(desc(ratio)))
wt <- wilcox.test(ko_df$abs_shift[ko_df$is_spliceosome],
                  ko_df$abs_shift[!ko_df$is_spliceosome])
cat(sprintf("\nAll KOs pooled: spliceosome |shift| vs other: median %.5f vs %.5f, Wilcoxon p = %.3g\n",
            median(ko_df$abs_shift[ko_df$is_spliceosome]),
            median(ko_df$abs_shift[!ko_df$is_spliceosome]), wt$p.value))
sink()
cat(readLines("results/virtual_ko_abs.txt"), sep = "\n")

p <- ggplot(summ, aes(reorder(ko, ratio), ratio, fill = type)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("Writer" = "#7B4EA3", "Eraser" = "#3D9B8C")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(x = NULL, y = "Relative sensitivity (spliceosome / other genes)",
       title = "Virtual KO: spliceosome genes are disproportionately sensitive",
       subtitle = "Absolute predicted shift ratio (GENIE3 unsigned GRN)")
ggsave("figures_final/fig_virtual_ko_abs.pdf", p, width = 7, height = 5)
cat("Done.\n")
