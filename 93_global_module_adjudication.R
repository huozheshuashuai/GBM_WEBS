# =============================================================================
# 93. KNCC module-level competing-hypothesis screen (whole proteome)
#     For each Hallmark module: delta_M ~ delta_prolif* + delta_neuron + delta_immune
#     (predictor scores exclude the module's own genes to avoid circularity)
#     Classes: A = fully explained by cell state | B = independent residual | C = mixed
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta   # proteins x 105 patients

# Hallmark modules
gmt <- readLines("data/raw/msigdb_hallmark.gmt")
hm <- lapply(gmt, function(l) { f <- strsplit(l, "\t")[[1]]; f[-(1:2)] })
names(hm) <- sub("\t.*", "", gmt)
hm <- lapply(hm, function(g) intersect(g, rownames(delta)))
hm <- hm[sapply(hm, length) >= 10]
cat("Hallmark modules with >=10 quantified proteins:", length(hm), "\n")

prolif_core <- intersect(c("MKI67","PCNA","TOP2A","MCM2","MCM3","MCM4","MCM5","MCM6",
                           "MCM7","AURKB","CCNB1","CDK1","BUB1","TTK"), rownames(delta))
neuron <- intersect(c("SYP","SNAP25","SYT1","NEFL","NEFM","RBFOX3","STMN2",
                      "GAP43","TUBB3","MAP2"), rownames(delta))
immune <- intersect(c("CD14","AIF1","FCGR3A","CD68","CD163","ITGAM","PTPRC",
                      "CSF1R","TMEM119","CX3CR1"), rownames(delta))

out <- lapply(names(hm), function(mn) {
  g <- hm[[mn]]
  d_m <- colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
  # predictors excluding this module's genes
  d_p <- colMeans(delta[setdiff(prolif_core, g), , drop = FALSE], na.rm = TRUE)
  d_n <- colMeans(delta[setdiff(neuron, g), , drop = FALSE], na.rm = TRUE)
  d_i <- colMeans(delta[setdiff(immune, g), , drop = FALSE], na.rm = TRUE)
  ok <- !is.na(d_m) & !is.na(d_p) & !is.na(d_n) & !is.na(d_i)
  if (sum(ok) < 50) return(NULL)
  fit <- lm(d_m[ok] ~ d_p[ok] + d_n[ok] + d_i[ok])
  cf <- summary(fit)$coefficients
  ci <- confint(fit)
  wt <- suppressWarnings(wilcox.test(d_m[ok], mu = 0))
  # patient direction consistency
  frac_neg <- mean(d_m[ok] < 0)
  data.frame(module = mn, n_genes = length(g),
             raw_median = median(d_m[ok]), raw_p = wt$p.value,
             adj_intercept = cf[1, 1], adj_lo = ci[1, 1], adj_hi = ci[1, 2], adj_p = cf[1, 4],
             beta_prolif = cf[2, 1], beta_prolif_p = cf[2, 4],
             beta_neuron = cf[3, 1], beta_neuron_p = cf[3, 4],
             beta_immune = cf[4, 1], beta_immune_p = cf[4, 4],
             r2 = summary(fit)$r.squared, frac_neg = frac_neg)
})
res <- do.call(rbind, out)
res$adj_q <- p.adjust(res$adj_p, "BH")
res$raw_q <- p.adjust(res$raw_p, "BH")
# classification: B = raw sig AND adjusted sig; A = raw sig, adjusted not; C = raw not sig
res$class <- ifelse(res$raw_q < 0.05 & res$adj_q < 0.05, "B (independent residual)",
             ifelse(res$raw_q < 0.05, "A (explained by cell state)", "C (n.s.)"))
res <- res[order(res$adj_intercept), ]
write.csv(res, "results/kncc/module_adjudication.csv", row.names = FALSE)

sink("results/kncc/module_adjudication.txt", split = TRUE)
cat("=== Module-level competing-hypothesis screen (n = 105 paired patients) ===\n\n")
cat("class A (cell-state explained):", sum(res$class == "A (explained by cell state)"), "\n")
cat("class B (independent residual):", sum(res$class == "B (independent residual)"), "\n")
cat("class C (not significant):", sum(res$class == "C (n.s.)"), "\n\n")
cat("--- class B modules (independent residual after cell-state adjustment) ---\n")
print(res[res$class == "B (independent residual)",
          c("module", "n_genes", "raw_median", "adj_intercept", "adj_q", "r2", "frac_neg")],
      row.names = FALSE)
cat("\n--- full table ---\n")
print(res[, c("module", "n_genes", "raw_median", "raw_q", "adj_intercept", "adj_q", "class")],
      row.names = FALSE)
sink()

# figure: raw vs adjusted per module
res$lab <- ifelse(res$class == "B (independent residual)", sub("HALLMARK_", "", res$module), NA)
p <- ggplot(res, aes(raw_median, adj_intercept, color = class)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.8, max.overlaps = 30) +
  scale_color_manual(values = c("A (explained by cell state)" = "grey60",
                                "B (independent residual)" = "#7B4EA3",
                                "C (n.s.)" = "grey85")) +
  theme_bw(base_size = 11) +
  labs(x = "Raw module delta (T2 - T1, median)",
       y = "Cell-state-adjusted effect (intercept)",
       title = "Which proteome modules change independently of cell state?",
       color = NULL)
ggsave("figures_kncc/global_module_adjudication.pdf", p, width = 8, height = 6)
cat("Done.\n")
