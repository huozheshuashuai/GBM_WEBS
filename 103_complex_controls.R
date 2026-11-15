# =============================================================================
# 103. Matched complex controls: is the spliceosome longitudinal decline
#      specific, or shared by other large protein complexes?
#      Complexes from ComplexPortal; match on size/abundance/variance
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(org.Hs.eg.db); library(AnnotationDbi) })
set.seed(42)

b <- readRDS("results/kncc/dual_omics_bundle.rds")
delta <- b$prot_delta
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

# ComplexPortal: extract named complexes
cp <- read.delim("data/raw/ComplexPortal/human_9606.tsv", sep = "\t", quote = "",
                 comment.char = "", check.names = FALSE)
namecol <- "Recommended name"
partcol <- colnames(cp)[grepl("Expanded participant", colnames(cp))]
extract_genes <- function(x) {
  ids <- sub("\\(.*", "", strsplit(x, "\\|")[[1]])
  up <- sub("-.*", "", ids)
  g <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(up),
                    keytype = "UNIPROT", columns = "SYMBOL"))
  unique(na.omit(g$SYMBOL))
}

complex_patterns <- list(
  "Spliceosome"        = NULL,  # from KEGG list
  "Ribosome 60S"       = "^60S|large ribosomal subunit",
  "Ribosome 40S"       = "^40S|small ribosomal subunit",
  "Proteasome 26S/20S" = "proteasome",
  "RNA polymerase II"  = "RNA polymerase II",
  "DNA replication (MCM/RF-C)" = "MCM|replication factor C|RF-C",
  "Mitoribosome"       = "mitochondrial ribosome|mitoribosom",
  "COPI/COPII vesicle" = "COPI|COPII|coatomer"
)

mods <- list(Spliceosome = intersect(splice, rownames(delta)))
for (nm in names(complex_patterns)[-1]) {
  hit <- grepl(complex_patterns[[nm]], cp[[namecol]], perl = TRUE, ignore.case = TRUE)
  if (any(hit)) {
    gs <- unique(unlist(lapply(cp[[partcol]][hit], extract_genes)))
    mods[[nm]] <- intersect(gs, rownames(delta))
  }
}
sizes <- sapply(mods, length)
print(sizes)

# per-complex paired stats
res <- t(sapply(names(mods), function(mn) {
  g <- mods[[mn]]
  if (length(g) < 5) return(c(n = length(g), med = NA, frac_neg = NA, p = NA))
  ds <- colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
  pd <- apply(delta[g, , drop = FALSE], 1, median, na.rm = TRUE)
  wt <- suppressWarnings(wilcox.test(ds, mu = 0))
  c(n = length(g), med = median(pd, na.rm = TRUE),
    frac_neg = mean(pd < 0, na.rm = TRUE), p = wt$p.value)
}))
res <- data.frame(complex = names(mods), res, row.names = NULL)
res$q <- p.adjust(res$p, "BH")
res <- res[order(res$med), ]

sink("results/kncc/complex_controls.txt", split = TRUE)
cat("=== Matched complex controls (KNCC paired T2-T1) ===\n\n")
print(res, row.names = FALSE)
cat("\nInterpretation: if all large complexes decline similarly,\n")
cat("the spliceosome decline is a generic feature, not complex-specific.\n")
sink()

# strict matched permutation for the spliceosome vs other complexes context:
# per-complex z against abundance/variance-matched random sets
expr_mean <- rowMeans(b$prot_samples, na.rm = TRUE)
expr_sd <- apply(b$prot_samples, 1, sd, na.rm = TRUE)
dec <- cut(expr_mean, breaks = quantile(expr_mean, probs = 0:10/10, na.rm = TRUE),
           include.lowest = TRUE, labels = FALSE)
names(dec) <- rownames(b$prot_samples)

matched_null <- function(g, n_iter = 1000) {
  g <- intersect(g, rownames(delta))
  pool <- setdiff(rownames(delta), g)
  out <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    cand <- pool[dec[pool] == dec[g[1]]]  # rough: match first-gene decile; refined below
    gs <- character(0)
    for (gj in g) {
      dj <- dec[gj]
      c2 <- setdiff(pool[dec[pool] == dj], gs)
      if (!length(c2)) c2 <- setdiff(pool, gs)
      pick <- sample(c2, 1)
      gs <- c(gs, pick)
    }
    out[i] <- median(apply(delta[gs, , drop = FALSE], 1, median, na.rm = TRUE), na.rm = TRUE)
  }
  out
}

cat("\nrunning matched nulls (1,000 per complex)...\n")
null_stats <- lapply(names(mods), function(mn) {
  g <- mods[[mn]]
  if (length(g) < 5) return(NULL)
  obs <- res$med[res$complex == mn]
  null <- matched_null(g)
  p <- (sum(abs(null) >= abs(obs)) + 1) / (length(null) + 1)
  data.frame(complex = mn, obs = round(obs, 3),
             null_median = round(median(null), 3), emp_p = signif(p, 3))
})
ns <- do.call(rbind, null_stats)
sink("results/kncc/complex_controls.txt", append = TRUE, split = TRUE)
cat("\n=== matched-null empirical p per complex ===\n")
print(ns, row.names = FALSE)
sink()
write.csv(res, "results/kncc/complex_controls.csv", row.names = FALSE)

res$complex <- factor(res$complex, levels = res$complex)
p <- ggplot(res, aes(complex, med, fill = med)) +
  geom_col() + coord_flip() +
  scale_fill_gradient2(low = "#7B4EA3", mid = "grey90", high = "#3D9B8C") +
  geom_hline(yintercept = 0) + theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "Median per-protein paired delta (T2 - T1, log2)",
       title = "Longitudinal change across protein complexes")
ggsave("figures_kncc/complex_controls.pdf", p, width = 7, height = 4.5)
cat("Done.\n")
