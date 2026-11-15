# =============================================================================
# 88. KNCC module 3: subcomplex-resolved spliceosome remodeling (T1->T2)
#     Module annotation: ComplexPortal 9606 (traceable), overlapping allowed
#     Key question: is the decline uniform or module-selective?
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2); library(org.Hs.eg.db);
                   library(AnnotationDbi); library(patchwork) })
set.seed(42)

obj <- readRDS("results/kncc/proteome_analysis_ready.rds")
pm <- obj$pm; cli <- obj$annot; delta <- obj$delta; paired_pts <- obj$paired
res <- read.csv("results/kncc/paired_differential_proteins.csv")

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
sp <- intersect(splice, rownames(pm))

# ---- parse ComplexPortal human complexes ----
cp <- read.delim("data/raw/ComplexPortal/human_9606.tsv", sep = "\t",
                 quote = "", comment.char = "", check.names = FALSE)
namecol <- "Recommended name"
partcol <- colnames(cp)[grepl("Expanded participant", colnames(cp))]

extract_genes <- function(x) {
  ids <- sub("\\(.*", "", strsplit(x, "\\|")[[1]])
  up <- sub("-.*", "", ids)
  g <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(up),
                    keytype = "UNIPROT", columns = "SYMBOL"))
  unique(na.omit(g$SYMBOL))
}

# spliceosome-related complex patterns (module assignment by curated name rules)
patterns <- list(
  "U1 snRNP"            = "^U1 (small nuclear|snRNP)|U1 snRNP",
  "U2 snRNP"            = "^U2 (small nuclear|snRNP)|U2 snRNP|SF3[ab] complex",
  "U4/U6 snRNP"         = "U4/U6(?!\\.U5)|U4/U6 di",
  "U5 snRNP"            = "^(U5 snRNP|U5 small nuclear)",
  "U4/U6.U5 tri-snRNP"  = "tri-snRNP|U4/U6\\.U5|U4/U6-U5",
  "NineTeen/Prp19"      = "NineTeen|NTC|Prp19|PRP19|CDC5L complex",
  "Catalytic spliceosome" = "[Cc]atalytic (step|spliceosome)|spliceosomal [BCP] complex|spliceosome [BC] "
)

modules <- list()
for (m in names(patterns)) {
  hit <- grepl(patterns[[m]], cp[[namecol]], perl = TRUE)
  if (any(hit)) {
    gs <- unique(unlist(lapply(cp[[partcol]][hit], extract_genes)))
    modules[[m]] <- intersect(gs, rownames(pm))
  }
}
# family modules (HUGO nomenclature-based, documented)
modules[["SR proteins (SRSF)"]]   <- intersect(grep("^SRSF", rownames(pm), value = TRUE), rownames(pm))
modules[["hnRNP proteins"]]       <- grep("^HNRNP", rownames(pm), value = TRUE)
modules[["DDX/DHX helicases"]]    <- intersect(grep("^(DDX|DHX)", rownames(pm), value = TRUE), sp)

# keep modules with >=3 quantified members, tag overlap with curated spliceosome set
mod_df <- data.frame(module = names(modules),
                     n = sapply(modules, length),
                     n_splice = sapply(modules, function(g) sum(g %in% sp)))
print(mod_df)

# ---- per-module stats ----
res_idx <- res[match(rownames(delta), res$gene), ]
prot_delta <- res$med[match(rownames(delta), res$gene)]

mod_stats <- lapply(names(modules), function(m) {
  g <- intersect(modules[[m]], rownames(delta))
  d <- prot_delta[rownames(delta) %in% g]
  # module score per patient = mean delta of members; paired Wilcoxon vs 0
  ds <- colMeans(delta[g, , drop = FALSE], na.rm = TRUE)
  wt <- suppressWarnings(wilcox.test(ds, mu = 0))
  data.frame(module = m, n = length(g),
             median_delta = median(d, na.rm = TRUE),
             mean_delta = mean(d, na.rm = TRUE),
             frac_neg = mean(d < 0, na.rm = TRUE),
             frac_fdr = mean(res$q[match(g, res$gene)] < 0.05, na.rm = TRUE),
             score_p = wt$p.value)
})
ms <- do.call(rbind, mod_stats)
ms$score_q <- p.adjust(ms$score_p, "BH")
ms <- ms[order(ms$median_delta), ]
print(ms)

# ---- KEY TEST: is remodeling module-selective? ----
# Kruskal-Wallis on per-protein delta across modules (proteins in >1 module counted once, first module)
prot2mod <- list()
for (m in names(modules)) for (g in modules[[m]]) if (is.null(prot2mod[[g]])) prot2mod[[g]] <- m
p2m <- unlist(prot2mod)
kd <- data.frame(gene = names(p2m), module = p2m,
                 delta = prot_delta[match(names(p2m), rownames(delta))])
kd <- kd[!is.na(kd$delta), ]
kw <- kruskal.test(delta ~ module, data = kd)
cat("\nKruskal-Wallis across modules (selectivity test): p =", signif(kw$p.value, 3), "\n")
# compare vs genome background
wt_bg <- wilcox.test(kd$delta, prot_delta[!rownames(delta) %in% names(p2m)])
cat("all-module proteins vs genome background: p =", signif(wt_bg$p.value, 3), "\n")

write.csv(ms, "results/kncc/subcomplex_module_stats.csv", row.names = FALSE)
write.csv(kd, "results/kncc/subcomplex_per_protein.csv", row.names = FALSE)

# ---- figures ----
ms$module <- factor(ms$module, levels = ms$module)
p_bar <- ggplot(ms, aes(module, median_delta, fill = median_delta)) +
  geom_col() + coord_flip() +
  scale_fill_gradient2(low = "#7B4EA3", mid = "grey90", high = "#3D9B8C") +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.1, size = 3) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = "Median paired delta (T2 - T1, log2)",
       title = "Subcomplex-resolved spliceosome remodeling",
       subtitle = paste0("module heterogeneity KW p = ", signif(kw$p.value, 2)))
ggsave("figures_kncc/m3_subcomplex_deltas.pdf", p_bar, width = 7, height = 5)

p_box <- ggplot(kd, aes(reorder(module, delta, median), delta)) +
  geom_boxplot(fill = "grey85", outlier.size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + theme_bw(base_size = 11) +
  labs(x = NULL, y = "Per-protein paired delta (T2 - T1, log2)",
       title = "Per-protein deltas by spliceosome module")
ggsave("figures_kncc/m3_subcomplex_boxplot.pdf", p_box, width = 7, height = 5)
cat("Done.\n")
