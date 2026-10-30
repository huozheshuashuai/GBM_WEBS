# =============================================================================
# Driver: Pharmacogenomics (DepMap 24Q2 CRISPR + GDSC2)
# =============================================================================
# Reproduces: WEBS of GBM cell lines, dependency correlations (MED1, MLLT1),
#             HDAC inhibitor (vorinostat/SAHA) sensitivity analysis
# Output: results/depmap_*.csv, results/gdsc2_*.csv
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })

WRITERS <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
ERASERS <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")
TARGET_GENES <- c(WRITERS, ERASERS, "MED1", "MLLT1", "SRSF1", "SF3B1", "PRPF8")

dir.create("results", showWarnings = FALSE)

# ---- 1. DepMap: WEBS of GBM cell lines ----
message("Loading DepMap expression (large file, be patient)...")
expr_raw <- read.csv("data/raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                     check.names = FALSE)
# keep default entries only; sample id = ModelID column
meta_cols <- c("", "SequencingID", "ModelConditionID", "ModelID",
               "IsDefaultEntryForMC", "IsDefaultEntryForModel")
if ("IsDefaultEntryForModel" %in% colnames(expr_raw))
  expr_raw <- expr_raw[expr_raw$IsDefaultEntryForModel %in% c("Yes", "TRUE", TRUE), ]
rownames(expr_raw) <- expr_raw$ModelID
expr <- expr_raw[, !colnames(expr_raw) %in% meta_cols]
colnames(expr) <- sub(" .*$", "", colnames(expr))  # "EP300 (8633)" -> "EP300"
rm(expr_raw); gc()
models <- read.csv("data/raw/DepMap/Model.csv", stringsAsFactors = FALSE)

webs_bulk <- function(mat_t) {
  # mat_t: samples x genes
  w <- intersect(WRITERS, colnames(mat_t)); e <- intersect(ERASERS, colnames(mat_t))
  data.frame(sample = rownames(mat_t),
             WEBS = rowMeans(mat_t[, w, drop = FALSE], na.rm = TRUE) -
                    rowMeans(mat_t[, e, drop = FALSE], na.rm = TRUE))
}
webs <- webs_bulk(expr)
write.csv(webs, "results/depmap_cellline_webs.csv", row.names = FALSE)

gbm_lines <- models$ModelID[grepl("glioblastoma|glioma|astrocytoma|gliosarcoma",
                                  models$OncotreePrimaryDisease, ignore.case = TRUE) |
                            grepl("glioblastoma|gliosarcoma",
                                  models$OncotreeSubtype, ignore.case = TRUE)]
gbm_idhwt <- models$ModelID[grepl("glioblastoma|gliosarcoma",
                                  models$OncotreeSubtype, ignore.case = TRUE)]
message("GBM/glioma lines: ", length(gbm_lines), " | IDH-wt GBM: ", length(gbm_idhwt))

# ---- 2. CRISPR dependency correlations (all GBM lines) ----
crispr <- read.csv("data/raw/DepMap/CRISPRGeneEffect.csv", row.names = 1,
                   check.names = FALSE)
colnames(crispr) <- sub(" .*$", "", colnames(crispr))
lines <- intersect(webs$sample, intersect(rownames(crispr), gbm_lines))
message("GBM lines with CRISPR + expression: ", length(lines))

w_sub <- webs$WEBS[match(lines, webs$sample)]
genes_c <- intersect(colnames(crispr), colnames(expr))

cor_all <- function(target_genes, lines_use) {
  wv <- webs$WEBS[match(lines_use, webs$sample)]
  res <- lapply(target_genes, function(g) {
    dep <- crispr[lines_use, g]
    if (all(is.na(dep))) return(NULL)
    ct <- suppressWarnings(cor.test(wv, dep, method = "spearman"))
    data.frame(gene = g, spearman_r = unname(ct$estimate), p = ct$p.value,
               frac_dependent = mean(dep < -1, na.rm = TRUE))  # Chronos < -1 = essential
  })
  do.call(rbind, res)
}

dep_targets <- cor_all(TARGET_GENES, lines)
dep_targets$q <- p.adjust(dep_targets$p, method = "BH")
dep_targets <- dep_targets %>% arrange(p)
write.csv(dep_targets, "results/depmap_webs_dependency.csv", row.names = FALSE)
print(dep_targets)

# ---- 3. GDSC2: HDAC inhibitor sensitivity vs WEBS ----
gdsc <- read.csv("data/raw/GDSC2/GDSC2_fitted_dose_response_24Jul22.csv",
                 stringsAsFactors = FALSE)
# map GDSC cell line -> DepMap ModelID via CellLineName in Model.csv
if ("CellLineName" %in% colnames(models)) {
  m2 <- models[, c("ModelID", "CellLineName")]
  colnames(m2)[2] <- "CELL_LINE_NAME"
} else {
  stop("No CellLineName column in Model.csv")
}

gdsc_w <- gdsc %>%
  inner_join(m2, by = "CELL_LINE_NAME") %>%
  mutate(WEBS = webs$WEBS[match(ModelID, webs$sample)]) %>%
  filter(!is.na(WEBS))

hdac_drugs <- gdsc_w %>%
  filter(grepl("vorinostat|SAHA|romidepsin|panobinostat|belinostat|entinostat|trichostatin",
               DRUG_NAME, ignore.case = TRUE)) %>%
  pull(DRUG_NAME) %>% unique()
message("HDAC inhibitors in GDSC2: ", paste(hdac_drugs, collapse = ", "))

hdac_res <- lapply(hdac_drugs, function(dg) {
  d <- gdsc_w %>% filter(DRUG_NAME == dg)
  if (nrow(d) < 10) return(NULL)
  ct <- suppressWarnings(cor.test(d$WEBS, d$LN_IC50, method = "spearman"))
  data.frame(drug = dg, n = nrow(d), spearman_r = unname(ct$estimate), p = ct$p.value)
})
hdac_res <- do.call(rbind, hdac_res)
if (!is.null(hdac_res)) {
  hdac_res$q <- p.adjust(hdac_res$p, method = "BH")
  write.csv(hdac_res, "results/gdsc2_hdac_sensitivity.csv", row.names = FALSE)
  print(hdac_res)
}

sink("results/pharmacogenomics_summary.txt")
cat("=== Pharmacogenomics summary ===\n\n")
cat("DepMap GBM lines analyzed:", length(lines), "\n\n")
cat("-- WEBS vs CRISPR dependency (Spearman; negative r = high WEBS more dependent) --\n")
print(dep_targets)
cat("\n(Paper: dependencies on MED1 and MLLT1/ENL in writer-dominant state)\n\n")
cat("-- HDAC inhibitors (GDSC2): WEBS vs LN_IC50 --\n")
if (!is.null(hdac_res)) print(hdac_res) else cat("No HDAC drugs with >=10 lines\n")
cat("\n(Paper: high-WEBS lines more sensitive to vorinostat => expect negative r)\n")
sink()

message("Done.")
