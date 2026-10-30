# =============================================================================
# 24. GDSC2 combination (Bliss) prediction + CMap-style drug repurposing
# =============================================================================
# (a) Bliss independence prediction from single-agent GDSC2 responses:
#     pairs of drugs where WEBS correlates with sensitivity (r < -0.2)
# (b) CMap-style reversal query: genes up in LHS (CPTAC limma) vs
#     drug-induced signatures approximated by GDSC2 drug-target annotation
# Output: results/gdsc2_combination.csv, results/cmap_candidates.csv
# =============================================================================

setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

# ---- (a) Bliss combination prediction ----
gdsc <- read.csv("data/raw/GDSC2/GDSC2_fitted_dose_response_24Jul22.csv",
                 stringsAsFactors = FALSE)
webs <- NULL
if (file.exists("results/depmap_cellline_webs.csv")) {
  webs <- read.csv("results/depmap_cellline_webs.csv")
} else {
  message("Recomputing cell-line WEBS (expression subset)...")
  expr_raw <- read.csv("data/raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                       check.names = FALSE, nrows = 5)
  meta_cols <- c("SequencingID", "ModelConditionID", "ModelID",
                 "IsDefaultEntryForMC", "IsDefaultEntryForModel")
  WRITERS <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
  ERASERS <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7")
  gene_cols <- colnames(expr_raw)[sub(" .*$", "", colnames(expr_raw)) %in% c(WRITERS, ERASERS)]
  expr <- read.csv("data/raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                   check.names = FALSE,
                   colClasses = c(rep("NULL", 1), rep("NULL", 2), "character",
                                  rep("NULL", 1), "NULL",  # keep ModelID col
                                  rep(NA, ncol(expr_raw) - 6)))
  stop("Run run_10_pharmaco.R first to generate results/depmap_cellline_webs.csv")
}

models <- read.csv("data/raw/DepMap/Model.csv", stringsAsFactors = FALSE)
m2 <- models[, c("ModelID", "CellLineName")]
colnames(m2)[2] <- "CELL_LINE_NAME"

d <- gdsc %>%
  inner_join(m2, by = "CELL_LINE_NAME") %>%
  mutate(WEBS = webs$WEBS[match(ModelID, webs$sample)]) %>%
  filter(!is.na(WEBS))

# per-drug WEBS correlation
drug_cor <- d %>%
  group_by(DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME) %>%
  filter(n() >= 30) %>%
  summarise(
    n = n(),
    spearman_r = suppressWarnings(cor(WEBS, LN_IC50, method = "spearman",
                                      use = "pairwise.complete.obs")),
    p = suppressWarnings(cor.test(WEBS, LN_IC50, method = "spearman")$p.value),
    mean_auc = mean(AUC, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(q = p.adjust(p, method = "BH"))
write.csv(drug_cor, "results/gdsc2_drug_webs_corr.csv", row.names = FALSE)

# candidate sensitizers in high-WEBS lines (r < -0.2)
sens <- drug_cor %>% filter(spearman_r < -0.2, q < 0.1) %>% arrange(spearman_r)
cat("Sensitizer candidates (r < -0.2, FDR<0.1):", nrow(sens), "\n")
print(head(sens, 15))

# Bliss prediction for pairs among top sensitizers
top_drugs <- head(sens$DRUG_NAME, 8)
if (length(top_drugs) >= 2) {
  wide <- d %>% filter(DRUG_NAME %in% top_drugs) %>%
    select(CELL_LINE_NAME, DRUG_NAME, AUC, WEBS) %>%
    tidyr::pivot_wider(names_from = DRUG_NAME, values_from = AUC,
                       values_fn = mean)
  pairs <- combn(top_drugs, 2, simplify = FALSE)
  bliss <- lapply(pairs, function(pr) {
    a <- wide[[pr[1]]]; b <- wide[[pr[2]]]
    ok <- !is.na(a) & !is.na(b)
    if (sum(ok) < 20) return(NULL)
    # viability proxy: 1 - AUC (higher = more killing); Bliss: Eab = Ea + Eb - Ea*Eb
    ea <- 1 - a[ok]; eb <- 1 - b[ok]
    bliss_expected <- ea + eb - ea * eb
    # synergy signal: correlation of expected combo efficacy with WEBS
    ct <- suppressWarnings(cor.test(bliss_expected, wide$WEBS[ok], method = "spearman"))
    data.frame(drug1 = pr[1], drug2 = pr[2], n = sum(ok),
               bliss_webs_rho = unname(ct$estimate), p = ct$p.value)
  })
  bliss <- do.call(rbind, bliss)
  bliss$q <- p.adjust(bliss$p, method = "BH")
  bliss <- bliss %>% arrange(desc(bliss_webs_rho))
  write.csv(bliss, "results/gdsc2_combination.csv", row.names = FALSE)
  cat("\nBliss combination predictions (high rho = combo works better in high-WEBS):\n")
  print(bliss)
}

# ---- (b) CMap-style: drug-target reversal ----
# LHS-up genes (CPTAC) that are targets of drugs with WEBS-sensitivity
de <- read.csv("results/cptac_limma_results.csv")
lhs_up <- de$protein[!is.na(de$adj.P.Val) & de$adj.P.Val < 0.05 & de$logFC > 0]

drug_targets <- gdsc %>% select(DRUG_NAME, PUTATIVE_TARGET) %>% distinct()
cmap <- drug_targets %>%
  filter(PUTATIVE_TARGET != "") %>%
  tidyr::separate_rows(PUTATIVE_TARGET, sep = ",\\s*") %>%
  filter(PUTATIVE_TARGET %in% lhs_up) %>%
  inner_join(drug_cor, by = "DRUG_NAME", suffix = c("", ".cor")) %>%
  filter(spearman_r < 0) %>%
  arrange(spearman_r)
write.csv(cmap, "results/cmap_candidates.csv", row.names = FALSE)

sink("results/translation_summary.txt")
cat("=== Translational analysis summary ===\n\n")
cat("Drugs tested:", nrow(drug_cor), "\n")
cat("WEBS-sensitizing drugs (r < -0.2, FDR < 0.1):", nrow(sens), "\n")
print(head(sens, 10))
cat("\nDrugs targeting LHS-up proteins and WEBS-sensitizing:\n")
print(head(cmap[, c("DRUG_NAME", "PUTATIVE_TARGET", "spearman_r", "q")], 15))
sink()
message("Done.")
