# =============================================================================
# 100. WEBS-free therapeutic analysis: spliceosome score vs DepMap dependency
#      and GDSC2 drug sensitivity (brain/GBM cell lines)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(42)

source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()

# ---------- DepMap ----------
crispr <- read.csv("data/raw/DepMap/CRISPRGeneEffect.csv", row.names = 1, check.names = FALSE)
colnames(crispr) <- sub(" \\(.*", "", colnames(crispr))
models <- read.csv("data/raw/DepMap/Model.csv")
expr_raw <- read.csv("data/raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                     row.names = 1, check.names = FALSE)
expr_raw <- expr_raw[expr_raw$IsDefaultEntryForModel == "Yes", ]
expr <- expr_raw[, !colnames(expr_raw) %in% c("SequencingID","ModelConditionID","ModelID",
                                              "IsDefaultEntryForMC","IsDefaultEntryForModel")]
rownames(expr) <- expr_raw$ModelID
expr <- expr[!is.na(rownames(expr)), ]
colnames(expr) <- sub(" \\(.*", "", colnames(expr))

gbm_lines <- models$ModelID[grepl("Brain", models$OncotreeLineage, ignore.case = TRUE)]
cat("brain lines in Model.csv:", length(gbm_lines), "\n")

common <- intersect(intersect(rownames(crispr), rownames(expr)), gbm_lines)
cat("brain lines with CRISPR + expression:", length(common), "\n")

sp_g <- intersect(splice, colnames(expr))
sps <- rowMeans(expr[common, sp_g, drop = FALSE])
names(sps) <- common

# correlate score with dependency of every gene (negative = more essential)
cor_dep <- function(g) {
  v <- crispr[common, g]
  ok <- !is.na(v)
  if (sum(ok) < 15 || sd(v[ok]) == 0) return(c(NA, NA))
  ct <- suppressWarnings(cor.test(sps[ok], v[ok], method = "spearman"))
  c(unname(ct$estimate), ct$p.value)
}
sp_dep <- intersect(splice, colnames(crispr))
rr <- t(sapply(sp_dep, cor_dep))
dep_df <- data.frame(gene = sp_dep, rho = rr[, 1], p = rr[, 2])
dep_df$q <- p.adjust(dep_df$p, "BH")
dep_df <- dep_df[order(dep_df$rho), ]
write.csv(dep_df, "results/reframe/depmap_spliceosome_dependency.csv", row.names = FALSE)

sink("results/reframe/depmap_gDSC2.txt", split = TRUE)
cat("=== DepMap: spliceosome score vs spliceosome-gene dependency (brain lines) ===\n")
cat("lines:", length(common), "| spliceosome genes tested:", nrow(dep_df), "\n")
cat("FDR<0.05:", sum(dep_df$q < 0.05, na.rm = TRUE),
    " | negative rho (more dependent in high-score):", sum(dep_df$rho < 0 & dep_df$q < 0.05, na.rm = TRUE), "\n")
cat("median rho across spliceosome genes:", round(median(dep_df$rho, na.rm = TRUE), 3), "\n")
cat("\ntop 15 (most negative rho):\n")
print(head(dep_df[, c("gene", "rho", "q")], 15), row.names = FALSE)

# one-sample test: are spliceosome-dependency correlations shifted negative?
wt <- wilcox.test(dep_df$rho, mu = 0)
cat("\none-sample Wilcoxon on rho (spliceosome deps): p =", signif(wt$p.value, 3), "\n")
sink()

# ---------- GDSC2 ----------
gdsc <- read.csv("data/raw/GDSC2/GDSC2_fitted_dose_response_24Jul22.csv")
# GDSC2 uses SANGER_MODEL_ID / COSMIC_ID; map via Model.csv SangerModelID
gdsc_brain <- gdsc[gdsc$SANGER_MODEL_ID %in% models$SangerModelID[models$ModelID %in% common], ]
map_sid <- setNames(models$ModelID, models$SangerModelID)
gdsc_brain$ModelID <- map_sid[gdsc_brain$SANGER_MODEL_ID]
gdsc_brain <- gdsc_brain[!is.na(gdsc_brain$ModelID), ]
gdsc_brain$sps <- sps[gdsc_brain$ModelID]

drug_res <- gdsc_brain %>%
  group_by(DRUG_NAME) %>%
  filter(n() >= 15) %>%
  summarise(n = n(),
            rho = suppressWarnings(cor(LN_IC50, sps, method = "spearman", use = "complete.obs")),
            p = tryCatch(suppressWarnings(cor.test(LN_IC50, sps, method = "spearman"))$p.value,
                         error = function(e) NA), .groups = "drop")
drug_res$q <- p.adjust(drug_res$p, "BH")
drug_res <- drug_res[order(drug_res$p), ]
write.csv(drug_res, "results/reframe/gdsc2_spliceosome_drugs.csv", row.names = FALSE)

sink("results/reframe/depmap_gDSC2.txt", append = TRUE, split = TRUE)
cat("\n=== GDSC2: spliceosome score vs drug sensitivity (brain lines) ===\n")
cat("drugs tested:", nrow(drug_res), "\n")
cat("FDR<0.1:", sum(drug_res$q < 0.1, na.rm = TRUE), "\n")
cat("negative rho = sensitivity in high-score lines\n")
cat("\ntop 15 drugs:\n")
print(head(drug_res[, c("DRUG_NAME", "n", "rho", "p", "q")], 15), row.names = FALSE)
sink()
cat("Done.\n")
