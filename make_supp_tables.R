# Package supplementary tables into a single xlsx
setwd("E:/KIMI/GBM_WEBS_Project")
if (!require("writexl", quietly = TRUE)) install.packages("writexl", repos = "https://cloud.r-project.org/")
library(writexl)

read_safe <- function(f) if (file.exists(f)) read.csv(f) else NULL

tables <- list(
  S_CPTAC_WEBS = read_safe("results/cptac_webs.csv"),
  S_CPTAC_limma = read_safe("results/cptac_limma_results.csv"),
  S_spliceosome = read_safe("results/cptac_spliceosome_proteins.csv"),
  S_acetyl_sites = read_safe("results/acetyl_sites_limma.csv"),
  S_acetyl_corr = read_safe("results/acetyl_site_webs_corr.csv"),
  S_bulk_survival = read_safe("results/bulk_survival_summary.csv"),
  S_signature_valid = read_safe("results/sixgene_validation_summary.csv"),
  S_bootstrap = read_safe("results/bootstrap_selection_frequency.csv"),
  S_SHAP = read_safe("results/shap_summary.csv"),
  S_pancancer = read_safe("results/pancancer_webs_summary.csv"),
  S_benchmark = read_safe("results/benchmarking_correlation.csv"),
  S_performance = read_safe("results/computational_performance.csv"),
  S_sensitivity_maxrank = read_safe("results/sensitivity_maxrank.csv"),
  S_sensitivity_loo = read_safe("results/sensitivity_leaveout.csv"),
  S_specificity = read_safe("results/specificity_control.csv"),
  S_spatial = read_safe("results/spatial_autocorrelation.csv"),
  S_MOFA2 = read_safe("results/mofa2_factor_webs_corr.csv"),
  S_DepMap = read_safe("results/depmap_webs_dependency.csv"),
  S_GDSC2_drugs = read_safe("results/gdsc2_drug_webs_corr.csv"),
  S_GDSC2_HDAC = read_safe("results/gdsc2_hdac_sensitivity.csv"),
  S_combinations = read_safe("results/gdsc2_combination.csv"),
  S_CMap = read_safe("results/cmap_candidates.csv"),
  S_CellChat = read_safe("results/cellchat_pathway_comparison.csv"),
  S_MIF_perm = read_safe("results/mif_permutation.csv"),
  S_docking = read_safe("docking/results/binding_energies.csv"),
  S_ablation = read_safe("results/ablation_sixgene.csv")
)
tables <- tables[!sapply(tables, is.null)]

dir.create("submission", showWarnings = FALSE)
write_xlsx(tables, "submission/Supplementary_Tables.xlsx")
cat("Supplementary_Tables.xlsx with", length(tables), "sheets\n")
