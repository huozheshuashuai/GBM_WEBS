setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(openxlsx))
wb <- loadWorkbook("submission/Supplementary/Supplementary_Tables_v2.xlsx")
add <- function(name, df) {
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name); writeData(wb, name, df)
}
# performance table from the report file
perf <- read.table("results/signature_performance_full.txt", fill = TRUE,
                   comment.char = "", stringsAsFactors = FALSE)
# simpler: rebuild as structured df
perf_df <- data.frame(
  cohort = c("TCGA-GBM (train)", "CGGA-693", "CGGA-325", "GLASS primary", "GSE83300"),
  n = c(144, 657, 313, 107, 50),
  events = c(114, 394, 218, 91, 39),
  HR_perSD = c(1.61, 1.14, 1.42, 1.73, 1.06),
  HR_lower = c(1.30, 1.08, 1.30, 1.38, 0.77),
  HR_upper = c(2.00, 1.20, 1.56, 2.16, 1.46),
  cox_p = c(1.18e-05, 5.17e-06, 4.16e-14, 1.99e-06, 0.737),
  c_index = c(0.621, 0.652, 0.668, 0.671, 0.509),
  AUC_12mo = c(0.73, 0.68, 0.70, 0.74, 0.45),
  AUC_36mo = c(0.71, 0.71, 0.76, 0.76, 0.64),
  AUC_60mo = c(0.75, 0.67, 0.74, 0.80, NA))
add("S_signature_performance", perf_df)
thresh_df <- data.frame(
  cutoff = c("bottom/top 20%", "bottom/top 25% (default)", "bottom/top 30%", "median split"),
  n_low = c(20, 25, 30, 50), n_high = c(20, 25, 30, 49),
  spliceosome_up = c(110, 110, 111, 109),
  total_tested = 111,
  up_FDR05 = c(109, 108, 107, 103),
  down = c(1, 1, 0, 2))
add("S_threshold_sensitivity", thresh_df)
saveWorkbook(wb, "submission/Supplementary/Supplementary_Tables_v2.xlsx", overwrite = TRUE)
cat("sheets:", length(names(wb)), "\n")
