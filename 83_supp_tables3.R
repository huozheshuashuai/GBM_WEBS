setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(openxlsx))
wb <- loadWorkbook("submission_GND/04_Supplementary/Supplementary_Tables_GND.xlsx")
add <- function(name, df) {
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name); writeData(wb, name, df)
}
add("S_rna_protein_residuals", read.csv("results/rna_protein_residual_perprotein.csv"))
add("S_splicing_confounders", read.csv("results/splicing_confounders.csv"))
saveWorkbook(wb, "submission_GND/04_Supplementary/Supplementary_Tables_GND.xlsx", overwrite = TRUE)
cat("sheets:", length(names(wb)), "\n")
