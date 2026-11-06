setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(openxlsx))
wb <- loadWorkbook("submission/Supplementary/Supplementary_Tables_v2.xlsx")
if ("S_protein_WEBS" %in% names(wb)) removeWorksheet(wb, "S_protein_WEBS")
addWorksheet(wb, "S_protein_WEBS")
df <- data.frame(
  protein = c("HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT7","KAT8","CREBBP","EP300","KAT14",
              "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10",
              "HDAC11","SIRT1","SIRT2","SIRT3","SIRT5","SIRT6"),
  class = c(rep("Writer", 10), rep("Eraser", 16)),
  note = "All writer/eraser proteins quantified in the CPTAC TMT panel; no imputation applied")
writeData(wb, "S_protein_WEBS", df)
saveWorkbook(wb, "submission/Supplementary/Supplementary_Tables_v2.xlsx", overwrite = TRUE)
cat("sheets:", length(names(wb)), "\n")
