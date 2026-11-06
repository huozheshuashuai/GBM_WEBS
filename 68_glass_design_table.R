setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(openxlsx))
wb <- loadWorkbook("submission/Supplementary/Supplementary_Tables_v2.xlsx")
if ("S_glass_design" %in% names(wb)) removeWorksheet(wb, "S_glass_design")
addWorksheet(wb, "S_glass_design")
df <- data.frame(
  analysis = c("Primary-tumor survival (signature/WEBS)",
               "Primary-recurrence WEBS change (paired)",
               "Primary WEBS vs recurrence interval",
               "First-pair-per-patient sensitivity"),
  n = c("107 primary tumors", "179 pairs", "213 pairs", "173 patients"),
  unit = c("patient", "pair", "pair (some patients contribute >1)",
           "independent patient"),
  inclusion = c("primary (TP) aliquots with survival data",
                "longitudinal TP-R1 RNA pairs with both aliquots",
                "longitudinal pairs with recorded surgical interval > 0",
                "first pair per patient only"),
  stringsAsFactors = FALSE)
writeData(wb, "S_glass_design", df)
saveWorkbook(wb, "submission/Supplementary/Supplementary_Tables_v2.xlsx", overwrite = TRUE)
cat("sheets:", length(names(wb)), "\n")
