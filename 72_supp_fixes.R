# 72. Supplementary table fixes: add S_ablation, full splicing events
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(openxlsx))

for (path in c("submission/Supplementary/Supplementary_Tables_v2.xlsx",
               "submission_GND/04_Supplementary/Supplementary_Tables_GND.xlsx")) {
  wb <- loadWorkbook(path)
  add <- function(name, df) {
    if (name %in% names(wb)) removeWorksheet(wb, name)
    addWorksheet(wb, name); writeData(wb, name, df)
  }
  # S_ablation (was missing)
  add("S_ablation", read.csv("results/ablation_sixgene.csv"))
  # full significant splicing events (replace top-500)
  sp <- read.csv("results/splicing_webs_corr.csv")
  sig <- sp[sp$q < 0.05 & abs(sp$rho) > 0.2, ]
  sig <- sig[order(-abs(sig$rho)), ]
  add("S_splicing_sig4387", sig)
  saveWorkbook(wb, path, overwrite = TRUE)
  cat(path, "->", length(names(wb)), "sheets\n")
}
