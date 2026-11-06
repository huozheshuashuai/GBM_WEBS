# 71. Finalize GND package: add missing supp tables + 5-folder structure
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(openxlsx))

for (path in c("submission/Supplementary/Supplementary_Tables_v2.xlsx",
               "submission_GND/Supplementary/Supplementary_Tables_GND.xlsx")) {
  wb <- loadWorkbook(path)
  add <- function(name, df) {
    if (name %in% names(wb)) removeWorksheet(wb, name)
    addWorksheet(wb, name); writeData(wb, name, df)
  }
  # Table 2 as supplementary sheet
  add("S_Table2_targets", data.frame(
    Tier = c(1, 2, 2, 2, 3, 3),
    Candidate = c("HDAC2", "SF3B1", "MED1", "SIRT2", "CLK3/CDK9/12/13", "KAT2A/EP300"),
    Key_evidence = c(
      "GDSC2 vorinostat q = 1.1e-4; docking -7.9 kcal/mol; strong spliceosome connectivity shift",
      "Spliceosome-complex up-regulation (107/111); DepMap nominal p = 0.038; docking -6.7 kcal/mol",
      "DepMap CRISPR r = -0.28, q = 0.049",
      "Signature member; docking AGK2 -8.5 kcal/mol",
      "Splicing-kinase abundance tracks WEBS (rho 0.49-0.56)",
      "Writer core; perturbation sensitivity"),
    Tool_compounds = c("Vorinostat, entinostat", "H3B-8800, E7107, pladienolide B",
                       "(preclinical)", "AGK2", "CLK inhibitors (preclinical)", "C646, MG149"),
    note = "Tier 1 = highest priority for experimental validation, not validated efficacy"))
  # GLASS design table
  add("S_glass_design", data.frame(
    analysis = c("Primary-tumor survival (signature/WEBS)",
                 "Primary-recurrence WEBS change (paired)",
                 "Primary WEBS vs recurrence interval",
                 "First-pair-per-patient sensitivity"),
    n = c("107 primary tumors", "179 pairs (148 patients)", "213 pairs (173 patients)", "173 patients"),
    unit = c("patient", "pair", "pair (some patients contribute >1)", "independent patient"),
    inclusion = c("primary (TP) aliquots with survival data",
                  "longitudinal TP-R1 RNA pairs with both aliquots quantified",
                  "longitudinal pairs with recorded surgical interval > 0",
                  "first pair per patient only")))
  saveWorkbook(wb, path, overwrite = TRUE)
  cat(path, "->", length(names(wb)), "sheets\n")
}
