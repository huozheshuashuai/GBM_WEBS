# one-off: dump HHS matrix to plain CSV for format-independent verification
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(readxl))
raw <- read_excel("E:/KIMI/data/HHS_Tatari2022/NIHMS1894660-supplement-Supplementary_Table_2.xlsx",
                  sheet = "Imputed values", col_names = FALSE, skip = 3)
colnames(raw) <- as.character(raw[1, ])
raw <- raw[-1, ]
write.csv(raw, "results/hhs/hhs_matrix_plain.csv", row.names = FALSE)
cat("written:", nrow(raw), "x", ncol(raw), "\n")
