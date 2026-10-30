# RPPA500-based LMS survival in TCGA-GBM (TCPA platform, samples x antibodies)
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival); library(survminer); library(dplyr) })

rppa <- read.delim("data/raw/LinkedOmics_TCGA/TCPA_TCGA_RPPA500.tsv",
                   check.names = FALSE, row.names = 1)
rppa <- as.matrix(rppa)  # samples x antibodies
lms_ab <- c("LDHA", "LDHB", "PKM2", "MCT4", "HIF1ALPHA", "ENOLASE1", "ENOLASE2", "GAPDH")
lms_ab <- intersect(lms_ab, colnames(rppa))
cat("LMS antibodies found:", length(lms_ab), "/", 8, "\n")

lms <- rowMeans(rppa[, lms_ab, drop = FALSE], na.rm = TRUE)
names(lms) <- rownames(rppa)

# match to clinical (sample IDs: TCGA-02-0003-01A vs TCGA.02.0003)
clin <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
short <- function(x) substr(gsub("[\\.]", "-", x), 1, 12)
clin$short <- short(clin$sample)
lms_df <- data.frame(short = short(names(lms)), LMS = as.numeric(lms))
clin <- clin %>% inner_join(lms_df, by = "short") %>% filter(!is.na(LMS))

cat("Matched TCGA-GBM samples with RPPA500 + survival:", nrow(clin), "\n")

q <- quantile(clin$LMS, c(0.25, 0.75))
clin$LMS_grp <- factor(ifelse(clin$LMS <= q[1], "Low", ifelse(clin$LMS >= q[2], "High", "Mid")),
                       levels = c("Low", "Mid", "High"))

sink("results/rppa_lms_summary.txt", split = TRUE)
cat("=== RPPA500 LMS survival (TCGA-GBM) ===\n")
cat("n =", nrow(clin), "\n\n")
cat("KM (Low/Mid/High):\n")
print(surv_pvalue(survfit(Surv(OS.time, OS) ~ LMS_grp, data = clin), data = clin))
cat("\nCox continuous LMS:\n")
print(summary(coxph(Surv(OS.time, OS) ~ LMS, data = clin))$coefficients)
cat("\nCox adjusted for age:\n")
print(summary(coxph(Surv(OS.time, OS) ~ LMS + age, data = clin))$coefficients)
sink()
