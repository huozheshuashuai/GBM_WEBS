setwd("E:/KIMI/GBM_WEBS_Project")
t1 <- data.frame(
  Cohort = c("KNCC (discovery)","HHS (validation)","CPTAC-GBM","GeoMx GSE254875","GLASS","GSE155434"),
  Type = c("Longitudinal proteome+phosphoproteome","Longitudinal proteome (FFPE, LFQ)","Cross-sectional proteome/PTM/RNA","Spatial WTA (79 ROIs)","Longitudinal RNA-seq","Longitudinal RNA-seq"),
  N = c("105 proteome pairs; 89 phospho pairs","40 pairs","99 tumors","11 patients","221 pairs","11 pairs"),
  Platform = c("TMT11","LFQ iBAQ (imputed)","TMT/MS","GeoMx WTA","RNA-seq TPM","RNA-seq"),
  Role = c("Discovery","Independent longitudinal validation","Cross-sectional concordance","Spatial localization","Transcript-level boundary","Transcript-level support")
)
write.csv(t1, "submission_CancerLetters/03_Tables/Table1_cohorts.csv", row.names = FALSE)
cat("Table 1 written\n")
