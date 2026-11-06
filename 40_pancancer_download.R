# Download all 33 TCGA cohorts via TCGAbiolinks (pan-cancer WEBS analysis)
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({
  library(TCGAbiolinks); library(SummarizedExperiment)
})

CANCERS <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD",
             "TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC",
             "TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC",
             "TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ",
             "TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA","TCGA-THYM",
             "TCGA-UCEC","TCGA-UCS","TCGA-UVM")

prepare_cohort <- function(project) {
  cache <- file.path("data/raw", project, paste0(project, "_prepared.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  dir.create(file.path("data/raw", project), recursive = TRUE, showWarnings = FALSE)
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  tryCatch({
    GDCdownload(query, directory = file.path("data/raw", project))
    se <- GDCprepare(query, directory = file.path("data/raw", project))
    st <- substr(colnames(se), 14, 15)
    keep_types <- if (project == "TCGA-SKCM") c("01", "06") else "01"
    se <- se[, st %in% keep_types]
    tpm <- assay(se, "tpm_unstrand")
    symbols <- rowData(se)$gene_name
    keep <- !is.na(symbols) & symbols != ""
    tpm <- tpm[keep, ]; symbols <- symbols[keep]
    o <- order(rowMeans(tpm), decreasing = TRUE)
    tpm <- tpm[o, ]; symbols <- symbols[o]
    tpm <- tpm[!duplicated(symbols), ]; rownames(tpm) <- symbols[!duplicated(symbols)]
    clin <- as.data.frame(colData(se))
    os_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death, clin$days_to_last_follow_up)
    os_event <- ifelse(clin$vital_status == "Dead", 1, 0)
    clinical <- data.frame(sample = colnames(tpm),
                           OS.time = as.numeric(os_time) / 30.44, OS = os_event)
    clinical <- clinical[!is.na(clinical$OS.time) & clinical$OS.time > 0, ]
    tpm <- tpm[, clinical$sample]
    saveRDS(list(tpm = tpm, clinical = clinical), cache)
    message(project, " OK: ", ncol(tpm), " samples")
    list(tpm = tpm, clinical = clinical)
  }, error = function(e) { message(project, " FAILED: ", e$message); NULL })
}

for (proj in CANCERS) prepare_cohort(proj)
cat("PAN-CANCER DOWNLOAD COMPLETE\n")
