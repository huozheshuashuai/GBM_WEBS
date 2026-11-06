setwd("E:/KIMI/GBM_WEBS_Project")
lo <- readRDS("data/raw/LinkedOmics_TCGA/tcga_gbm_clinical.rds")
cat("cols:", paste(head(colnames(lo), 10), collapse = ", "), "\n")
idcol <- colnames(lo)[1]
cat("id col:", idcol, "head:", head(as.character(lo[[1]]), 3), "\n")
xena <- read.delim(gzfile("E:/KIMI/data/Xena_GDC/TCGA-GBM.star_tpm.tsv.gz"),
                   check.names = FALSE, nrows = 2)
gdc_ids <- substr(colnames(xena)[-1], 1, 12)
lo_ids <- unique(substr(gsub("[.]", "-", as.character(lo[[1]])), 1, 12))
cat("LinkedOmics patients:", length(lo_ids), " GDC patients:", length(unique(gdc_ids)), "\n")
cat("overlap:", length(intersect(lo_ids, unique(gdc_ids))), "\n")
