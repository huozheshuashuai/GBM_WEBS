# Trace: why do 29 (RDS cache) and 61 (raw zip) give different six-gene results?
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(survival) })
PAPER_COEFS <- c(HIF1A = -0.284, SIRT2 = -0.176, MMP2 = 0.262,
                 ITGA5 = 0.230, HDAC1 = -0.220, SPP1 = 0.146)

# Path A: RDS cache (used by 29)
eA <- readRDS("data/raw/CGGA/cgga325_expr.rds")
cA <- readRDS("data/raw/CGGA/cgga325_clinical.rds")
cat("RDS expr: dim", dim(eA), "range", round(range(as.matrix(eA[, 1:5]), na.rm = TRUE), 2), "\n")
cat("RDS rownames head:", head(rownames(eA), 5), "\n")
cat("RDS sample col head:", head(cA$sample, 3), "\n")

# Path B: raw txt (used by 61)
ex <- read.delim("data/raw/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
                 check.names = FALSE, nrows = 5)
cat("raw cols head:", colnames(ex)[1:4], "\n")
cat("raw first col head:", head(ex[[1]], 5), "\n")

# full raw
ex <- read.delim("data/raw/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",
                 check.names = FALSE)
cat("raw dim:", dim(ex), "\n")
emB <- as.matrix(sapply(ex[, -1], as.numeric))
rownames(emB) <- sub("\\..*$", "", ex[[1]])
cat("rownames emB head:", head(rownames(emB), 5), "\n")
cat("genes matched to PAPER_COEFS in B:", intersect(names(PAPER_COEFS), rownames(emB)), "\n")
cat("genes matched in A:", intersect(names(PAPER_COEFS), rownames(eA)), "\n")

# risk scores both ways (same samples)
g <- intersect(names(PAPER_COEFS), rownames(eA))
rsA <- as.vector(t(as.matrix(eA[g, ])) %*% PAPER_COEFS[g])
names(rsA) <- colnames(eA)

emB2 <- log2(emB + 1)
gB <- intersect(names(PAPER_COEFS), rownames(emB2))
rsB <- as.vector(t(emB2[gB, ]) %*% PAPER_COEFS[gB])
names(rsB) <- colnames(emB2)

# match samples
idA <- names(rsA)
idB <- gsub("\\..*$", "", names(rsB))
common <- intersect(idA, idB)
cat("common samples:", length(common), "\n")
rsB_matched <- rsB[match(common, idB)]
rsA_matched <- rsA[match(common, idA)]
cat("cor(rsA, rsB):", cor(rsA_matched, rsB_matched, use = "complete.obs"), "\n")

# Cox on each
clinA <- cA[match(common, cA$sample), ]
coxA <- summary(coxph(Surv(OS.time, OS) ~ rsA_matched, data = clinA))
cat("Path A (RDS): Cox p =", coxA$coefficients[1, 5], "\n")
# for path B use same clinical
coxB <- summary(coxph(Surv(clinA$OS.time, clinA$OS) ~ rsB_matched))
cat("Path B (raw+log2): Cox p =", coxB$coefficients[1, 5], "\n")
