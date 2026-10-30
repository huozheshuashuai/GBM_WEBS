# =============================================================================
# snATAC-seq Analysis (GSE138794, 6 samples)
# =============================================================================
# Methods: Signac -> quality filtering -> gene activity matrix -> 
#          LSI dimensionality reduction -> UMAP -> cisTopic for TF topic modeling
#          DiffBind for differential chromatin accessibility
# =============================================================================

library(Signac)
library(Seurat)
library(cisTopic)

# ---- 1. Create Chromatin Assay ----
# Requires: fragment file (.tsv.gz), peak matrix, cell metadata

create_snatac_object <- function(counts_matrix, metadata, fragments_path, genome = "hg38") {
  
  # Create chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = counts_matrix,
    sep = c(":", "-"),
    genome = genome,
    fragments = fragments_path,
    min.cells = 10,
    min.features = 200
  )
  
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  return(obj)
}

# ---- 2. QC ----
run_atac_qc <- function(obj) {
  obj <- NucleosomeSignal(object = obj)
  obj <- TSSEnrichment(object = obj, fast = FALSE)
  
  # Filter based on QC metrics
  obj <- subset(
    x = obj,
    subset = nCount_peaks > 3000 &
      nCount_peaks < 20000 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )
  
  return(obj)
}

# ---- 3. Gene Activity Matrix ----
compute_gene_activity <- function(obj, annotation) {
  # annotation: EnsDb or TxDb object
  # Get gene annotations
  gene.activities <- GeneActivity(obj, annotation = annotation)
  
  # Add as new assay
  obj[["RNA"]] <- CreateAssayObject(counts = gene.activities)
  obj <- NormalizeData(obj, assay = "RNA")
  
  return(obj)
}

# ---- 4. Dimensionality Reduction ----
run_atac_dimred <- function(obj) {
  obj <- FindTopFeatures(obj, min.cutoff = "q0")
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  
  # Use LSI components (exclude first component = sequencing depth)
  obj <- RunUMAP(object = obj, reduction = "lsi", dims = 2:30)
  obj <- FindNeighbors(object = obj, reduction = "lsi", dims = 2:30)
  obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
  
  return(obj)
}

# ---- 5. Transfer WEBS from scRNA-seq ----
transfer_webs_atac <- function(atac_obj, rna_obj) {
  # Find anchors
  transfer.anchors <- FindTransferAnchors(
    reference = rna_obj,
    query = atac_obj,
    reference.assay = "RNA",
    query.assay = "RNA",  # gene activity
    reduction = "cca",
    dims = 1:30
  )
  
  # Transfer WEBS
  predicted <- TransferData(
    anchorset = transfer.anchors,
    refdata = rna_obj$WEBS,
    dims = 1:30
  )
  
  atac_obj <- AddMetaData(atac_obj, metadata = predicted)
  atac_obj$predicted_WEBS <- predicted$predicted.id
  
  # Also classify
  q1 <- quantile(atac_obj$predicted_WEBS, 0.25, na.rm = TRUE)
  q3 <- quantile(atac_obj$predicted_WEBS, 0.75, na.rm = TRUE)
  atac_obj$WEBS_group <- ifelse(
    atac_obj$predicted_WEBS <= q1, "LLS",
    ifelse(atac_obj$predicted_WEBS >= q3, "LHS", "LDTS")
  )
  
  return(atac_obj)
}

# ---- 6. cisTopic TF Analysis ----
run_cisTopic <- function(obj, n_topics = 30) {
  # cisTopic requires specific input format
  # Extract count matrix from chromatin assay
  counts <- GetAssayData(obj, assay = "peaks", layer = "counts")
  
  # Run cisTopic
  cisTopicObject <- createcisTopicObject(counts, project = "GBM_ATAC")
  cisTopicObject <- runModels(cisTopicObject, topic = c(10, 20, 30, 40, 50), 
                               seed = 123, nCores = 4)
  cisTopicObject <- selectModel(cisTopicObject)
  
  # Get cell-topic distributions
  cellTopics <- cisTopicObject@cell.classification
  
  return(cisTopicObject)
}

# ---- 7. Differential Accessibility ----
find_diff_peaks <- function(obj, group.by = "WEBS_group", ident.1 = "LHS", ident.2 = "LLS") {
  da_peaks <- FindMarkers(
    object = obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min.pct = 0.05,
    test.use = "LR",
    latent.vars = "peak_region_fragments"
  )
  
  return(da_peaks)
}

# ---- 8. Motif Enrichment ----
run_motif_analysis <- function(obj, da_peaks, pvalueCutoff = 0.05) {
  library(chromVAR)
  library(motifmatchr)
  library(JASPAR2020)
  
  # Get PWMs from JASPAR
  pwm_list <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM")
  )
  
  # Add motif information
  obj <- AddMotifs(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm_list)
  
  # Test motif enrichment
  enriched.motifs <- FindMotifs(
    object = obj,
    features = rownames(da_peaks)[da_peaks$p_val_adj < pvalueCutoff]
  )
  
  return(enriched.motifs)
}
