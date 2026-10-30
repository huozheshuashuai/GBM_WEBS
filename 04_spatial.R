# =============================================================================
# Spatial Transcriptomics Analysis (10x Visium)
# =============================================================================
# Cohorts: GSE194329 (6 slices), 5 patients
# Methods: CreateSeuratObject(counts, images) -> SCTransform -> cluster
#          Reference-mapping via scRNA-seq atlas -> RCTD for decomposition
# =============================================================================

library(Seurat)
library(SeuratData)

# ---- 1. Load Visium Data ----
load_visium <- function(spatial_dir, sample_id) {
  # spatial_dir should contain: filtered_feature_bc_matrix/ + spatial/
  obj <- Load10X_Spatial(
    data.dir = spatial_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample_id,
    filter.matrix = TRUE
  )
  return(obj)
}

# ---- 2. Quality Control & Preprocessing ----
preprocess_spatial <- function(obj) {
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
  obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
  obj <- FindClusters(obj, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
  return(obj)
}

# ---- 3. Reference Mapping (RCTD / Seurat Transfer) ----
# Paper uses RCTD for cell-type decomposition
# Alternative: Seurat's TransferData (reference = scRNA-seq atlas)

deconvolute_rctd <- function(spatial_obj, sc_reference) {
  # RCTD is available from CRAN: install.packages("spacexr")
  # But requires careful formatting
  # Here is the Seurat alternative:
  
  # Find anchors between spatial and scRNA-seq
  anchors <- FindTransferAnchors(
    reference = sc_reference,
    query = spatial_obj,
    normalization.method = "SCT",
    dims = 1:30
  )
  
  # Transfer cell type labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = sc_reference$cell_type,
    dims = 1:30
  )
  
  spatial_obj <- AddMetaData(spatial_obj, metadata = predictions)
  
  return(spatial_obj)
}

# ---- 4. WEBS on Spatial ----
compute_webs_spatial <- function(spatial_obj, sc_reference) {
  # Predict WEBS from spatial data using reference mapping
  anchors <- FindTransferAnchors(
    reference = sc_reference,
    query = spatial_obj,
    normalization.method = "SCT",
    dims = 1:30
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = sc_reference$WEBS,
    dims = 1:30
  )
  
  spatial_obj$predicted_WEBS <- predictions$predicted.id
  
  # Also compute WEBS from expression (approximate)
  source("tools/WEBS/webs.R")
  spatial_obj <- webs_compute(spatial_obj, assay = "SCT", slot = "data")
  
  return(spatial_obj)
}

# ---- 5. Cross-Slice Analysis ----
# Compare WEBS trajectories across slices
# Identify high-WEBS hotspots (spatial clusters)

analyze_spatial_patterns <- function(spatial_obj_list) {
  library(ggplot2)
  
  results <- lapply(names(spatial_obj_list), function(id) {
    obj <- spatial_obj_list[[id]]
    
    data.frame(
      slice = id,
      mean_webs = mean(obj$WEBS, na.rm = TRUE),
      max_webs = max(obj$WEBS, na.rm = TRUE),
      n_hotspots = sum(obj$WEBS > quantile(obj$WEBS, 0.9, na.rm = TRUE))
    )
  })
  
  do.call(rbind, results)
}
