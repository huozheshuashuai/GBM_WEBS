# =============================================================================
# Slingshot Pseudotime Trajectory Analysis
# =============================================================================
# Method: PCA > UMAP > slingshot with region_type as starting clusters
# Starting clusters: "Core" region as putative origin
# =============================================================================

library(slingshot)
library(SingleCellExperiment)

run_slingshot <- function(obj, start_cluster = "Core") {
  
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(obj)
  
  # Use UMAP embeddings for trajectory
  reducedDims(sce)$UMAP <- Embeddings(obj, "umap")
  
  # Run slingshot with region_type as clusters
  sce <- slingshot(sce, 
                   clusterLabels = "Region_Type",
                   reducedDim = "UMAP",
                   start.clus = start_cluster,
                   allow.breaks = FALSE)
  
  # Extract pseudotime
  pseudotime <- slingPseudotime(sce)
  
  # Add back to Seurat
  obj$pseudotime <- pseudotime[, 1]  # First lineage
  
  # Fit GAM to test WEBS ~ pseudotime (per cell type)
  library(mgcv)
  
  cell_types <- unique(obj$cell_type)
  gam_results <- lapply(cell_types, function(ct) {
    cells <- WhichCells(obj, expression = cell_type == ct)
    if (length(cells) < 50) return(NULL)
    
    df <- data.frame(
      WEBS = obj$WEBS[cells, 1],
      pseudotime = obj$pseudotime[cells]
    )
    
    # GAM with cubic regression spline
    fit <- gam(WEBS ~ s(pseudotime, bs = "cr"), data = df)
    
    data.frame(
      cell_type = ct,
      p_value = summary(fit)$s.pv,
      r_sq = summary(fit)$r.sq,
      dev_expl = summary(fit)$dev.expl
    )
  })
  
  gam_results <- do.call(rbind, gam_results)
  print(gam_results)
  
  return(list(obj = obj, sce = sce, gam_results = gam_results))
}
