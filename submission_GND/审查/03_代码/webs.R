# =============================================================================
# WEBS (Writer-Eraser Balance Score) - Standalone Implementation
# =============================================================================
# Based on: Chen et al. "Spliceosome activation defines an epigenetically 
#           hyperactive and therapeutically targetable malignant cell state 
#           in glioblastoma"
#
# Methods: WEBS = UCell(writer genes) - UCell(eraser genes)
#          per-cell computation using rank-based signature scores
# =============================================================================

library(Seurat)
library(UCell)

#' Compute Writer-Eraser Balance Score (WEBS)
#'
#' @param object Seurat object or expression matrix (genes x cells)
#' @param writer_genes Vector of writer gene symbols
#' @param eraser_genes Vector of eraser gene symbols
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "data" for log-normalized)
#' @param name Name for the resulting metadata column
#' @return Seurat object with WEBS stored in metadata, or vector if matrix input
#'
#' @examples
#' # Writers: HATs/KATs with lactyltransferase activity
#' writers <- c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP")
#' # Erasers: HDACs and sirtuins with delactylase activity
#' erasers <- c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2", 
#'              "SIRT3", "SIRT6", "SIRT7")
#' seurat_obj <- webs_compute(seurat_obj, writers, erasers)
#'
webs_compute <- function(object, 
                         writer_genes = c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP"),
                         eraser_genes = c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2",
                                          "SIRT3", "SIRT6", "SIRT7"),
                         assay = "RNA",
                         slot = "data",
                         name = "WEBS") {
  
  if (inherits(object, "Seurat")) {
    # Check which genes are present
    all_genes <- rownames(object[[assay]])
    writers_present <- intersect(writer_genes, all_genes)
    erasers_present <- intersect(eraser_genes, all_genes)
    
    message("Writer genes present: ", length(writers_present), "/", length(writer_genes))
    if (length(setdiff(writer_genes, writers_present)) > 0) {
      message("Missing writers: ", paste(setdiff(writer_genes, writers_present), collapse = ", "))
    }
    
    message("Eraser genes present: ", length(erasers_present), "/", length(eraser_genes))
    if (length(setdiff(eraser_genes, erasers_present)) > 0) {
      message("Missing erasers: ", paste(setdiff(eraser_genes, erasers_present), collapse = ", "))
    }
    
    # Compute UCell scores
    object <- AddModuleScore_UCell(object, 
                                    features = list(Writers = writers_present),
                                    assay = assay, slot = slot,
                                    name = "_WriterUCell")
    object <- AddModuleScore_UCell(object,
                                    features = list(Erasers = erasers_present),
                                    assay = assay, slot = slot,
                                    name = "_EraserUCell")
    
    # WEBS = Writer UCell - Eraser UCell (stored as plain numeric vector)
    object[[name]] <- as.numeric(object[["Writers_WriterUCell"]][, 1] -
                                 object[["Erasers_EraserUCell"]][, 1])
    
    # Store component scores
    object[["Writer_UCell"]] <- object[["Writers_WriterUCell"]]
    object[["Eraser_UCell"]] <- object[["Erasers_EraserUCell"]]
    
    # Clean up intermediate columns
    object[["Writers_WriterUCell"]] <- NULL
    object[["Erasers_EraserUCell"]] <- NULL
    
    return(object)
    
  } else {
    stop("Input must be a Seurat object. For matrices, use webs_compute_matrix().")
  }
}

#' Classify cells into WEBS states
#'
#' @param object Seurat object with WEBS computed
#' @param webs_col Column name containing WEBS values
#' @param method Stratification method: "iqr" (default), "tertile", "median"
#' @return Seurat object with WEBS_group metadata
#'
webs_classify <- function(object, webs_col = "WEBS", method = "iqr") {
  webs <- object[[webs_col]][, 1]
  
  if (method == "iqr") {
    # Default: LLS (bottom 25%), LDTS (middle 50%), LHS (top 25%)
    q1 <- quantile(webs, 0.25)
    q3 <- quantile(webs, 0.75)
    
    group <- ifelse(webs <= q1, "LLS",
                    ifelse(webs >= q3, "LHS", "LDTS"))
    group <- factor(group, levels = c("LLS", "LDTS", "LHS"))
    
  } else if (method == "tertile") {
    t1 <- quantile(webs, 1/3)
    t2 <- quantile(webs, 2/3)
    group <- ifelse(webs <= t1, "LLS",
                    ifelse(webs >= t2, "LHS", "LDTS"))
    group <- factor(group, levels = c("LLS", "LDTS", "LHS"))
    
  } else if (method == "median") {
    med <- median(webs)
    group <- ifelse(webs < med, "LLS", "LHS")
    group <- factor(group, levels = c("LLS", "LHS"))
  }
  
  object[["WEBS_group"]] <- group
  return(object)
}

#' Plot WEBS distributions
#'
#' @param object Seurat object
#' @param group.by Grouping variable
#' @return ggplot object
webs_plot <- function(object, group.by = "cell_type") {
  library(ggplot2)
  library(ggpubr)
  
  df <- data.frame(
    WEBS = object[["WEBS"]][, 1],
    Group = object[[group.by]][, 1]
  )
  
  p <- ggplot(df, aes(x = Group, y = WEBS, fill = Group)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    stat_compare_means(method = "kruskal.test") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "WEBS Distribution", x = group.by, y = "WEBS")
  
  return(p)
}

#' Compute WEBS for bulk expression matrix (TPM/FPKM)
#'
#' @param expr_matrix Gene expression matrix (genes x samples)
#' @param writer_genes Writer gene symbols
#' @param eraser_genes Eraser gene symbols
#' @return Data frame with WEBS, WriterScore, EraserScore per sample
webs_compute_bulk <- function(expr_matrix,
                               writer_genes = c("EP300", "KAT2A", "KAT2B", "KAT5", "KAT8", "CREBBP"),
                               eraser_genes = c("HDAC1", "HDAC2", "HDAC3", "SIRT1", "SIRT2",
                                                "SIRT3", "SIRT6", "SIRT7")) {
  
  all_genes <- rownames(expr_matrix)
  writers_present <- intersect(writer_genes, all_genes)
  erasers_present <- intersect(eraser_genes, all_genes)
  
  # For bulk: use mean log2 expression (consistent with Methods description)
  writer_score <- colMeans(log2(expr_matrix[writers_present, ] + 1), na.rm = TRUE)
  eraser_score <- colMeans(log2(expr_matrix[erasers_present, ] + 1), na.rm = TRUE)
  
  webs <- writer_score - eraser_score
  
  data.frame(
    sample = colnames(expr_matrix),
    WEBS = webs,
    WriterScore = writer_score,
    EraserScore = eraser_score
  )
}

#' Protein-level WEBS computation (CPTAC)
#'
#' @param protein_matrix Protein abundance matrix (proteins x samples)
#' @param writer_proteins Protein IDs for writers
#' @param eraser_proteins Protein IDs for erasers
#' @return Data frame with protein-level WEBS
webs_compute_protein <- function(protein_matrix,
                                  writer_proteins = c("HAT1", "KAT2A", "KAT2B", "KAT5", "KAT6A",
                                                      "KAT7", "KAT8", "CREBBP", "EP300", "KAT14"),
                                  eraser_proteins = c("HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5",
                                                      "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10",
                                                      "HDAC11", "SIRT1", "SIRT2", "SIRT3", "SIRT5", "SIRT6")) {
  
  all_proteins <- rownames(protein_matrix)
  writers_present <- intersect(writer_proteins, all_proteins)
  erasers_present <- intersect(eraser_proteins, all_proteins)
  
  writer_score <- colMeans(protein_matrix[writers_present, ], na.rm = TRUE)
  eraser_score <- colMeans(protein_matrix[erasers_present, ], na.rm = TRUE)
  
  webs <- writer_score - eraser_score
  
  data.frame(
    sample = colnames(protein_matrix),
    WEBS = webs,
    WriterScore = writer_score,
    EraserScore = eraser_score
  )
}
