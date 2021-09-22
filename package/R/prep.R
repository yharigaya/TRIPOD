#' Get a metacell matrix
#'
#' This function takes a Seurat object as an input and returns a matrix
#' containing normalized RNA expression or chromatin accessibility per
#' metacell.
#'
#' @param object a Seurat object.
#' @param cluster.name a character string specifying the cell clusters on which
#' metacells are based. This must be one of the column names in the meta data
#' of the Seurat object.
#' @param assay a character string specifying the assay. This should be
#' "RNA" and "ATAC" for RNA expression and chromatin accessibility, respectively.
#'
#' @return a matrix
#' @import Seurat Matrix GenomicRanges
#' @export
getMetacellMatrix <- function(
	object,
  cluster.name = "seurat_clusters",
	assay) {
	if (all(assay != c("RNA", "ATAC"))) {
    stop('The assay argument should be either "RNA" or "ATAC".')
	}
	if (!(cluster.name %in% colnames(object@meta.data))) {
    stop("The cluster.name argument does not match the meta data.")
	}
  clust.levels <- levels(object@meta.data[[cluster.name]])
  assay.matrix <- GetAssay(object = object, assay = assay)
  nrow <- length(clust.levels)
  ncol <- nrow(assay.matrix)
  metacell <- matrix(nrow = nrow,
    ncol = ncol)
  rownames(metacell) <- paste0("metacell_", clust.levels)
  colnames(metacell) <- rownames(assay.matrix)
  for (i in 1:nrow(metacell)) {
    metacell[i, ] <- apply(assay.matrix[, clust.levels == (i - 1)], 1, sum)
    # adjust for total read counts for each metacell
    metacell[i, ] <- metacell[i, ]/sum(metacell[i, ])*10^6
  }
  return(metacell)
}

#' Optimize cluster resolutions
#'
#' This function takes a Seurat object as input and returns a data frame
#' containing cluster numbers at specified resolutions.
#'
#' @param object a Seurat object.
#' @param assay.name a character string specifying the assay name. This must
#' be one of "SCT", "ATAC", or "WNN".
#' @param graph.name a character string specifying the graph name. This must be
#' one of the graphs stored in the Seurat object.
#' @param algorithm an integer representing a clustering algorithm. See
#' {\code{\link{FindClusters}}} in the Seurat package.
#' @param resolutions a numerical vector specifying resolutions to be examined.
#' @param min.num an integer setting a threshold for the minimum number of
#' cells per cluster.
#' @param ... further arguments to be passed to {\code{\link{FindClusters}}}.
#'
#' @return a data frame with three columns:
#' \describe{
#'   \item{resolusion}{the resolution examined}
#'   \item{num_clusters}{the number of clusters}
#'   \item{num_below}{the number of clusters with fewer than
#'   the threshold number of single cells}
#' }
#' @import Seurat
#' @export
optimizeResolution <- function(
	object, assay.name, graph.name = NULL,
	algorithm = NULL, resolutions = seq(10, 35, 5), min.num, ...
) {
	if (all(assay.name != c("SCT", "ATAC", "WNN"))) {
    stop('The assay.name argument must be one of "SCT", "ATAC", or "WNN".')
	}
  # object <- FindNeighbors(object, dims = 1:30)
  num.clusters <- num.below <- rep(NA, length(resolutions))
  for (i in 1:length(resolutions)){
    res <- resolutions[i]
    if (assay.name %in% c("SCT", "ATAC")) {
    	# DefaultAssay(object) <- assay.name
    	if (is.null(graph.name)) graph.name <- paste0(assay.name, "_snn")
      if (is.null(algorithm)) algorithm <- 1
      # snn.name <- paste0(assay.name, "_snn_res.", res)
    } else if (assay.name == "WNN") {
    	if (is.null(graph.name)) graph.name <- "wsnn"
      if (is.null(algorithm)) algorithm <- 3
      # snn.name <- paste0("wsnn_res.", res)
    }
    object <- FindClusters(object, graph.name = graph.name,
    	resolution = res, algorithm = algorithm, verbose = FALSE, ...)
    meta.col.name <- paste0(graph.name, "_res.", res)
    num.clusters[i] <- length(levels(object@meta.data[[meta.col.name]]))
    num.below[i] <- sum(table(object@meta.data[[meta.col.name]]) < min.num)
  }
  results <- data.frame(resolution = resolutions,
                        num_clusters = num.clusters,
                        num_below = num.below)
  return(results)
}

#' Get cell types for metacells
#'
#' @param object a Seurat object.
#' @param celltype.col.name a character string specifying the name of the
#' meta data column in the Seurat object containing cell types.
#' @param cluster.col.name a character string specifying the name of the
#' meta data column in the Seurat object containing cluster numbers.
#'
#' @return a character vector
#' @export
getCellTypeForMetacell <- function(
	object, metacell.matrix,
	celltype.col.name = "celltype", cluster.col.name = "seurat_clusters"
) {
  tmp <- table(unlist(object[[celltype.col.name]]),
  	unlist(object[[cluster.col.name]]))
  metacell.celltype <- rep(NA, ncol(tmp))
  for (i in 1:length(metacell.celltype)) {
    tmp.i_1 <- tmp[, colnames(tmp) == as.character(i-1)]
    metacell.celltype[i] <- names(tmp.i_1)[which.max(tmp.i_1)]
  }
  return(metacell.celltype)
}

#' Get a mapping between cell types and colors for single cells
#'
#' @param object a Seurat object.
#' @param reduction a character string specifying a dimension reduction method.
#' @param celltype.col.name a character string specifying the name of the meta
#' data column in the Seurat object containing cell types.
#'
#' @return a data frame
#' @export
getColorsForSingleCells <- function(
	object, reduction, celltype.col.name
) {
	# get the corresponding color for each cell type from Seurat
	p <- Seurat::DimPlot(e18, reduction = reduction, label = TRUE,
		group.by = celltype.col.name)
  # use ggplot_build to deconstruct the ggplot object
  pbuild <- ggplot2::ggplot_build(p)
  # get the color palette by Seurat
  pdata <- pbuild$data[[1]]
  pdata <- cbind(object[[celltype.col.name]], pdata)
  sc.color.map <- pdata[, 1:2]
  colnames(sc.color.map) <- c("celltype", "color")
  return(sc.color.map)
}

#' Get a mapping between cell types and colors for meta cells
#'
#' @param metacell.celltype a character string
#' @param sc.color.map a data frame
#'
#' @return a data frame
#' @export
getColorsForMetacells <- function(metacell.celltype, sc.color.map) {
  metacell.celltype.col <- rep(NA, length(metacell.celltype))
  for (i in 1:length(metacell.celltype)) {
    metacell.celltype.col[i] <-
    	sc.color.map$color[min(which(sc.color.map$celltype == metacell.celltype[i]))]
  }
  metacell.color.map <- data.frame(
  	celltype = metacell.celltype,
	  color = metacell.celltype.col
  )
}

#' Get a mapping between cell types and colors
#'
#' @param ordered.celltypes a character string
#' @param metacell.color.map a data frame
#'
#' @return a data frame
#' @export
getColors <- function(ordered.celltypes, metacell.color.map) {
  color.map <- unique(metacell.color.map)
  color.map <- color.map[match(ordered.celltypes, color.map$celltype), ]
  rownames(color.map) <- NULL
  return(color.map)
}

