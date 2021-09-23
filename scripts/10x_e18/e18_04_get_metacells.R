# load packages
library(Seurat)
library(Signac)

# set directories
dir.in <- "data"
dir.out <- "output"
dir.fig <- "figures"
dir.r <- "functions"

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "e18.intersect.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path); rm(path)

# get rna meta-cells by clustering cells on basis of their scRNA-seq profiles
DefaultAssay(e18) <- "SCT"
e18 <- FindNeighbors(e18, reduction = "pca", dims = 1:30)

# optimize clustering resolution
num.clusters <- optimizeResolution(
	object = e18,
	graph.name = "SCT_snn",
	assay.name = "SCT",
	resolutions = seq(10, 35, 5),
	min.num = 20
)
# write the result to a csv file
file <- "num_clusters.csv"
path <- file.path(dir.out, file)
write.csv(num.clusters, path); rm(path)

res <- 15
e18 <- FindClusters(e18, graph.name = "SCT_snn",
    	resolution = res, verbose = FALSE)

# reorder the factor levels
tmp <- as.character(e18$SCT_snn_res.15)
levels.15 <- as.character(sort(as.numeric(levels(e18$SCT_snn_res.15))))
tmp <- factor(tmp, levels = levels.15)
e18$SCT_snn_res.15 <- tmp
e18@meta.data$seurat_clusters <- e18@meta.data$SCT_snn_res.15
e18@meta.data$SCT_snn_res.15 <- NULL

# get metacell matrices
metacell.rna <- matrix(nrow = length(levels(e18$seurat_clusters)),
	ncol = nrow(e18@assays$RNA))
rownames(metacell.rna) <- paste("metacell_", levels(e18$seurat_clusters), sep = "")
colnames(metacell.rna) <- rownames(e18@assays$RNA)
for(i in 1:nrow(metacell.rna)){
  metacell.rna[i,] <- apply(e18@assays$RNA@counts[, e18$seurat_clusters == (i-1)], 1, sum)
  # adjust for library size for each metacell
  metacell.rna[i,] <- metacell.rna[i,]/sum(metacell.rna[i,])*10^6 
}

metacell.peak <- matrix(nrow = length(levels(e18$seurat_clusters)),
	ncol = nrow(e18@assays$ATAC))
rownames(metacell.peak) <- paste("metacell_", levels(e18$seurat_clusters), sep = "")
colnames(metacell.peak) <- rownames(e18@assays$ATAC)
for(i in 1:nrow(metacell.peak)){
  metacell.peak[i,] <- apply(e18@assays$ATAC@counts[, e18$seurat_clusters == (i-1)], 1, sum)
  # adjust for library size for each metacell
  metacell.peak[i,] <- metacell.peak[i,]/sum(metacell.peak[i,])*10^6
}

# remove metacell clusters with fewer than 20 cells
remove <- as.vector(table(e18$seurat_clusters)) < 20
names.remove <- names(table(e18$seurat_clusters))[remove]
metacell.rna <- metacell.rna[!remove, ]
metacell.peak <- metacell.peak[!remove, ]

file <- "metacell.rna.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.rna, path); rm(path)

file <- "metacell.peak.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.peak, path); rm(path)

# remove cells in the removed metacell clusters from the Seurat object
e18 <- e18[, !(e18$seurat_clusters %in% names.remove)]
e18$seurat_clusters <- droplevels(e18$seurat_clusters)

file <- "e18.metacell.rds"
path <- file.path(dir.out, file)
saveRDS(e18, path); rm(path)

# assigning cell types and colors for visualization
metacell.celltype <- getCellTypeForMetacell(e18,
	celltype.col.name = "celltype", cluster.col.name = "seurat_clusters")
sc.color.map <- getColorsForSingleCells(
	object = e18, reduction = "umap.rna", celltype.col.name = "celltype")
metacell.color.map <- getColorsForMetacells(
	metacell.celltype = metacell.celltype, sc.color.map = sc.color.map)
color.map <- getColors(ordered.celltypes = levels(e18$celltype),
	metacell.color.map = metacell.color.map)

file <- "sc.color.map.rds"
path <- file.path(dir.out, file)
saveRDS(sc.color.map, path); rm(path)

file <- "metacell.color.map.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.color.map, path); rm(path)

file <- "color.map.rds"
path <- file.path(dir.out, file)
saveRDS(color.map, path); rm(path)

# obtain top 3000 highly variable genes based on SCT
DefaultAssay(e18) <- "SCT"
head(VariableFeatures(e18))
length(VariableFeatures(e18))
hvg <- VariableFeatures(e18)

file <- "hvg.rds"
path <- file.path(dir.out, file)
saveRDS(hvg, path); rm(path)
