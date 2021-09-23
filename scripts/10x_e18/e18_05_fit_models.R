# load packages
library(GenomicRanges)
library(nbpMatching)
library(BiocParallel)

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
file <- "transcripts.gr.rds"
path <- file.path(dir.out, file)
transcripts.gr <- readRDS(path); rm(path)

file <- "peaks.gr.rds"
path <- file.path(dir.out, file)
peak.gr <- readRDS(path); rm(path)

file <- "motifxTF.rds"
path <- file.path(dir.out, file)
motifxTF <- readRDS(path); rm(path)

file <- "peakxmotif.rds"
path <- file.path(dir.out, file)
peakxmotif <- readRDS(path); rm(path)

file <- "metacell.rna.rds"
path <- file.path(dir.out, file)
metacell.rna <- readRDS(path); rm(path)

file <- "metacell.peak.rds"
path <- file.path(dir.out, file)
metacell.peak <- readRDS(path); rm(path)

file <- "hvg.rds"
path <- file.path(dir.out, file)
hvg <- readRDS(path); rm(path)

file <- "metacell.color.map.rds"
path <- file.path(dir.out, file)
metacell.color.map <- readRDS(path); rm(path)

# get a list of target genes of interest
tfs.neuro <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Tbr1")
tfs.glio <- c("Olig2", "Sox10", "Nkx2-2", "Sox9", "Nfia", "Ascl1")
genes <- unique(c(hvg[1:1000], tfs.neuro, tfs.glio))
file <- "genes.rds"
path <- file.path(dir.out, file)
saveRDS(genes, path) ; rm(path)

ext.upstream <- ext.downstream <- 2e5

xymats.list <- bplapply(
	genes,
	getXYMatrices,
  ext.upstream = ext.upstream,
  transcripts.gr = transcripts.gr,
  peaks.gr = peaks.gr,
  metacell.rna = metacell.rna,
  metacell.peak = metacell.peak,
  peakxmotif = peakxmotif,
  motifxTF = motifxTF,
  metacell.celltype = metacell.color.map$celltype,
  metacell.celltype.col = metacell.color.map$color
)
names(xymats.list) <- genes
file <- "xymats.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.list, path); rm(path)

# fit marginal models
xymats.m.list <- bplapply(
	xymats.list,
  fitModel,
	model.name = "marginal"
)
names(xymats.m.list) <- genes
file <- "xymats.m.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.m.list, path); rm(path)

# fit conditional models
xymats.c.list <- bplapply(
	xymats.list,
  fitModel,
	model.name = "conditional"
)
names(xymats.c.list) <- genes
file <- "xymats.c.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.c.list, path); rm(path)

# fit interaction models
xymats.i.list <- bplapply(
	xymats.list,
  fitModel,
	model.name = "interaction"
)
names(xymats.i.list) <- genes
file <- "xymats.i.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.i.list, path); rm(path)

# run TRIPOD matching Xt
xymats.tripod.Xt.list <- bplapply(
	xymats.list,
	fitModel,
	model.name = "TRIPOD",
	match.by = "Xt"
)
names(xymats.tripod.Xt.list) <- genes
file <- "xymats.tripod.Xt.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.tripod.Xt.list, path); rm(path)

# run TRIPOD matching Yj
xymats.tripod.Yj.list <- bplapply(
	xymats.list,
	fitModel,
	model.name = "TRIPOD",
	match.by = "Yj"
)
names(xymats.tripod.Yj.list) <- genes
file <- "xymats.tripod.Yj.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.tripod.Yj.list, path); rm(path)
