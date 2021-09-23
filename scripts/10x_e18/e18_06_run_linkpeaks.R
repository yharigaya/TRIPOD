# load packages
library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)

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
file <- "e18.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path) ; rm(path)

file <- "genes.rds"
path <- file.path(dir.out, file)
genes <- readRDS(path) ; rm(path)

# run linkpeaks
DefaultAssay(e18) <- "ATAC"
e18 <- RegionStats(
	object = e18, 
	genome = BSgenome.Mmusculus.UCSC.mm10)
e18 <- LinkPeaks(
	object = e18, 
	peak.assay = "ATAC", 
	expression.assay = "SCT",
  pvalue_cutoff = 1, 
  score_cutoff = 0,
  distance = 2e+05, 
  method = "pearson", 
  genes.use = genes)
links <- Links(e18)
# correct p-values for two-sided tests
links$pvalue <- 2*pnorm(abs(links$zscore), lower.tail = FALSE)
file <- "links.rds"
path <- file.path(dir.out, file)
saveRDS(links, path); rm(path)

# set FDR < 0.01
fdr.thresh <- 0.01
links$adj <- p.adjust(links$pvalue, method = "BH")
links.pos <- links[links$score > 0 & links$adj < fdr.thresh]
links.neg <- links[links$score < 0 & links$adj < fdr.thresh]

# create data frames (gene, peak, coef, pval, adj)
columns <- c("gene", "peak", "coef", "pval", "adj")

pos.df <- as.data.frame(mcols(links.pos))
pos.df <- pos.df[, c(2, 3, 1, 5, 6)]
colnames(pos.df) <- columns

neg.df <- as.data.frame(mcols(links.neg))
neg.df <- neg.df[, c(2, 3, 1, 5, 6)]
colnames(neg.df) <- columns

# save to csv files
file <- paste0("linkpeaks_pos_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(pos.df, path); rm(path)

file <- paste0("linkpeaks_neg_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(neg.df, path); rm(path)

# # get sets of peak-gene pair strings
# linkpeaks.pos.strings <- unique(apply(pos.df, 1, convertPairToString, col.1 = 1, col.2 = 2))
# linkpeaks.neg.strings <- unique(apply(neg.df, 1, convertPairToString, col.1 = 1, col.2 = 2))
# c(length(linkpeaks.pos.strings), length(linkpeaks.neg.strings))
# # [1] 900   5
# 
# file <- "linkpeaks.pos.0.01.strings.rds"
# path <- file.path(dir.model, file)
# saveRDS(linkpeaks.pos.strings, path)
# 
# file <- "linkpeaks.neg.0.01.strings.rds"
# path <- file.path(dir.model, file)
# saveRDS(linkpeaks.neg.strings, path)
# 
# ## set FDR < 0.1
# fdr.thresh <- 0.1
# links.bh <- benjaminiHochsbergVectorAdjust(links$pvalue.2, fdr.thresh = fdr.thresh)
# table(links.bh$which.reject)
# # 
# # FALSE  TRUE 
# # 25794  1766
# identical(links.bh$pval.c.adj < fdr.thresh, links.bh$which.reject)
# # [1] TRUE
# links$adj <- links.bh$pval.c.adj
# links.pos <- links[links$score > 0 & links.bh$which.reject]
# links.neg <- links[links$score < 0 & links.bh$which.reject]
# 
# ## output csv files
# # gene, peak, coef, pval, adj
# pos.df <- as.data.frame(mcols(links.pos))
# pos.df <- pos.df[, c(2, 3, 1, 6, 7)]
# colnames(pos.df) <- c("gene", "peak", "coef", "pval", "adj")
# 
# neg.df <- as.data.frame(mcols(links.neg))
# neg.df <- neg.df[, c(2, 3, 1, 6, 7)]
# colnames(neg.df) <- c("gene", "peak", "coef", "pval", "adj")
# 
# # get sets of peak-gene pair strings
# linkpeaks.pos.strings <- unique(apply(pos.df, 1, convertPairToString, col.1 = 1, col.2 = 2))
# linkpeaks.neg.strings <- unique(apply(neg.df, 1, convertPairToString, col.1 = 1, col.2 = 2))
# c(length(linkpeaks.pos.strings), length(linkpeaks.neg.strings))
# # [1] 1728   38
# 
# # save to files
# file <- "linkpeaks.pos.0.1.csv"
# path <- file.path(dir.model, file)
# write.table(pos.df, path, quote = F, sep = ",")
# # pos.df <- read.csv(path)
# 
# file <- "linkpeaks.neg.0.1.csv"
# path <- file.path(dir.model, file)
# write.table(neg.df, path, quote = F, sep = ",")
# # neg.df <- read.csv(path)
# 
# file <- "linkpeaks.pos.0.1.strings.rds"
# path <- file.path(dir.model, file)
# saveRDS(linkpeaks.pos.strings, path)
# 
# file <- "linkpeaks.neg.0.1.strings.rds"
# path <- file.path(dir.model, file)
# saveRDS(linkpeaks.neg.strings, path)

