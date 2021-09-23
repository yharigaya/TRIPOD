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
file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "xymats1.list.rds"
path <- file.path(dir.out, file)
xymats1.list <- readRDS(path); rm(path)

file <- "xymats4.list.rds"
path <- file.path(dir.out, file)
xymats4.list <- readRDS(path); rm(path)

file <- "xymats5X.list.rds"
path <- file.path(dir.out, file)
xymats5X.list <- readRDS(path); rm(path)

file <- "xymats5Y.list.rds"
path <- file.path(dir.out, file)
xymats5Y.list <- readRDS(path); rm(path)

# set FDR < 0.01
fdr.thresh <- 0.01

# set the sign coefficient
sign <- "positive"; sign.str <- "pos"
# sign <- "negative"; sign.str <- "neg"

# get peak-gene pairs with significantly positive alpha coefficients
# from model 1 (marginal Yg ~ Xt)
xymats1.alpha.pos.df <- getPeakGenePairs(
	xymats.list = xymats.list.list$xymats1,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "1a")
file <- paste0("xymats1.alpha.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats1.alpha.pos.df, path, sep = ",", quote = FALSE)

# get TF-gene pairs with significantly positive beta coefficients
# from model 1 (marginal Yg ~ Yj)
xymats1.beta.pos.df <- getTFGenePairs(
	xymats.list = xymats.list.list$xymats1,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "1a")
file <- paste0("xymats1.beta.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats1.beta.pos.df, path, sep = ",", quote = FALSE)

# get trios with significantly positive gamma coefficients 
# from model 4 (interaction)
xymats4.gamma.pos.df <- getTrios(
	xymats.list = xymats.list.list$xymats4, 
	fdr.thresh = fdr.thresh, 
	sign = sign,
  model.name = "4a")
file <- paste0("xymats4.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats4.gamma.pos.df, path, sep = ",", quote = FALSE)

# get trios with significantly positive beta coefficients
# from model 5X level 1 (TRIPOD level 1 matching Xt)
xymats5X.beta.pos.df <- getTrios(
	xymats.list = xymats.list.list$xymats5X, 
	fdr.thresh = fdr.thresh, 
	sign = sign,
  model.name = "5a",
	level = 1)
file <- paste0("xymats5X.beta.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5X.beta.pos.df, path, sep = ",", quote = FALSE)

# get trios with significantly positive gamma coefficients
# from model 5X level 2 (TRIPOD level 2 matching Xt)
xymats5X.gamma.pos.df <- getTrios(
	xymats.list = xymats.list.list$xymats5X, 
	fdr.thresh = fdr.thresh, 
	sign = sign,
  model.name = "5a",
	level = 2)
file <- paste0("xymats5X.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5X.gamma.pos.df, path, sep = ",", quote = FALSE)

# get trios with significantly positive alpha coefficients
# from model 5Y level 1 (TRIPOD level 1 matching Yj)
xymats5Y.alpha.pos.df <- getTrios(
	xymats.list = xymats.list.list$xymats5Y, 
	fdr.thresh = fdr.thresh, 
	sign = sign,
  model.name = "5a",
	level = 1)
file <- paste0("xymats5Y.alpha.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5Y.alpha.pos.df, path, sep = ",", quote = FALSE)

# get trios with significantly positive gamma coefficients
# from model 5Y level 2 (TRIPOD level 2 matching Yj)
xymats5Y.gamma.pos.df <- getTrios(
	xymats.list = xymats.list.list$xymats5Y, 
	fdr.thresh = fdr.thresh, 
	sign = sign,
  model.name = "5a",
	level = 2)
file <- paste0("xymats5Y.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5Y.gamma.pos.df, path, sep = ",", quote = FALSE)

