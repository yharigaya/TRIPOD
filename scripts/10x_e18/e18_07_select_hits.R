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
file <- "xymats.m.list.rds"
path <- file.path(dir.out, file)
xymats.m.list <- readRDS(path); rm(path)

file <- "xymats.i.list.rds"
path <- file.path(dir.out, file)
xymats.i.list <- readRDS(path); rm(path)

file <- "xymats.tripod.Xt.list.rds"
path <- file.path(dir.out, file)
xymats.tripod.Xt.list <- readRDS(path); rm(path)

file <- "xymats.tripod.Yj.list.rds"
path <- file.path(dir.out, file)
xymats.tripod.Yj.list <- readRDS(path); rm(path)

# set FDR < 0.01
fdr.thresh <- 0.01

file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)
xymats <- xymats$Sox10

# set the sign of coefficients
sign <- "positive"
sign.str <- "pos"

# get a data frame of hits from the marginal model
xymats.alpha.m.pos.df <- getPeakGenePairs(
	xymats.list = xymats.m.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "marginal")
file <- paste0("xymats_alpha_m_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.alpha.m.pos.df, path); rm(path)

xymats.beta.m.pos.df <- getTFGenePairs(
	xymats.list = xymats.m.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "marginal")
file <- paste0("xymats_beta_m_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.beta.m.pos.df, path); rm(path)

# get a data fame of hits from the interaction model
xymats.i.pos.df <- getTrios(
	xymats.list = xymats.i.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "interaction"
)
file <- paste0("xymats_i_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.i.pos.df, path); rm(path)

# get a data fame of hits from TRIPOD level 1 matching Xt
xymats.tX1.pos.df <- getTrios(
	xymats.list = xymats.tripod.Xt.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 1
)
file <- paste0("xymats_tX1_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tX1.pos.df, path); rm(path)

# get a data fame of hits from TRIPOD level 2 matching Xt
xymats.tX2.pos.df <- getTrios(
	xymats.list = xymats.tripod.Xt.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 2
)
file <- paste0("xymats_tX2_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tX2.pos.df, path); rm(path)

# get a data fame of hits from TRIPOD level 1 matching Yj
xymats.tY1.pos.df <- getTrios(
	xymats.list = xymats.tripod.Yj.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 1
)
file <- paste0("xymats_tY1_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tY1.pos.df, path); rm(path)

# get a data fame of hits from TRIPOD level 2 matching Yj
xymats.tY2.pos.df <- getTrios(
	xymats.list = xymats.tripod.Yj.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 2
)
file <- paste0("xymats_tY2_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tY2.pos.df, path); rm(path)

# set the sign of coefficients
sign <- "negative"
sign.str <- "neg"

# get a data frame of hits from the marginal model
xymats.alpha.m.neg.df <- getPeakGenePairs(
	xymats.list = xymats.m.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "marginal")
file <- paste0("xymats_alpha_m_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.alpha.m.neg.df, path); rm(path)

xymats.beta.m.neg.df <- getTFGenePairs(
	xymats.list = xymats.m.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "marginal")
file <- paste0("xymats_beta_m_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.beta.m.neg.df, path); rm(path)

# get a data fame of hits from the interaction model
xymats.i.neg.df <- getTrios(
	xymats.list = xymats.i.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "interaction"
)
file <- paste0("xymats_i_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.i.neg.df, path); rm(path)

# get a data fame of hits from TRIPOD level 1 matching Xt
xymats.tX1.neg.df <- getTrios(
	xymats.list = xymats.tripod.Xt.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 1
)
file <- paste0("xymats_tX1_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tX1.neg.df, path); rm(path)

# get a data fame of hits from TRIPOD level 2 matching Xt
xymats.tX2.neg.df <- getTrios(
	xymats.list = xymats.tripod.Xt.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 2
)
file <- paste0("xymats_tX2_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tX2.neg.df, path); rm(path)

# get a data fame of hits from TRIPOD level 1 matching Yj
xymats.tY1.neg.df <- getTrios(
	xymats.list = xymats.tripod.Yj.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 1
)
file <- paste0("xymats_tY1_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tY1.neg.df, path); rm(path)

# get a data fame of hits from TRIPOD level 2 matching Yj
xymats.tY2.neg.df <- getTrios(
	xymats.list = xymats.tripod.Yj.list,
	fdr.thresh = fdr.thresh,
	sign = sign,
  model.name = "TRIPOD",
	level = 2
)
file <- paste0("xymats_tY2_", sign.str, "_", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(xymats.tY2.neg.df, path); rm(path)
