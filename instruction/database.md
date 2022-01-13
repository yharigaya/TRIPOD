## Annotated motifs in addition to the JASPAR database?
  
Motif databases other than the JASPAR database can be used with
`TRIPOD` after the data is converted to
[the JASPAR format](http://jaspar.genereg.net/docs/). 
This document shows an example script 
for running `TRIPOD` on the mouse embryonic brain data with
[the HOCOMOCO database](https://hocomoco11.autosome.ru/downloads_v11).

The mouse motif data can be downloaded as follows.

```
curl -O https://hocomoco11.autosome.ru/final_bundle/hocomoco11/\
 full/MOUSE/mono/HOCOMOCOv11_full_MOUSE_mono_jaspar_format.txt
curl -O https://hocomoco11.autosome.ru/final_bundle/hocomoco11/\
 full/MOUSE/mono/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv
```

The input `Seurat` object, which has been processed from multiome
data, can be downloaded [here](https://www.dropbox.com/s/4afi9rp4t5d5km0/e18.chromvar.rds?dl=0).
Note that `<path/to/seurat>`, `<path/to/motif>`, and `<path/to/anno>`
need to be replaced with paths to the Seurat file, motif file, and
annotation file, respectively.

```
library(tidyverse)
library(Seurat)
library(Signac)
library(TRIPOD)
library(TFBSTools)

# read in the Seurat object
seurat <- <path/to/seurat> %>% readRDS()

pfm.hocomoco <- <path/to/motif> %>% readLines()

formatPFM <- function(x, nucleotide) {
    if (!grepl("^>", x)) {
        x <- gsub("\t", " ", x)
	x <- paste0(nucleotide, " [", x, " ]")
    }
    x
}

input <- pfm.hocomoco
output <- rep(NA, length(input))
for (i in 1:length(input)) {
    if (i %% 5 == 1) {
        output[i] <- input[i]
    } else if (i %% 5 == 2) {
  	output[i] <- formatPFM(input[i], "A")
    } else if (i %% 5 == 3) {
        output[i] <- formatPFM(input[i], "C")
    } else if (i %% 5 == 4) {
        output[i] <- formatPFM(input[i], "G")
    } else if (i %% 5 == 0) {
        output[i] <- formatPFM(input[i], "T")
    }
}

# save the data to a file
dir.out <- "hocomoco"
if (!dir.exists(dir.out)) dir.create(dir.out)
file.out <- "pfm.hocomoco.full.mouse.txt"
file.path(dir.out, file.out) %>%
    write.table(output, file = ., quote = FALSE, row.names = FALSE, col.names = FALSE)

# read back in the file
pfm.set.hocomoco <- file.path(dir.out, file.out) %>%
    readJASPARMatrix(matrixClass = "PFM")

# get TF names
anno <- <path/to/anno> %>%
    read.delim()

for (i in 1:length(pfm.set.hocomoco)) {
    pfm.set.hocomoco[[i]]@name <- anno$Transcription.factor[i]
}

# read in the Seurat object
seurat <- <path/to/seurat> %>% readRDS()

# create a motif object
motif.matrix <- CreateMotifMatrix(
    features = granges(seurat),
    pwm = pfm.set.hocomoco,
    genome = "mm10",
    use.counts = FALSE
)
motif.object <- CreateMotifObject(
    data = motif.matrix,
    pwm = pfm.set.hocomoco
)
```


