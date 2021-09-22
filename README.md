# TRIPOD
Detecting Transcriptional Regulatory Relationships in Single-Cell RNA and Chromatin Accessibility Multiomic Data

## Author
Yuriko Harigaya, Nancy R. Zhang, Yuchao Jiang

## Maintainer
Yuriko Harigaya <harigaya@email.unc.edu>

## Description
TRIPOD is a statistical framework for detecting three-way regulatory relationships between a cis-regulatory region, a transcription factor, and a target gene using single-cell multiomic data. The main functionality of this package is as follows.

* Infers trio regulatory relationships using robust nonparametric models
* Builds RNA prediction models from ATAC-seq data
* Identifies cell-types in which a trio regulatory relationship is active by estimating the influence of a data point in linear regression

## Installation
```r
install.packages("devtools")
devtools::install_github("yharigaya/TRIPOD/package")
```

## Vignettes
* [Preprocessing vignette](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/preprocessing.html)
* [TRIPOD vignette](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/TRIPOD.html)


##  Manuscript code

* [10X Genomics PBMC](https://github.com/yharigaya/TRIPOD/tree/main/scripts/10x_pbmc/)
* [10X Genomics mouse embryonic brain](https://github.com/yharigaya/TRIPOD/tree/main/scripts/10x_e18/)
* [SHARE-seq mouse skin](https://github.com/yharigaya/TRIPOD/tree/main/scripts/share_seq_skin/)
* [SNARE-seq adult mouse brain](https://github.com/yharigaya/TRIPOD/tree/main/scripts/snare_seq_mbrain/)

## Reference
Yuriko Harigaya, Zhaojun Zhang, Chi-Yun Wu, Kevin Z. Lin, Hongpan Zhang, Chongzhi Zang, Nancy R. Zhang, Yuchao Jiang. Nonparametric Interrogation of Transcriptional Regulation in Single-Cell RNA and Chromatin Accessibility Data. ***bioRxiv***, 2021.
