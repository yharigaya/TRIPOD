# TRIPOD
Detecting transcriptional regulatory relationships in single-cell RNA and chromatin accessibility multiomic data

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
devtools::install_github("yharigaya/TRIPOD")
```

## Vignettes & manuscript code
* [Preprocessing vignette](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/preprocessing.html)
* [TRIPOD vignette](http://htmlpreview.github.io/?https://github.com/yharigaya/TRIPOD/blob/main/vignettes/TRIPOD.html)
* [Manuscript code](https://github.com/yharigaya/TRIPOD_manuscript)

## Manuscript
Yuriko Harigaya, Zhaojun Zhang, Chi-Yun Wu, Kevin Z. Lin, Hongpan Zhang, Chongzhi Zang, Nancy R. Zhang, Yuchao Jiang. Nonparametric Interrogation of Transcriptional Regulation in Single-Cell RNA and Chromatin Accessibility Data. ***bioRxiv***, 2021.
