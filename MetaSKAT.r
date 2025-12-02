#!/usr/bin/env Rscript
wd <- getwd()

library(packrat)
packrat::init("~/runs/sajanth/SKAT")
library(SKAT)
library(MetaSKAT)
require(methods)
suppressPackageStartupMessages(require(Xmisc)) #Basically, kept getting The following object is masked from 'package:base': dir.exists

#-----------------------------------------------
parser <- ArgumentParser$new()
parser$add_argument('--gene', type='character', help='Gene name eg. VPS11', required=TRUE)
parser$add_argument('--cohort', type='character', nargs='+', help='Cohort names', required=TRUE)
args <- parser$get_args()
#-----------------------------------------------
setwd(wd)
gene <- args$gene
cohort_list <- args$cohort
File.Meta.SKATO = paste0("METASKAT/", gene, "/", gene, ".results.metaskato")
set.seed(2025)

File.Mat.vec <- character(length(cohort_list))
File.SetInfo.vec <- character(length(cohort_list))

for (i in seq_along(cohort_list)) {
    cohort <- cohort_list[i]
    File.Mat.vec[i] <- paste0(cohort, "/analysis/SKAT/", gene, "/", cohort, ".", gene, ".mat")
    File.SetInfo.vec[i] <- paste0(cohort, "/analysis/SKAT/", gene, "/", cohort, ".", gene, ".MInfo")
}


Cohort.Info <- Open_MSSD_File_2Read(File.Mat.vec, File.SetInfo.vec)
out.metaskato <- MetaSKAT_MSSD_ALL(Cohort.Info, method="optimal") # apparently optimal is equivalent to optimal.adj in this case
write.table(out.metaskato, file= File.Meta.SKATO, col.names = TRUE, row.names = FALSE)
