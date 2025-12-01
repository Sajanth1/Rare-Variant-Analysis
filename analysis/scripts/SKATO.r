#!/usr/bin/env Rscript
wd <- getwd()

library(packrat)
packrat::init("~/runs/sajanth/SKAT")
library(SKAT)
library(MetaSKAT)

require(methods)
suppressPackageStartupMessages(require(Xmisc)) #Basically, kept getting The following object is masked from ‘package:base’: dir.exists
options(warning.length = 8000)

parser <- ArgumentParser$new()
#### Setup script doc and help messages
parser$add_usage('R.SKAT [options]')
parser$add_description('Setup and run R SKAT for one gene one cohort one depth')
parser$add_argument('--h', type='logical', action='store_true', help='Print the help page')
parser$add_argument('--help', type='logical', action='store_true', help='Print the help page')
parser$add_argument('--cohort', type='character', help='Cohort name eg. AMP_PD', required=TRUE)
parser$add_argument('--gene', type='character', help='Gene name eg. VPS11', required=TRUE)
args <- parser$get_args()

#################################
setwd(wd)
cohort <- args$cohort
gene <- args$gene

File.Bed   = paste0("SKAT/", gene, "/", cohort, ".", gene, ".bed")
File.Bim   = paste0("SKAT/", gene, "/", cohort, ".", gene, ".bim")
File.Fam   = paste0("SKAT/", gene, "/", cohort, ".", gene, ".fam")
File.Cov   = paste0("covar_", cohort,".txt")
File.SetID = paste0("SKAT/", cohort, ".", gene, ".SETID")

File.SSD   = paste0("SKAT/", gene, "/", cohort, ".", gene, ".SSD")
File.Info  = paste0("SKAT/", gene, "/", cohort, ".", gene, ".info")
File.Mat <- paste0("SKAT/", gene, "/", cohort, ".", gene, ".mat")
File.SetInfo <- paste0("SKAT/", gene, "/", cohort, ".", gene, ".MInfo")
File.Results.SKATO  = paste0("SKAT/", gene, "/", cohort, ".", gene, ".results.skato")



Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info) #Merges based on 2nd column of both .bim and SETID

SSD.INFO<-Open_SSD(File.SSD, File.Info)


if(file.exists(File.Cov)){
    message("Analysis with covar")
    FAM <-Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE,  cov_header=TRUE) #Only merges if sample present IN BOTH files, so adjusted fam is sufficient for sub-analyses
    y <-FAM$Phenotype
    N.Sample <-length(y)

    Sex <-FAM$Sex.y
    Age <-FAM$Age
    PC1 <-FAM$pc1
    PC2 <-FAM$pc2
    PC3 <-FAM$pc3
    PC4 <-FAM$pc4
    PC5 <-FAM$pc5
    # PC6 <-FAM$pc6
    # PC7 <-FAM$pc7
    # PC8 <-FAM$pc8
    # PC9 <-FAM$pc9
    # PC10 <-FAM$pc10

    set.seed(2025) #small sample adjustment (n<2000) uses random sampling
    out_type <- "D"
    obj<-SKAT_Null_Model(y ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5, out_type=out_type)

} else {
    message("Analysis without covar")
    FAM<-Read_Plink_FAM(File.Fam, Is.binary = TRUE)
    y<-FAM$Phenotype
    N.Sample<-length(y)
    obj<-SKAT_Null_Model(y ~ 1, out_type="D")
}


# Perform SKAT-O analysis
message("Performing SKAT-O analysis...\n")
out.skato<-SKATBinary_Robust.SSD.All(SSD.INFO, obj, method="SKATO")
#out.skato.burden<-SKATBinary.SSD.All(SSD.INFO, obj, method="Burden")

# Write results
write.table(out.skato$results, file=File.Results.SKATO, col.names = TRUE, row.names = FALSE)

# Save Metadata for METASKAT
re1<-Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.Mat, File.SetInfo, N.Sample)

# Print warnings
w <- warnings()
if (length(w) > 10) {
  print(w)
}
