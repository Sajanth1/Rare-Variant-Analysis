#!/usr/bin/env Rscript
#library(packrat)
#packrat::init("/home/sajanth/scratch/AMP_PD/SKAT")

library(stats)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)
cohort <- args[1]

# Load your data files
all_genes <- read.table(paste0("SKAT/", cohort, ".results.skato"), sep = " ", header = TRUE)

#Reformat table to have _type as 2nd column; can modify to pivot table where all p values are columns in same row (=gene)
extract_components <- function(data) {
  data <- separate(data, SetID, into = c("Gene", "Variant_Type"), sep = "_", remove = FALSE)
  return(data)
}

# Apply the extract_components function to the data
all_genes <- extract_components(all_genes)
#print(all_genes)

# Adjust the p-values using the FDR method
all_genes$adjusted_P.value <- p.adjust(all_genes$P.value, method = "fdr")
# Reorder columns to place adjusted_P.value in the fourth position
all_genes <- all_genes[, c(1:4, ncol(all_genes), 5:(ncol(all_genes) - 1))]
#print(all_genes)
write.table(all_genes, paste0("annotation/", cohort, ".results.skato.adjusted.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
