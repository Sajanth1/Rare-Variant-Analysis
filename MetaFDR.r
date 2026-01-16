#!/usr/bin/env Rscript
#library(packrat)
#packrat::init("/home/sajanth/scratch/AMP_PD/SKAT")

library(stats)
library(tidyr)


# Load your data files
all_genes <- read.table(paste0("METASKAT/results.metaskato"), sep = " ", header = TRUE)

#Reformat table to have _type as 2nd column; can modify to pivot table where all p values are columns in same row (=gene)
extract_components <- function(data) {
  data <- separate(data, SetID, into = c("Gene", "Variant_Type"), sep = "_", remove = FALSE)
  return(data)
}

# Apply the extract_components function to the data
all_genes <- extract_components(all_genes)

# Adjust the p-values using the FDR method
all_genes$adjusted_P.value <- p.adjust(all_genes$p.value, method = "fdr")

write.table(all_genes, paste0("METASKAT/results.metaskato.adjusted.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
