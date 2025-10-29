#!/bin/bash

# Define the cohort
cohort="$1"

# Get the list of genes
genelist=("${@:2}")

# Create subsets of bfiles for each gene (feeding entire bed file seems to create issues with SKAT)
for gene in "${genelist[@]}"; do
  #rm -r "SKAT/$gene" "logs/SKAT_${gene}.log"
  mkdir "SKAT/$gene"
  plink2 --bfile $cohort --double-id --extract <(cut -f2 "SKAT/$cohort.$gene.SETID") --make-bed --out "SKAT/$gene/$cohort.$gene" 
done


# Loop through each gene and run the R script
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3 r/3.5.2

for gene in "${genelist[@]}"; do
  Rscript scripts/SKATO.r --cohort "$cohort" --gene "$gene" 2>&1 | tee "logs/SKAT_${gene}.log"
done
