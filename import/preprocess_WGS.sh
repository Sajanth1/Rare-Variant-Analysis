#!/bin/bash
#---------
# Version (30/04/2026): ancestry mgmt
#
#--------

cohort=GP2
bed_file=/config/KYN16.GRCh38.bed

#GP2's R11 existing QC for vcf: • Left-aligned and normalized indels • Multiallelic split • Filtered by 'PASS' flag

# workflow:
## 1. subset vcf file by samples & bed file range
## 2. label variants with DP<15 or GQ<25 or imbalanced heterozygous with missing genotype
## 3. convert to plink file and do missingness and mac 1 filtering
## 4. rename ID column of the .bim file to 1:bp:A2:A1. (it was chr1:bp:SG/IG)
## 5. convert back to vcf file

#wb resource mount
# -> env: mamba create -n import -c bioconda -c conda-forge bcftools=1.16 plink2=2.00a5.10
eval "$(mamba shell hook --shell bash)" ; mamba activate import


for ancestry in "AFR" "AJ" "AMR" "CAS" "EAS" "EUR" "MDE" "SAS" "CAH"; do

  printf "\n\n\n\nProcessing ${ancestry}...\n\n"
  name="${ancestry}/${ancestry}_${cohort}"
  keep_samples="../covar/${ancestry}_sampleIDs_GP2.txt"
  mkdir -p ${ancestry}/vcfs


  # 1a. Subset
  if [[ -f "${ancestry}/vcfs_rerun.txt" ]]; then
    input_list="${ancestry}/vcfs_rerun.txt"
  else
    ls /config/workspace/gp2_tier2_eu_release11/wgs/deepvariant_joint_calling/vcfs/${ancestry}/*.vcf.gz | grep -vE "chrX|chrY" | sort -V > "${ancestry}/vcfs_list.txt"
    input_list="${ancestry}/vcfs_list.txt"
  fi

  while read -r file; do
    echo "subsetting $(basename "$file")"
    out=${ancestry}/vcfs/$(basename "$file")
    
    bcftools view --threads 4 --samples-file ${keep_samples} -Oz -o ${out} -R ${bed_file} $file --force-samples ; tabix ${out}
  done < $input_list


  # 1b. Merge
  echo "merging vcfs"
  ls ${ancestry}/vcfs/*.vcf.gz | sort -V > ${ancestry}/merge_list.txt
  bcftools concat --naive -f ${ancestry}/merge_list.txt --threads 10 -Oz -o ${name}_merged.vcf.gz && echo 'concatenation done' ; tabix ${name}_merged.vcf.gz



  # 2. DP < 15 or GQ < 25 or imbalanced heterozygous genotypes will be marked as missing in order; used VAF column (pre-created by GP2) is defined as "The fraction of reads with alternate allele (nALT/nSumAll)">
  echo "label bad variants missing ./."
  bcftools +setGT ${name}_merged.vcf.gz -Ou --threads 10 \
    -- -t q -n "." -i 'FMT/DP<15' \
  | \
  bcftools +setGT - -Ou --threads 10 \
    -- -t q -n "." -i 'FMT/GQ<25' \
  | \
  bcftools +setGT - -Ou --threads 10 \
    -- -t q -n "." -i 'GT="het" & FMT/VAF < 0.15' \
  | \
  bcftools +setGT - -Ou --threads 10 \
    -- -t q -n "." -i 'GT="het" & FMT/VAF > 0.85' \
  | \
  bcftools view -Oz -o ${name}_DP25_GQ25_AB.vcf.gz
  ## if you wanna see the changes
  # bcftools query -f '%CHROM\t%POS[\t%GT:%DP:%GQ:%AD:%VAF]\n' $name.subset.DP25.GQ25.vcf.gz| head|cut -f1-17 -d":"


  # 3. convert to bfiles and do 5% missingness and mac 1 filtering
  echo "converting vcf to plink files"
  plink2 --vcf ${name}_DP25_GQ25_AB.vcf.gz --vcf-half-call m --geno 0.05 --mac 1 --make-bed --out ${name}_DP25_GQ25_AB_MISS95


  # 4. reformat the ID column from chr1:bp:SG/IG to 1:bp:ref:alt (A2:A1) as well as the chr column
  awk 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1); $2= $1 ":" $4 ":" $6 ":" $5; print}' ${name}_DP25_GQ25_AB_MISS95.bim > ${name}_DP25_GQ25_AB_MISS95.renamed.bim ; mv ${name}_DP25_GQ25_AB_MISS95.renamed.bim ${name}_DP25_GQ25_AB_MISS95.bim

  # 5. convert back to vcf file
  plink2 --bfile ${name}_DP25_GQ25_AB_MISS95 --export vcf bgz id-paste=iid --out ${name} ; tabix ${name}.vcf.gz


  rm -r ${ancestry}/vcfs
  rm ${name}_merged.vcf.gz*
  rm ${name}_DP25_GQ25_AB.vcf.gz
  rm ${name}_DP25_GQ25_AB_MISS95*

done
