#!/bin/bash
file=$1
keep_samples=$2
bed_file=$3

module load StdEnv/2020 gcc/9.3.0 bcftools/1.16

name=$(basename $file .vcf.gz)

# workflow
## 1. subset vcf file by samples
## 2. label variants with DP<25/GQ<25 with missing genotype
## 3. split multiallelic variants - otherwise plink can't handle the multiallelic variants. you can remove them if you want
## 4. convert to plink file (also removes INFO field)
## 5. do missingness filtering
## 6. rename ID column of the .bim file to chr:bp:A2:A1
## 7. convert back to vcf file

# you can delete this step if your file is subset
echo "subsetting samples"
bcftools view --threads 10 --samples-file ${keep_samples} -Oz -o ${name}.subset.vcf.gz -R ${bed_file} $file --force-samples 
tabix ${name}.subset.vcf.gz

# DP < 25/ GQ < 25 genotypes will be marked as missing in order
echo "label bad variants missing ./."
bcftools +setGT ${name}.subset.vcf.gz -Ou --threads 10 \
  -- -t q -n "." -i 'FMT/DP<25' \
| \
bcftools +setGT - -Ou --threads 10 \
  -- -t q -n "." -i 'FMT/GQ<25' \
|  \
bcftools view -Oz -o $name.subset.DP25.GQ25.vcf.gz
## if you wanna see the changes
# bcftools query -f '%CHROM\t%POS[\t%GT:%DP:%GQ]\n' $name.subset.DP25.GQ25.vcf.gz| head|cut -f1-10 -d":"

# split multiallelic variants
echo "splitting multiallelic variants"
bcftools norm --threads 5 -m -both -Oz -o $name.subset.DP25.GQ25.split.vcf.gz $name.subset.DP25.GQ25.vcf.gz

# convert to plink file and do 5% missingness filtering
## plink2 can handle multiallelic variants? 
echo "converting vcf to plink files"
if [[ "$name" == *"cX"* ]]; then
  # plink2 can't well handle chrX data without --split-par
  echo "chr X detected"
  plink2 --vcf $name.subset.DP25.GQ25.split.vcf.gz --vcf-half-call m --split-par b38 --update-sex sex_for_plink.txt --make-bed --out $name.subset.DP25.GQ25.split
  awk '$5 == 1' $name.subset.DP25.GQ25.split.fam > males.txt
  awk '$5 == 2' $name.subset.DP25.GQ25.split.fam > females.txt
  # For females (diploid X)
  plink2 --bfile $name.subset.DP25.GQ25.split --keep females.txt --missing --out chrX_females

  # For males (hemizygous X)
  plink2 --bfile $name.subset.DP25.GQ25.split --keep males.txt --missing --out chrX_males
  awk '$5 < 0.05 { print $2 }' chrX_females.vmiss > chrX_pass_f.txt
  awk '$5 < 0.05 { print $2 }' chrX_males.vmiss > chrX_pass_m.txt
  # Keep only SNPs present in both lists
  sort chrX_pass_f.txt chrX_pass_m.txt | uniq -d > chrX_pass_snps.txt

  plink2 --bfile $name.subset.DP25.GQ25.split \
  --extract chrX_pass_snps.txt \
  --make-bed \
  --out $name.subset.DP25.GQ25.MISS95.split
else
  plink2 --vcf $name.subset.DP25.GQ25.split.vcf.gz --vcf-half-call m --make-bed --out $name.subset.DP25.GQ25.split
  plink2 --bfile $name.subset.DP25.GQ25.split --geno 0.05 --make-bed --out $name.subset.DP25.GQ25.MISS95.split

fi

# reformat the ID column from chr:bp:SG/IG to chr:bp:ref:alt (A2:A1)
awk 'BEGIN{OFS="\t"} {$2= $1 ":" $4 ":" $6 ":" $5; print}' $name.subset.DP25.GQ25.MISS95.split.bim > $name.subset.DP25.GQ25.MISS95.split.renamed.bim
mv $name.subset.DP25.GQ25.MISS95.split.renamed.bim $name.subset.DP25.GQ25.MISS95.split.bim

# convert back to vcf file
plink2 --bfile $name.subset.DP25.GQ25.MISS95.split --export vcf bgz id-paste=iid --out $name.subset.DP25.GQ25.MISS95.split

tabix $name.subset.DP25.GQ25.MISS95.split.vcf.gz


rm $name.subset.DP25.GQ25.MISS95.split.bim
rm $name.subset.DP25.GQ25.MISS95.split.bed
rm $name.subset.DP25.GQ25.MISS95.split.fam
rm $name.subset.DP25.GQ25.MISS95.split.log
rm $name.subset.DP25.GQ25.split.bim
rm $name.subset.DP25.GQ25.split.bed
rm $name.subset.DP25.GQ25.split.fam
rm $name.subset.DP25.GQ25.split.log
rm $name.subset.DP25.GQ25.split.vcf.gz
rm $name.subset.DP25.GQ25.vcf.gz
rm ${name}.subset.vcf.gz
rm ${name}.subset.vcf.gz.tbi