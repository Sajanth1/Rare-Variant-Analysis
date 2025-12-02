#!/bin/bash
#SBATCH --account=?
#SBATCH --time=8:00:00           # time (DD-HH:MM:SS)
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G                  # memory per node
#SBATCH --job-name=rare_variant_analysis
#SBATCH --error=logs/job.%x-%j.err
#SBATCH --output=logs/job.%x-%j.out
#-------------------------------------------------------------------------------------
# Version 26/08/2025: Added rename option if not subsetting (works for both cohorts)
#-------------------------------------------------------------------------------------

# USER INPUT: Specify cohort, UCSC bed file, covar and VCF file (named vcf_preprocess despite the fact that it must be processed via sample-level and variant-level QC first)
cohort=UKBB
bed=~/path/to/genelist.bed

# covar MUST be named in the following format: "covar_COHORT.txt" where COHORT is the name specified above. It is assumed that your covar has at least the following columns: "FID IID Sex Age PC1 PC2 PC3 PC4 PC5 Status" where Status is coded as 1=control, 2=case. If not, update SKATO.r script accordingly.
covar=~/path/to/covar_UKBB.txt 
num_controls=100000  # Number of controls to sample or big number to keep all

vcf_preprocess=~/path/to/UKBB_preprocess.vcf.gz # Option 1
bfile_preprocess=~/path/to/AMP_PD_preprocess # Option 2
mapfile -t genelist < <(awk '{print $5}' $bed)



# 0. Set up env
# module purge ; source ~/.bashrc
module load StdEnv/2023 plink/2.00a5.8 perl/5.36.1 gcc/13.3 mariadb/11.5.0 htslib/1.22.1 db/18.1.32 python/3.13.2
pip install pandas numpy scipy --no-index   # Please create virtual env if you don't feel comfortable downloading these packages usersystem-wide

set -e
mkdir -p annotation logs SKAT


# # 1. Prepare Dataset

if [[ $cohort == "AMP_PD" ]]; then
    printf "\n\nPreparing AMP_PD dataset...\n\n"
    cp -n $covar ./ 

    # Get bfiles && keep only rare variants (UCSC bed file is 0-based)
    plink2 --bfile $bile_preprocess --max-maf 0.01 --rm-dup force-first \
    --double-id --extract bed0 $bed --pheno $covar --pheno-name Status \
    --keep $covar  \
    --make-bed --out AMP_PD_allctrl  
    
    awk 'BEGIN{OFS="\t"} {$2 = $1 ":" $4 ":" $6 ":" $5; print}' AMP_PD_allctrl.bim > tmp && mv tmp AMP_PD_allctrl.bim   # Add unique tag (chr:start:ref:alt) Keep in mind that .bim file is only bfile that contains variant IDs
    plink2 --bfile AMP_PD_allctrl --recode vcf bgz --out AMP_PD_all # Create VCF
    bcftools view -G --threads 4 -Oz -o AMP_PD.vcf.gz AMP_PD_all.vcf.gz ; tabix AMP_PD.vcf.gz  # Create smaller vcf (by removing genotype columns)


elif [[ $cohort == "UKBB" ]]; then
    printf "\n\nPreparing UKBB dataset...\n\n"
    cp -n $covar ./

    # Create Sex file for plink
    awk 'BEGIN{FS=OFS="\t"} FNR==1{print "FID", "IID", "Sex"} FNR!=1{print $1, $2, $7}' $covar > sex_UKBB.txt # construct your sex file such that FID, IID, Sex are the columns

    # Get bfiles && keep only rare variants (keep in mind that multiallelic variants were already split during QC) && update sex & case/ctrl in .fam (--keep to subset)
    plink2 --memory 60000 --vcf $vcf_preprocess --max-maf 0.01 --vcf-half-call m --rm-dup force-first \
    --double-id -update-sex sex_UKBB.txt --pheno $covar --pheno-name Status \
    --keep $covar \
    --make-bed --out UKBB_allctrl

    plink2 --memory 60000 --bfile UKBB_allctrl --recode vcf bgz --out UKBB_all  # Recreate VCF
    bcftools view -G --threads 4 -Oz -o UKBB.vcf.gz UKBB_all.vcf.gz ; tabix UKBB.vcf.gz  # Create smaller vcf (by removing genotype columns)
fi   


# 1b. Adapt bfiles to each sub-analysis
name="${cohort}_allctrl"  # bfile name
awk '$6 == 1' $name.fam | shuf -n $num_controls > logs/sampled_controls.txt # Randomize controls
awk '$6 == 1' $name.fam | grep -v -f logs/sampled_controls.txt > logs/rm_controls.txt # Get controls to remove
plink2 --memory 60000 --bfile $name --double-id --remove logs/rm_controls.txt --make-bed --out $cohort # Update bfiles


#Clean up
mv *.log  logs/
mv *_all* logs/



# 2. Annotate
printf "\n\nAnnotating variants with VEP...\n\n"
bash scripts/VEP114_annotation.sh $cohort &> logs/VEP_annotation_${cohort}.log


# 3. Set Ids: Group variants for each gene
printf "\n\nGrouping variants...\n\n"
python3 scripts/vep_setid_prep.py $cohort "${genelist[@]}"


# 4. Run SKAT-O on each gene
printf "\n\nRunning SKAT-O...\n\n" 
bash scripts/run_multiskato.sh $cohort "${genelist[@]}" 

#Merge results
module load StdEnv/2023 plink/2.00a5.8 python/3.13.2 r-bundle-bioconductor/3.20 r/4.4.0
rm -f SKAT/${cohort}.results.skato ; find "SKAT/" -type f -name "*.results.skato" -print0 | xargs -0 awk 'FNR==1 && NR!=1 {next} {print}' | { read -r header; echo "$header"; sort -k1,1; } > "SKAT/${cohort}.results.skato.raw"
awk 'FNR==1 || ($5 >= 2 && $5 != "NA" && $2 != 1) {print}' "SKAT/${cohort}.results.skato.raw" > "SKAT/${cohort}.results.skato" # Added filter for sets with less than 2 variants or SKAT outputs exactly 1 or NA
Rscript scripts/FDR.r $cohort # FDR correction


# 5. Single variant Analysis
printf "\n\nRunning single variant analysis...\n\n"
bash scripts/run_singlevariant.sh $cohort

printf "\n\nDone! Your results are in annotation/ \n(Set analysis = *.results.skato.adjusted, single variant ORs = *.single_variants.csv) \nand intermediate files in SKAT/"
