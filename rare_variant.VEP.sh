#!/bin/bash
#SBATCH --account=def-grouleau
#SBATCH --time=8:00:00           # time (DD-HH:MM:SS)
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G                  # memory per node
#SBATCH --job-name=rare_variant_analysis
#SBATCH --error=logs/job.%x-%j.err
#SBATCH --output=logs/job.%x-%j.out


#---------------------------------------------------------------------------
# Version 26/08/2025 (Rorqual version): Added rename option if not subsetting (works for both cohorts)
# Starting files: rare_variant.VEP.sh, bed

# TO IMPROVE
# Move covar specification to top? (along with sex_UKBB?)
# Covar is inaccurate after 1b (contains non-filtered controls, actually already inaccurate since sampleIDs missing after data extraction)
# SAIGE-GENE and META-SAIGE instead of SKAT package, better for biobank data and case/ctrl imbalances

# Note: In 1b., control numbers should be adjusted (verify via logs/cohort.info)
#
# Check out SKAT_CommonRare (later)
# Cool way of waiting for previous job to finish: while true; do n=$(squeue -h -j JOBID? | grep -c JOBID?); if [[ "$n" != 1 ]]; then sbatch rare_variant.VEP.sh; break; else echo "waiting..."; fi; sleep 300; done
#---------------------------------------------------------------------------

# USER INPUT: Specify cohort and UCSC bed file 
cohort=UKBB
bed=Sajanth44.GRCh38.bed
vcf_raw=~/links/scratch/DYT-PD/UKBB/import/UKBB_raw.vcf.gz
mapfile -t genelist < <(awk '{print $5}' $bed)
# genelist=("AOPEP")


# 0. Set up env
# module purge ; source ~/.bashrc
module load StdEnv/2023 plink/2.00a5.8 perl/5.36.1 gcc/13.3 mariadb/11.5.0 htslib/1.22.1 db/18.1.32 python/3.13.2
pip install pandas numpy scipy --no-index   # Please create virtual env if you don't feel comfortable downloading these packages usersystem-wide

set -e
# cp -r -n ~/runs/sajanth/2025/DYT-PD/DYT-PD/UKBB/analysis/scripts ./
mkdir -p annotation logs SKAT


# # 1. Prepare Dataset

if [[ $cohort == "AMP_PD" ]]; then
    printf "\n\nPreparing AMP_PD dataset...\n\n"
    cp -n ~/runs/sajanth/2025/AMP_PD_data/covar_AMP_PD.txt ./ # Get covar

    # Get bfiles && keep only rare variants (UCSC bed file is 0-based)
    plink2 --bfile ~/runs/sajanth/2025/AMP_PD_data/AMP_PD_FILTERED_ALL_CHR --max-maf 0.01 --rm-dup force-first \
    --double-id --extract bed0 $bed --pheno covar_AMP_PD.txt --pheno-name Status \
    --keep covar_AMP_PD.txt  \
    --make-bed --out AMP_PD_allctrl  
    
    awk 'BEGIN{OFS="\t"} {$2 = $1 ":" $4 ":" $6 ":" $5; print}' AMP_PD_allctrl.bim > tmp && mv tmp AMP_PD_allctrl.bim   # Add unique tag (ref:alt) Keep in mind that .bim file is only bfile that contains variant IDs
    plink2 --bfile AMP_PD_allctrl --recode vcf bgz --out AMP_PD_all # Create VCF
    bcftools view -G --threads 4 -Oz -o AMP_PD.vcf.gz AMP_PD_all.vcf.gz ; tabix AMP_PD.vcf.gz  # Create smaller vcf (by removing genotype columns)


elif [[ $cohort == "UKBB" ]]; then
    printf "\n\nPreparing UKBB dataset...\n\n"
    cp -n ../covar/covar_UKBB.txt ./ # Get covar

    # Get bfiles && keep only rare variants (keep in mind that multiallelic variants were already split during QC) && update sex & case/ctrl in .fam (--keep to subset)
    plink2 --memory 60000 --vcf $vcf_raw --max-maf 0.01 --vcf-half-call m --rm-dup force-first \
    --double-id -update-sex ../covar/sex_UKBB.txt --pheno covar_UKBB.txt --pheno-name Status \
    --keep covar_UKBB.txt \
    --make-bed --out UKBB_allctrl

    plink2 --memory 60000 --bfile UKBB_allctrl --recode vcf bgz --out UKBB_all  # Recreate VCF
    bcftools view -G --threads 4 -Oz -o UKBB.vcf.gz UKBB_all.vcf.gz ; tabix UKBB.vcf.gz  # Create smaller vcf (by removing genotype columns)
fi   


# 1b. Adapt bfiles to each sub-analysis
num_controls=33840  # Number of controls to sample or big number to keep all
name="${cohort}_allctrl"  # bfile name
# awk '$6 == 1' $name.fam | shuf -n $num_controls > logs/sampled_controls.txt # Randomize controls
# awk '$6 == 1' $name.fam | grep -v -f logs/sampled_controls.txt > logs/rm_controls.txt # Get controls to remove
plink2 --memory 60000 --bfile $name --double-id --remove logs/rm_controls.txt --make-bed --out $cohort # Update bfiles

# #Rename if not subsetting controls
# for file in ${cohort}_allctrl.*; do mv "$file" "${cohort}${file#${cohort}_allctrl}"; done

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
