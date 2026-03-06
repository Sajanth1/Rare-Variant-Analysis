#!/bin/bash

#---------------------------------------------------------------------------
# USER INPUT:
cohort=
bed=
mapfile -t genelist < <(awk '{print $5}' $bed)
# genelist=("AOPEP") # manually select gene(s)

vcf_preprocessed=
covar=
cp -n $covar ./ # Get covar

sample_controls= # y or n
num_controls= # Number of controls to sample
#---------------------------------------------------------------------------

# 0. Set up env: mamba env for tools and containers for VEP
#sudo apt update && sudo apt install -y software-properties-common &&  sudo add-apt-repository -y ppa:apptainer/ppa && sudo apt update && sudo apt install -y apptainer  #Install apptainer
#mamba create -n analysis -c bioconda -c conda-forge bcftools=1.16 plink=1.90b6.21 plink2=2.00a5.10 python=3.12.3 pandas numpy scipy r-base=4.4 r-tidyr=1.3.2
#mamba create -n SKAT -c bioconda -c conda-forge r-base=4.3 r-skat=2.2.5 r-metaskat=0.90 r-argparse
unset PYTHONPATH ; eval "$(mamba shell hook --shell bash)" ; mamba activate analysis

set -e
# cp scripts from github
mkdir -p annotation logs SKAT


# 1. Prepare Dataset
if [[ $cohort == "AMP_PD" ]]; then
    # bfile pathway
    printf "\n\nPreparing AMP_PD dataset...\n\n"

    # Get bfiles && keep only rare variants (UCSC bed file is 0-based)
    plink2 --bfile ~/.../AMP_PD_data/ --max-maf 0.01 --mac 1 --rm-dup force-first \
    --double-id --extract bed0 $bed --pheno covar_AMP_PD.txt --pheno-name STATUS \
    --keep covar_AMP_PD.txt  \
    --make-bed --out AMP_PD_allctrl  
    
    awk 'BEGIN{OFS="\t"} {$2 = $1 ":" $4 ":" $6 ":" $5; print}' AMP_PD_allctrl.bim > tmp && mv tmp AMP_PD_allctrl.bim   # Add unique tag (ref:alt) Keep in mind that .bim file is only bfile that contains variant IDs


elif [[ $cohort == "UKBB" || $cohort == "GP2" ]]; then
    # VCF pathway
    printf "\n\nPreparing ${cohort} dataset...\n\n"

    # Get bfiles && keep only rare variants (note: multiallelic variants were already split during QC) && update sex & case/ctrl in .fam (covar must be: FID IID choice1 choice2 AGE SEX PC1-10 STATUS)
    plink2 --memory 60000 --vcf $vcf_preprocessed --max-maf 0.01 --mac 1 --vcf-half-call m --rm-dup force-first \
    --double-id --update-sex covar_${cohort}.txt col-num=6  --pheno covar_${cohort}.txt --pheno-name STATUS \
    --keep covar_${cohort}.txt \
    --make-bed --out ${cohort}_allctrl
fi   


# 1b. Adapt bfiles to each sub-analysis
if [[ $sample_controls == "y" ]]; then
    name="${cohort}_allctrl"  # bfile name
    # awk '$6 == 1' $name.fam | shuf -n $num_controls > logs/sampled_controls.txt # Randomize controls
    # awk '$6 == 1' $name.fam | grep -v -f logs/sampled_controls.txt > logs/rm_controls.txt # Get controls to remove
    plink2 --memory 60000 --bfile $name --double-id --remove logs/rm_controls.txt --mac 1 --make-bed --out $cohort # Update bfiles

elif [[ $sample_controls == "n" ]]; then
    for file in ${cohort}_allctrl.*; do mv "$file" "${cohort}${file#${cohort}_allctrl}"; done #Rename if not subsetting controls
fi


# 1c. (re)Create VCF
plink2 --memory 60000 --bfile $cohort --recode vcf bgz --out ${cohort}_all
bcftools view -G --threads 4 -Oz -o ${cohort}.vcf.gz ${cohort}_all.vcf.gz ; tabix ${cohort}.vcf.gz  # Create smaller vcf (by removing genotype columns)


#Clean up
mv *.log  logs/
mv *_all* logs/



# 2. Annotate
printf "\n\nAnnotating variants with VEP...\n\n"
bash scripts/VEP_annotation.sh $cohort &> logs/VEP_annotation_${cohort}.log


# 3. Set Ids: Group variants for each gene
printf "\n\nGrouping variants...\n\n"
python3 scripts/vep_setid_prep.py $cohort "${genelist[@]}"


# 4. Run SKAT-O on each gene
printf "\n\nRunning SKAT-O...\n\n" 
bash scripts/run_multiskato.sh $cohort "${genelist[@]}" 

#Merge results
rm -f SKAT/${cohort}.results.skato ; find "SKAT/" -type f -name "*.results.skato" -print0 | xargs -0 awk 'FNR==1 && NR!=1 {next} {print}' | { read -r header; echo "$header"; sort -k1,1; } > "SKAT/${cohort}.results.skato.raw"
awk 'FNR==1 || ($5 >= 2 && $5 != "NA") {print}' "SKAT/${cohort}.results.skato.raw" > "SKAT/${cohort}.results.skato" # Added filter for sets with less than 2 variants or SKAT outputs NA
Rscript scripts/FDR.r $cohort # FDR correction


# 5. Single variant Analysis
printf "\n\nRunning single variant analysis...\n\n"
bash scripts/run_singlevariant.sh $cohort

printf "\n\nDone! Your results are in annotation/ \n(Set analysis = *.results.skato.adjusted, single variant ORs = *.single_variants.csv) \nand intermediate files in SKAT/"
