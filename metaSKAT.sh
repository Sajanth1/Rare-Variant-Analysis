#!/bin/bash
#SBATCH --account=def-grouleau
#SBATCH --time=2:0:00           # time (DD-HH:MM:SS)
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                  # memory per node
#SBATCH --job-name=metaskat_analysis
#SBATCH --error=METASKAT/logs/job.%x-%j.err
#SBATCH --output=METASKAT/logs/job.%x-%j.out

# Run from scratch with all cohorts as children dir
set -e
module load StdEnv/2020 nixpkgs/16.09 gcc/7.3 r/3.5.2

#--------- USER INPUT ---------
cohorts=("AMP_PD" "UKBB")   # Specify your cohorts in list format here
bed=~/path/to/genelist.bed
mapfile -t genelist < <(awk '{print $5}' $bed)
#------------------------------

for gene in "${genelist[@]}"; do
    mkdir -p "METASKAT/$gene" METASKAT/logs
    Rscript metaSKAT.r --cohort "${cohorts[@]}" --gene "$gene" 2>&1 | tee "METASKAT/$gene/MetaSKAT.$gene.log"
    
    # Parallelize per gene:
    # command="module load StdEnv/2020 nixpkgs/16.09 gcc/7.3 r/3.5.2 && Rscript MetaSKAT.r --gene $gene 2>&1 | tee 'METASKAT/$gene/MetaSKAT.$gene.log' && echo 'Done!'"
    # sbatch -c 4 --mem=64g -t 1:0:0 --wrap "$command" --account=def-grouleau --out metaskat_$gene.out
done

#Merge results & FDR
module load StdEnv/2023 r/4.4.0 r-bundle-bioconductor/3.20
find "METASKAT/" -type f -name "*.results.metaskato" -print0 | xargs -0 awk 'FNR==1 && NR!=1 {next} {print}' | { read -r header; echo "$header"; sort -k1,1; } | grep -wv "NA" | tr -d '"' > "METASKAT/results.metaskato.raw"


# ---------------------------------------
# Merge SKAT variant numbers across cohorts; this all assumes that decision about whether variants are tested is done at cohort-lvl (ie. monomorphic & missingness filters done only at cohort lvl)
echo "SetID AMP_PD_all AMP_PD_tested UKBB_all UKBB_tested total_all total_tested" > METASKAT/logs/set_lengths.txt
join -a 1 -a 2 -e NA -o 0,1.2,1.3,2.2,2.3 \
  <(tail -n +2 AMP_PD/analysis/SKAT/AMP_PD.results.skato.raw | awk '{print $1, $4, $5}' | sort) \
  <(tail -n +2 UKBB/analysis/SKAT/UKBB.results.skato.raw | awk '{print $1, $4, $5}' | sort) \
  >> METASKAT/logs/set_lengths.txt

awk 'NR==1{print; next} {$6 = ($2=="NA"?0:$2) + ($4=="NA"?0:$4); $7 = ($3=="NA"?0:$3) + ($5=="NA"?0:$5); print}' METASKAT/logs/set_lengths.txt | tr -d '"' > tmp && mv tmp METASKAT/logs/set_lengths.txt # Sum up while counting NA as 0

# In reality, anything can reduce to 1 as filtration criteria are stricter than SKAT. MANUAL CHECK needed
awk '$7 <= 1 {print $1}' METASKAT/logs/set_lengths.txt > METASKAT/logs/singleton_and_empty_sets.txt
grep -wv -f METASKAT/logs/singleton_and_empty_sets.txt "METASKAT/results.metaskato.raw" | grep -wv "1" > "METASKAT/results.metaskato" ## skip: cp METASKAT/results.metaskato.raw METASKAT/results.metaskato
# ---------------------------------------

Rscript metaFDR.r

#Verify: variants from each cohort, totals ; raw p-values after metaSKAT, adjusted p-values of sets that were kept.
join -a 1 METASKAT/logs/set_lengths.txt <(join -a 1 METASKAT/results.metaskato.raw METASKAT/results.metaskato.adjusted.txt -o 0,0,1.2,2.5) > METASKAT/manualcheck.txt

echo "Done! Your results are in METASKAT/results.metaskato.adjusted"
