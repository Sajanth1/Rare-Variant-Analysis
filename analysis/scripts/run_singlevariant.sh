#!/bin/bash

# Version 11/8/2025

cohort="$1"

cohort_cases=$(awk '$6 ==2 {print}' ${cohort}.fam | wc -l)
cohort_ctrls=$(awk '$6 ==1 {print}' ${cohort}.fam | wc -l)

printf "Cohort: $cohort\nCases: $cohort_cases\nControls: $cohort_ctrls" > logs/cohort.info

#Run associations
plink --memory 60000 --bfile ${cohort} --assoc fisher --counts --freq --out annotation/${cohort} --output-chr M
plink --memory 60000 --bfile ${cohort} --assoc fisher --out annotation/${cohort}_MAF --output-chr M

mv annotation/*.log logs/
mv "logs/${cohort}.log" "logs/${cohort}_counts.log"

#Show Variants
python3 scripts/vep_singlevariants.py $cohort $cohort_cases $cohort_ctrls
