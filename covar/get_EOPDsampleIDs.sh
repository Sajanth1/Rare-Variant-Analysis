#!/bin/bash

#Run from covar
module load StdEnv/2020 gcc/9.3.0 r-bundle-bioconductor/3.14 r/4.1.2 rust/1.47.0
# cargo install xsv # (v0.13.0)
# export PATH="$HOME/.cargo/bin:$PATH"


# This script assumes that covar headers are ordered this way: "eid" "20002" "41270" "20111" "20110" "20107" "22001" "22006" "22009" "22000" "34" "21022" "22021" "22019" "22027" "22189" "40001" "40002" "26260" "26261" "42032" "42033" "53"
ukb="/home/sajanth//scratch/DYT-PD.run.early/covar/PD_covar.csv"


## PD patients

#xsv headers $ukb | grep '\s20002-' (Self-report) gives col numbers 2-137
self_report=$(xsv select 1,2-137 $ukb | grep -w '1262' | awk -F',' '{print $1}') #Won't look at eIDs since they are all 7-digit

#xsv headers $ukb | grep '\s41270-' (ICD10) gives 138
ICD10=$(xsv select 1,138 $ukb | grep -w 'G20' | awk -F',' '{print $1}') #-w considers G20|X12 as separate words

#xsv headers $ukb | grep -E '\s40001-|\s40002-' (1ry & 2ry deaths) gives 200-231
death=$(xsv select 1,200-231 $ukb | grep -w 'G20' | awk -F',' '{print $1}') 

#Merge
PD=$(echo "$self_report" "$ICD10" "$death" | tr ' ' '\n' | grep -v "eID" | sort -u) #Using "$_" preserves /n separators BUT when concatenating, echo adds spaces in b/w (NOT newlines)
echo "$PD" > PD.txt



## Early onset PD
#old: xsv headers $ukb | grep -E '\s34-|\s42032-' (YearOfBirth, algorithm PD report date) gives 194,234; notes: 1900-01-01 represents "Date is unknown" (turns out that no one has this in this case)
# PD=$(xsv select 1,194,234 $ukb | awk -F',' '$3 != ""' | grep -v 1900-01-01 | awk -F',' 'FNR !=1{print $1}')
# echo "$PD" > PD.txt

#Compute age at diagnosis and get early-onset; notes: edge case from rounding down to YYYY= turns 51 on that year but had Dx before birthday (that's why <=51 and not 50)
EOPD=$(xsv select 1,194,234 $ukb | awk -F',' '$3 != ""' | grep -v 1900-01-01 | awk -F',' 'FNR !=1 {split($3, date_parts, "-");print $1 "," $2 "," date_parts[1]}'| \
  awk -F',' '{age = $3 - $2;if (age <= 51) print $1}')

echo "$EOPD" > EOPD.txt


## Proxies
#Select proxy cases from father, mother, sibling, and exclude PD cases

#xsv headers $ukb | grep -E '\s20111-|\s20110-|\s20107-' gives 139-150
all_proxy=$(xsv select 1,139-150 $ukb | grep -w -P '(?<!-)11(?!\d)' | awk -F',' '{print $1}'| sort -u) #"perl" regex to only select 11 (and not -11) See UKBB glossary for possible values
proxy=$(echo "$all_proxy" | grep -v -F -x -f <(echo "$PD")) #-F means no regex, -x is whole line match, -f feeds file of patterns
echo "$proxy" > proxy.txt


## Controls
cases=$(echo "$PD" "$proxy" |tr ' ' '\n')
ctrl=$(xsv select 1 $ukb | grep -v -F -x -f <(echo "$cases"))
echo "$ctrl" > ctrl.txt

#---------------------------------------------------------------------------------------------

Rscript sampleIDs.R $ukb


#For Log
mkdir -p logs
echo "PD cases: $(wc -l < PD.txt)" > logs/sampleIDs.log
echo "EOPD cases: $(wc -l < EOPD.txt)" >> logs/sampleIDs.log
echo "PD proxies: $(wc -l < proxy.txt)" >> logs/sampleIDs.log
echo "Controls: $(wc -l < ctrl.txt)" >> logs/sampleIDs.log

printf "\nPost sample-level QC: \n" >> logs/sampleIDs.log
echo "PD cases: $(($(wc -l < ukbb_PD_case_covar.txt) -1))" >> logs/sampleIDs.log
echo "EOPD cases: $(($(wc -l < ukbb_EOPD_case_covar.txt) -1))" >> logs/sampleIDs.log
echo "PD proxies: $(($(wc -l < ukbb_PD_proxy_covar.txt) -1))" >> logs/sampleIDs.log
echo "Controls: $(($(wc -l < ukbb_PD_control_covar.txt) -1))" >> logs/sampleIDs.log

# #Clean up
# rm PD.txt proxy.txt ctrl.txt EOPD.txt

#----------------------------------------------------------------------------------------------
#Merge for WGS extraction (from UKB_RAP platform) while skipping headers

awk -F',' 'FNR==1 {next} {print $1}' \
  ukbb_PD_case_covar.txt \
  ukbb_PD_proxy_covar.txt \
  ukbb_PD_control_covar.txt | sort -u > Sajanth_sampleIDs
