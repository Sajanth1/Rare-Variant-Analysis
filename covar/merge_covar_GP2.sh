#!/bin/bash
#---------
# Version (/05/2026): managing new age definition
#
#--------

#mamba create -n covar python=3.12.3 pandas numpy pathlib
eval "$(mamba shell hook --shell bash)" ; mamba activate covar


# Generate ancestry-specific covars
#wb resource mount
python3 covar_prep_GP2.py &> covar.log

# Select covars (other ancestries: AAC FIN)
for ancestry in "AFR" "AJ" "AMR" "CAS" "EAS" "EUR" "MDE" "SAS" "CAH"; do

    # Remove controls w/ AAO or Dx as age source
    awk 'BEGIN{OFS="\t"} NR==1{print; next} ( $4 >=3 || $17 ==2 ) {print $0}' ${ancestry}_covar_GP2.txt > tmp && mv tmp ${ancestry}_covar_GP2.txt  # See covar_prep_GP2.py for mapping


    #-----------OPTIONS------------------

    # # Select EOPD
    #awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ( ($4 =="age_of_onset" && $5 <=51) || $17 ==1 ) {print $0}' ${ancestry}_covar_GP2.txt > tmp && mv tmp ${ancestry}_covar_GP2.txt

    # Remove samples w/ missing age covars
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($5!="NA") {print}' ${ancestry}_covar_GP2.txt > tmp && mv tmp ${ancestry}_covar_GP2.txt

    # Remove AMP-PD samples
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($3!="1") {print}' ${ancestry}_covar_GP2.txt > tmp && mv tmp ${ancestry}_covar_GP2.txt

    # Remove Brainbank samples?

    #-----------------------------------

    awk 'BEGIN{FS=OFS="\t"} FNR!=1{print $1}' ${ancestry}_covar_GP2.txt > ${ancestry}_sampleIDs_GP2.txt #samples to extract for each ancestry

done


