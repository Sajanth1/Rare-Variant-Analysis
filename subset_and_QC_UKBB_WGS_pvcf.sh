#!/bin/bash
# Version: 25/06/2025 (no longer asks for genelist and does not attempt download)


# if you wanna provide a bed file. name it with the format ${pathname}.GRCh38.bed
#gene_list is simply a text file with gene name on each line
# gene_list txt file of gene names on each line
# output output folder for bed files
# pathname name of the task
# keep_samples, sample file with sample ID on each line 

# dx upload sex_for_plink.txt --path temp/


# Run from import/script
read output pathname keep_samples <<< $@
echo "output:$output; pathway:$pathname; keep_samples:$keep_samples"
echo "Do not forget to upload bed and sampleID files to dna nexus (in temp)"

if [[ -z "${pathname}.GRCh38.bed" ]]; then
  echo "Error: bed file is missing."
fi

if [[ -z "$output" ]]; then
  echo "Error: output is missing."
fi

if [[ -z "$pathname" ]]; then
  echo "Error: pathname is missing."
fi

if [[ -z "$keep_samples" ]]; then
  echo "Error: keep_samples is missing. you have to submit a sample ID file. this is mandatory in the latest update."
fi


module load StdEnv/2020 scipy-stack/2020a python/3.8.10
# before you run this command
## make sure your dx is setup
## check this to select project https://documentation.dnanexus.com/user/projects/project-navigation
# dx select project-GvFxJ08J95KXx97XFz8g2X2g
export PATH="$HOME/.local/bin:$PATH"

echo "confirm the bed file is ok"
echo "continue?(y/n)"
read answer 
if [[ "$answer" == "y" ]];then 
    # identify batch files covering the genes based on the bed file
    if [ ! -f ${pathname}.batch.txt ];then
        echo "identifying batch files associated with the genes"
        bash identify_batch_files.sh $pathname.GRCh38.bed
    else 
        echo "batch file found"
    fi
    echo "confirm the batch file is ok"
    echo "continue?(y/n)"
    read answer
    if [[ "$answer" == "y" ]];then 
        dx mkdir -p genes/$pathname/
        out=genes/$pathname/
        echo "bed file:$pathname.GRCh38.bed $out batch file:${pathname}.batch.txt, keep_samples:$keep_samples"
        bash call_variant_WGS.sh $pathname.GRCh38.bed ${pathname}.batch.txt $keep_samples $out 
        
        # the logic here can be adjusted. it is somewhat confusing (basically exits call_variant_WGS.sh once submits first job as test. After, reruns same script but the while loop runs in that script runs instead)
        echo "if you chose yes and submitted one job to test. please enter y to rerun the program after you confirm the result (y/n), otherwise, type n to quit the program"
        read answer
        if [[ "$answer" == "y" ]];then 
            bash call_variant_WGS.sh $pathname.GRCh38.bed ${pathname}.batch.txt $keep_samples $out
        else
            echo "quit reruning. you need to look into log on UKB RAP and see what is going on. check your bed and sample ID file if anything went wrong"
            exit 1
        fi
    else 
        echo "batch file not good, check your bed files and re-run identify_batch_files.sh"
        exit 42
    fi 
    
else
    echo "bed file not good, check your gene names and re-run this"
    rm $pathname.GRCh38.bed
    exit 42
fi 


# History
# Cem's 418 genes (ML predictions)
# bash subset_and_QC_UKBB_WGS_pvcf.sh . Cem418 Sajanth_sampleIDs

# UPDATE: no longer asks for genelist nor does it attempt download

# Sajanth's 13 genes (Dystonia-PD)
# bash call_variant_UKBB_WGS_pvcf.part1.sh Sajanth13.genelist.txt . Sajanth13 ~/scratch/UKBB/WGS_partial/import/WGS Sajanth_sampleIDs

# Sajanth's 44 genes (full Dystonia-PD)
# bash call_variant_UKBB_WGS_pvcf.part1.sh Sajanth44.genelist.txt . Sajanth44 ~/scratch/UKBB/WGS_full/import/WGS Sajanth_sampleIDs
