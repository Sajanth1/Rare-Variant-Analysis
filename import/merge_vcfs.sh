#!/bin/bash
# run from import

module load StdEnv/2020 gcc/9.3.0 bcftools/1.16

## USER INPUT ##
vcf_folder="WGS"
cohort="UKBB"
out_name="${cohort}_preprocess.vcf.gz"
#--------------------------------

# Name your vcfs in such a way that they can be sorted (modify sort command as needed). Example:
ls ${vcf_folder}/*.vcf.gz | sort -t'_' -k2.2n -k3.2n > script/merge.list #numerical sort based on second (-c*) then third field (_b*) after header

echo "Here are the sorted batch files:"
cat script/merge.list
echo "Please confirm the sorting is correct (y/n):"
answer=y
read answer

if [[ "$answer" == "y" ]]; then 

    # Count the number of files
    file_count=$(wc -l < script/merge.list)

    # Time requested: 1 minute per file + 1 hour base
    time_minutes=$((file_count + 60))
    time_hours=$((time_minutes / 60))
    time_remaining_minutes=$((time_minutes % 60))
    time_requested=$(printf "%02d:%02d:00" $time_hours $time_remaining_minutes)

    # Fixed memory allocation
    mem_requested=30

    # Create the command with multi-threading
    command="bcftools concat --naive-force -f script/merge.list --threads 10 -Oz -o ${out_name} && echo 'concatenation done' && tabix ${out_name}"

    # Submit the job
    sbatch -c 2 --mem=${mem_requested}g -t ${time_requested} \
        --wrap "$command" \
        --account=def-grouleau \
        --output=logs/${cohort}_WGS_MERGE \
        --job-name=merge_vcfs

    echo "Submitted job for ${file_count} files, memory: ${mem_requested}G, time: ${time_requested}"
else
    echo "Redo the sorting for merge list and resubmit the job."
fi
