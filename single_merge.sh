#!/bin/bash
# # run from import


# # concatenate files on DNAnexus

# ## 1. identify files in the target directory
# ## 2. create commands 
# ### --input files
# ### --output directory
# ### --use bcftools
# ### --create a concatenation text file of file names
# # Your DNAnexus project ID
# project_src="project-GvFxJ08J95KXx97XFz8g2X2g"
# module load StdEnv/2020 scipy-stack/2020a python/3.8.10
# export PATH="~/.local/bin/dx:$PATH"
# # The source folder *within your DNAnexus project* where the .vcf.gz files are stored
# vcf_src_folder=${project_src}:$1
# job_name=$(basename $vcf_src_folder)
# # The destination folder *within your DNAnexus project* where the final merged VCF will be saved
# destination_folder="/concat"
# # Name for the final output merged VCF file
# out_name="ukb23374_project_${job_name}.vcf.gz"
# # Name for the job (will be visible in the DNAnexus UI)
# merge_list=script/${job_name}_merge.list
# # Instance type for the job (you can change this based on needs, e.g., mem1_ssd1_v2_x4)
# instance_type=mem1_ssd1_v2_x8
# echo "instance type will be ${instance_type}, it can only work with files with total size < 200gb storage. if you are concatenating a large amount of files. please go to https://documentation.dnanexus.com/developer/api/running-analyses/instance-types and see which instance fits"


# # Get folder size: 
# # dx describe genes/Cem418/* | grep Size | awk '{ if ($3 == "MB") { total += $2 * 1024 * 1024 } else { total += $2 } } END { printf "Total: %.2f GB\n", total / 1024 / 1024 / 1024}'

# # =======================================================

# # --- 1. Create the local merge.list file ---
# echo "Creating merge.list from the target directory..."
# echo "this steps requires the vcf.gz files in the folder to have the format of ukb223374_c[chr]_b[batch_number]......"
# echo "if this is not the case, please create your own merge list file in import/script/ with name of [folder_name]_merge.list and run this script"
# # This assumes your local directory structure mirrors the project.
# # It finds all .vcf.gz files and writes their DNAnexus paths to merge.list.
# # If your local paths are different, you must construct the list differently.
# if [ ! -f ${merge_list} ];then 
# dx ls ${vcf_src_folder}/*.vcf.gz > ${merge_list}.temp
# cat "${merge_list}.temp" | sort -t'_' -k2.2,2n -k3.2,3n > ${merge_list}
# rm ${merge_list}.temp
# else 
# echo "merge list file found in script/${merge_list}, skipped"
# echo "you should confirm this merge list file is correct"
# fi 


# echo "Contents of merge list: (head)"
# head $merge_list
# echo "line number of merge list"
# echo "$(wc -l $merge_list)"

# # --- 2. Upload the merge.list file to the DNAnexus project ---
# echo "Uploading merge list to ${project_id}/concat..."
# if dx ls "${project_id}:concat/${job_name}_merge.list" > /dev/null 2>&1; then
#     echo "merge list file found on UKB RAP"
# else 
#     echo "uploading the merge list file"
#     dx upload $merge_list --path "${project_id}:concat/" --brief > /dev/null
# fi

# # The --brief flag outputs just the file ID, which we suppress with > /dev/null since we don't need it here.

# # --- 3. Build the dx run command ---
# echo "Constructing the dx run command..."

# # Start the base command
# icmd="dx run app-swiss-army-knife"

# # Add instance type and priority
# icmd="$icmd --instance-type \"$instance_type\""
# icmd="$icmd --priority high"

# # Add the INPUT for the merge.list file (uploaded in the previous step)
# icmd="$icmd -iin=\"${project_src}:concat/${job_name}_merge.list\""

# # Add an INPUT for EVERY VCF file listed in merge.list
# # This is crucial: the job needs explicit permission to access each file.
# while read vcf_path; do
#   icmd="$icmd -iin=\"${project_src}:$1/${vcf_path}\""
# done < $merge_list

# # Add job naming and destination
# icmd="$icmd --name=\"${job_name}_concatenation\""
# icmd="$icmd --destination=\"${project_src}:${destination_folder}\""

# # Add the iCMD with the bcftools concat command
# # The files are downloaded to the job's working directory.
# # Their names will be the same as their names in the project.
# # bcftools concat will read the list of paths from the merge.list file.
# icmd="$icmd -icmd=\"bcftools concat --naive-force -f ${job_name}_merge.list --threads 10 -Oz -o ${out_name} && tabix -p vcf ${out_name}\""

# # --- 4. Print and execute the command ---
# echo "Launching job with the following command:"
# echo "----------------------------------------"
# echo "$icmd"
# echo "----------------------------------------"


# echo "please check if you are concatenating the right folder, the files in the merge list file is in order and then confirm the number of files"
# echo "continue?(y/n)" 
# read answer 
# if [[ "$answer" == "y" ]]; then
# eval $icmd
# fi


# echo "Processing....Check the job status on DNANexus"
# echo "you can also watch the job by dx watch job-ID(job can be found when you submit, just scroll up"



#-------------------------- LOCAL -----------------------------
module load StdEnv/2020 gcc/9.3.0 bcftools/1.16

ls WGS/*.vcf.gz | sort -t'_' -k2.2n -k3.2n > script/merge.list #numerical sort based on second (-c*) then third field (_b*) after header

echo "Here are the sorted batch files:"
cat script/merge.list
echo "Please confirm the sorting is correct (y/n):"
answer=y
read answer

if [[ "$answer" == "y" ]]; then 
    # Define output file name
    out_name="UKBB_raw.vcf.gz"

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
        --output=logs/UKB_WGS_MERGE \
        --job-name=merge_vcfs

    echo "Submitted job for ${file_count} files, memory: ${mem_requested}G, time: ${time_requested}"
else
    echo "Redo the sorting for merge list and resubmit the job."
fi
