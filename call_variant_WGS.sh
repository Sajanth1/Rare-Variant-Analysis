#!/bin/bash
# example: bash call_variant_WGS.sh COMMANDER23.GRCh38.bed COMMANDER23.batch.txt  /lustre03/project/6004655/COMMUN/runs/lang/Commander/data/Filter_ID.txt genes/COMMANDER23_QC
bed_file_local=$1
batch=$2
sample_file_local=$3
destination_folder=$4
project_src="project-GvFxJ08J95KXx97XFz8g2X2g"
script_file=SK_preprocess_WGS.sh
bed_file1=$(basename $bed_file_local)
sample_file=$(basename $sample_file_local)
instance_type=mem2_ssd1_v2_x8 #default is mem1_ssd1_v2_x4 you can do mem1_ssd1_v2_x8 if there is an oom killed error (enough for most, largest require mem2). or check this page for other instance https://20779781.fs1.hubspotusercontent-na1.net/hubfs/20779781/Product%20Team%20Folder/Rate%20Cards/BiobankResearchAnalysisPlatform_Rate%20Card_Current.pdf
bed_file_name=$(basename $bed_file1 .GRCh38.bed)

if dx ls "temp/${bed_file1}" > /dev/null 2>&1; then
    echo "bed file found on UKB RAP"
else 
    echo "uploading the bed file"
    dx upload $bed_file_local --path temp/
fi
if dx ls "temp/${sample_file}" > /dev/null 2>&1; then
    echo "sample file found on UKB RAP"
else 
    dx upload $sample_file_local --path temp/
fi
first_file=$(head -n 1 "$batch")
file_name=$(basename "$first_file" .vcf.gz)
output_file="${file_name}.subset.DP25.GQ25.MISS95.split.vcf.gz"
if dx ls "${destination_folder}/${output_file}" > /dev/null 2>&1; then
  echo "Output for first job $output_file already exists. Skipping first line."
  echo "we are going to submitting jobs for the remaining files"
else
  read -p "Do you want to submit ONE job to test? (y/n): " choice
  if [[ "$choice" == "y" ]]; then
    # Submit only the first job
    first_file=$(head -n 1 "$batch")
    file_name=$(basename "$first_file")
    job_name="${bed_file_name}.${file_name}"
    dx run app-swiss-army-knife \
      --instance-type "$instance_type" \
      --priority high \
      -y \
      -iin="${first_file}" \
      -iin="${first_file}.tbi" \
      -iin="${project_src}:temp/sex_for_plink.txt" \
      -iin="${project_src}:temp/${sample_file}" \
      -iin="${project_src}:temp/${bed_file1}" \
      -iin="${project_src}:scripts/${script_file}" \
      --name="${job_name}" \
      --destination="${destination_folder}" \
      --brief \
      -icmd="bash ${script_file} ${file_name} ${sample_file} ${bed_file1}"

    exit 0
  fi
fi

while read -r file; do
  file_name=$(basename "$file" .vcf.gz)
  output_file="${file_name}.subset.DP25.GQ25.MISS95.split.vcf.gz"
  file_name=$(basename "$file")
  if ! dx ls "${destination_folder}/${output_file}" > /dev/null 2>&1; then
    job_name="${bed_file_name}.${file_name}"
    dx run app-swiss-army-knife \
      --instance-type "$instance_type" \
      --priority high \
      -y \
      -iin="${file}" \
      -iin="${file}.tbi" \
      -iin="${project_src}:temp/sex_for_plink.txt" \
      -iin="${project_src}:temp/${sample_file}" \
      -iin="${project_src}:temp/${bed_file1}" \
      -iin="${project_src}:scripts/${script_file}" \
      --name="${job_name}" \
      --destination="${destination_folder}" \
      --brief \
      -icmd="bash ${script_file} ${file_name} ${sample_file} ${bed_file1}"
  else 
    echo "output file ${output_file} found, skipped"
  fi 
done < "$batch"
