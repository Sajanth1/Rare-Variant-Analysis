#!/bin/bash
# Run from import

# index the files
module load StdEnv/2020 gcc/9.3 bcftools/1.16

#------------------------------------------------------------------
# For multiple vcfs in a directory:
dir="WGS"

for file in $dir/*.vcf.gz;do
   command="tabix $file"
   if [ ! -f $file.tbi ];then
       sbatch -c 2 --mem=30g -t 00:03:00 --out logs/$file.tabix.out --account=def-grouleau --wrap "$command"
   else 
       echo "index file found, skipped"
   fi
done


#------------------------------------------------------------------
## For single file:
name="UKBB.vcf.gz"

command="tabix $name"
sbatch -c 2 --mem=16g -t 6:00:00 --out logs/$name.tabix --account=def-grouleau --wrap "$command"
