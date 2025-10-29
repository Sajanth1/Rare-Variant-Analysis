#!/bin/bash

# Run from import

# index the file
# concatenate them
module load StdEnv/2020 gcc/9.3 bcftools/1.16

#for file in WGS/*.reheadered.vcf.gz;do
#    command="tabix $file"
#    if [ ! -f $file.tbi ];then
#        sbatch -c 2 --mem=30g -t 00:03:00 --out logs/$file.tabix.out --account=def-grouleau --wrap "$command"
#    else 
#        echo "index file found, skipped"
#    fi
#done


#------------------------------------------------------------------
# Index final vcf
command="tabix UKBB.vcf.gz"
sbatch -c 2 --mem=16g -t 6:00:00 --out logs/UKBB.vcf.gz.tabix --account=def-grouleau --wrap "$command"
