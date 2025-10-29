#!/bin/bash
# Use perl virtual environment
# virtual environment path: ~/soft/packages/ensembl-vep_env/perl_local_lib
# For more options running Vep see:
# https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html

#Run from scratch/dir
cohort=$1
VCF="${cohort}.vcf.gz"
OUTPUT=./annotation/$cohort.VEPannotated

# Load dependencies
# module load StdEnv/2023 perl/5.36.1 gcc/13.3 mariadb/11.5.0 htslib/1.22.1 db/18.1.32
perl -Mlocal::lib=~/links/projects/def-grouleau/COMMON/soft/packages/ensembl-vep_env/perl_local_lib
eval $(perl -Mlocal::lib=~/links/projects/def-grouleau/COMMON/soft/packages/ensembl-vep_env/perl_local_lib)


# Setting environmental variables
export VEP_DIR=~/links/projects/def-grouleau/COMMON/soft/packages/ensembl-vep_env/ensembl-vep
CACHE=${VEP_DIR}/.vep/cache
PLUGINS_DIR=${VEP_DIR}/.vep/Plugins
FASTA=${CACHE}/homo_sapiens/114_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

#ClinVar custom annotation
CLINVAR_VCF=~/links/projects/def-grouleau/COMMON/data/VEP/clinvar/clinvar_20240407.GRCh38.vcf.gz
CLINVAR_FIELDS="ALLELEID,CLNDN,CLNREVSTAT,CLNSIG,CLNSIGCONF"

#Plugin data files
ALPHAMISSENSE_DATA=${PLUGINS_DIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz
CADD_SNV=${PLUGINS_DIR}/CADD/whole_genome_SNVs.tsv.gz # This is version 1.7
CADD_INDEL=${PLUGINS_DIR}/CADD/gnomad.genomes.r4.0.indel.tsv.gz


#Run VEP
$VEP_DIR/vep -i ${VCF} -o ${OUTPUT} --force_overwrite --safe \
    --fork 4 --cache $CACHE --buffer_size 5000 \
    --format vcf --offline  --fasta ${FASTA} --xref_refseq --assembly GRCh38 \
    --pick_allele \
    --show_ref_allele --af_gnomadg --hgvs --symbol \
    --custom ${CLINVAR_VCF},ClinVar,vcf,exact,0,${CLINVAR_FIELDS} \
	--dir_plugins ${PLUGINS_DIR} --plugin AlphaMissense,file=${ALPHAMISSENSE_DATA} --plugin CADD,${CADD_SNV},${CADD_INDEL}


# Deactivate perl environment
eval $(perl -Mlocal::lib=--deactivate,~/links/projects/def-grouleau/COMMON/soft/packages/ensembl-vep_env/perl_local_lib)
