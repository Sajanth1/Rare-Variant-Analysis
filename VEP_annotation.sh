#!/bin/bash
# Note: Currently apptainer fetches version 115 (GENCODE V49)

#Run from scratch/dir
cohort=$1
VCF="${cohort}.vcf.gz"
OUTPUT=./annotation/$cohort.VEPannotated


## Install VEP
#apptainer pull --name vep.sif docker://ensemblorg/ensembl-vep
#mkdir $HOME/vep_data
#apptainer exec vep.sif INSTALL.pl -c $HOME/vep_data -a cf -s homo_sapiens -y GRCh38 # Get cache and fasta

# Setting environmental variables
VEP_DIR=$HOME/vep_data
FASTA=${VEP_DIR}/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
PLUGINS_DIR=${VEP_DIR}/Plugins


# Install Clinvar (custom annotation); this gets latest version (s.t. don't have to keep updating link but check version for paper)
#curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi --output-dir "$VEP_DIR"
CLINVAR_VCF=${VEP_DIR}/clinvar.vcf.gz
CLINVAR_FIELDS="ALLELEID,CLNDN,CLNREVSTAT,CLNSIG,CLNSIGCONF"

# Install Plugins
mkdir -p "${PLUGINS_DIR}/CADD"

CADD_SNV=${PLUGINS_DIR}/CADD/whole_genome_SNVs.tsv.gz
#wget -nv -P "${PLUGINS_DIR}/CADD" https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz #CADD 1.7 
#wget -nv -P "${PLUGINS_DIR}/CADD" https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi

CADD_INDEL=${PLUGINS_DIR}/CADD/gnomad.genomes.r4.0.indel.tsv.gz
#wget -nv -P "${PLUGINS_DIR}/CADD" https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz #CADD 1.7 (indels)
#wget -nv -P "${PLUGINS_DIR}/CADD" https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz.tbi

ALPHAMISSENSE=${PLUGINS_DIR}/AlphaMissense_hg38.tsv.gz
#wget -nv -P $PLUGINS_DIR https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz #AM (August 3, 2023)
#tabix -s 1 -b 2 -e 2 -f -S 1 "$ALPHAMISSENSE"



#Run VEP
apptainer exec --bind $HOME --pwd $HOME/analysis vep.sif \
    vep -i ${VCF} -o ${OUTPUT} --dir ${VEP_DIR} --force_overwrite --safe \
        --fork 4 --cache --buffer_size 5000 \
        --format vcf --offline  --fasta ${FASTA} --xref_refseq --assembly GRCh38 \
        --pick_allele \
        --show_ref_allele --af_gnomadg --hgvs --symbol \
        --custom ${CLINVAR_VCF},ClinVar,vcf,exact,0,${CLINVAR_FIELDS} \
        --dir_plugins ${PLUGINS_DIR} --plugin AlphaMissense,file=${ALPHAMISSENSE} --plugin CADD,${CADD_SNV},${CADD_INDEL}

