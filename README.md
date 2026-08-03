# Rare-Variant-Analysis

Run gene-based rare variant analysis on genomics data from large biobanks (here, AMP-PD and UKBB). This branch is specifically for the GP2 cohort. Relevant paper: [INSERT]

<em>For any questions: sajanth.kanagasingam@mail.mcgill.ca</em> (I would be more than happy to help!)
<br/> <br/>

SECTIONS:
Covar - prepare your covariate files here <br/>
Import - process raw vcf data here <br/>
Analysis - all filtration and analysis happens here! <br/>

*Note: package management is done via mamba and apptainer. Also, unless otherwise stated, all scripts are run with working directory set to "import", "covar", or "analysis" depending on which part of the workflow you are operating in. 

<br/><br/>

## 1. Covar

Covariates file MUST be named in the following format: "covar_COHORT.txt" where COHORT is the name of your dataset (eg. covar_UKBB.txt). 

If generating covar under same conditions as we did, simply run merge_covar_GP2.sh to create covar files for all ancestries.


## 2. Import
Before starting a project that uses WGS data, one must always create a bed file that contains the genomic coordinates of the genes to be studied. A template ("import/Sajanth44.GRCh38.bed") is provided to aid you in the process of creating one. Your bed file must be formatted in the same way (CHR START END ENST GENE_NAME)
<br/><br/>

For the WGS data itself, the main analysis script ("analysis/rare_variant.VEP.sh") has two file-type options to analyse pre-processed WGS data. Feel free to mix-and-match, but the provided script assumes that bfiles are provided for AMP_PD and a vcf.gz file is provided for UKBB and GP2.

Option 1: bfiles (ie. plink1.9 format or plink2 with --make-bed). 

Option 2: vcf.gz. If you have multiple vcf.gz, please index, sort and merge them beforehand. Template index and merge scripts are provided in "import".
<br/><br/>

*Note: a template preprocessing script for vcf.gz files is provided (preprocess_gp2.sh). Methodology is further described in the paper associated with this repo.
<br/><br/>




## 3. Analysis

As simple as running rare_variant.VEP.sh after modifying "USER INPUT" section with your paths! Make sure to have VEP installed beforehand (suggestion: run all downloads beforehand then comment out download commands)!

Now, this pipeline has quite a bit of flexibility in that you can directly modify the variant sets and their definitions (vep_setid_prep.py line 146+), the covariates to be included in SKAT-O (SKATO.r line 48), p-value correction method (FDR.r line 25), etc.
<br/><br/>
