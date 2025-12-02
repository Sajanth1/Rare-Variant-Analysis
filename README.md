# Rare-Variant-Analysis

Run gene-based rare variant analysis on genomics data from large biobanks (here, AMP-PD and UKBB). This pipeline includes pre-processing of bfiles or vcfs, as well as single-variant analysis and meta-analysis. 

<em>For any questions: sajanth.kanagasingam@mail.mcgill.ca</em> (I would be more than happy to help!)
<br/> <br/>


Please organize your directory tree in the following manner (so that the relative paths work):
![schema1](https://github.com/Sajanth1/Rare-Variant-Analysis/blob/main/Media/Structure.png)

Import - download raw vcf or bfile data here <br/>
Covar - prepare your covariate file here <br/>
Analysis - all processing and analysis happens here! <br/>

*Note: package management is done via module load (standard for shared compute clusters). Also, unless otherwise stated, all scripts are run with working directory set to "import", "covar", or "analysis" depending on which part of the workflow you are operating in. 

*Jobs will be submitted using slurm. If you wish to run these scripts locally, simply remove "sbatch" wrapper around some of the commands.
<br/><br/>


## 1. Import
The main analysis script ("rare_variant.VEP.sh") has two file-type options to analyse pre-processed WGS data. Feel free to mix-and-match, but the provided script assumes that bfiles are provided for AMP_PD and a vcf.gz file is provided for UKBB.

Option 1: vcf.gz. If you have multiple vcf.gz, please index, sort and merge them beforehand. Template index and merge scripts are provided in "import".

Option 2: bfiles (ie. plink1.9 format or plink2 with --make-bed). 
<br/>

*Note: a template preprocessing script for vcf.gz files is provided (preprocess_WGS.sh). Methodology is further described in the paper associated with this repo.
<br/><br/>


## 2. Covar

Covariates file MUST be named in the following format: "covar_COHORT.txt" where COHORT is the name specified above. 

It is assumed that your covar has at least the following columns: "FID IID Sex Age PC1 PC2 PC3 PC4 PC5 Status", where Status is coded as 1=control, 2=case. If not, update SKATO.r script accordingly (in the analysis section). An example covar (for UKBB) is also provided. 
<br/><br/>


## 3. Analysis

As simple as running rare_variant.VEP.sh after modifying "USER INPUT" section with your paths! Make sure to have VEP installed beforehand and to update the paths in VEP114_annotation.sh.

Now, there is also quite a bit of flexibility in that you can directly modify the variant sets and their definitions (vep_setid_prep.py line 146+), the covariates to be included in SKAT-O (SKATO.r line 48), p-value correction method (FDR.r line 25), etc.
<br/><br/>


## 4. Meta-Analysis
*metaSKAT.sh should be run from a dir that contains all cohorts as sub-directories.

Each set's variants will be aggregated across cohorts and meta-analyzed. As usual, provide appropriate info in "USER INPUT" section of metaSKAT.sh. Here's a schema to help you visualize the meta-analysis process.
![schema2](https://github.com/Sajanth1/Rare-Variant-Analysis/blob/main/Media/Meta_Analysis.png)

*Note: one must maintain the directory architecture created by this pipeline for each cohort in order for the scripts to be able to fetch the appropriate files.
