# Rare-Variant-Analysis (IN-PROGRESS)
Run rare variant analysis on genomics data from large biobanks (here, AMP-PD and UKBB). This pipeline includes pre-processing of bfiles or vcfs, as well as single-variant analysis and meta-analysis.


Please organize your directory tree in the following manner (so that the relative paths work):
![schema1](https://github.com/Sajanth1/Rare-Variant-Analysis/blob/main/Media/Schema1.png)

Import - download raw vcf or bfile data here <br/>
Covar - prepare your covariate file here <br/>
Analysis - all processing and analysis happens here! <br/>


## 1. Import
Option 1: vcf.gz. If you have multiple vcf.gz, please order, index and merge them beforehand. A template merge script is provided.

Option 2: bfiles (ie. plink1.9 format or plink2 with --make-bed). 

## 2. Covar



## 3. Analysis
