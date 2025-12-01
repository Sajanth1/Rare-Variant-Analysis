#!/bin/bash

# Add Status columns to each group && choose which subsets are used

# awk 'BEGIN{FS=OFS=","} FNR==1{print $0, "Status"} FNR!=1{print $0, 2}' ukbb_PD_case_covar.txt > tmp && mv tmp ukbb_PD_case_covar_eval.txt # PD cases
awk 'BEGIN{FS=OFS=","} FNR==1{print $0, "Status"} FNR!=1{print $0, 2}' ukbb_EOPD_case_covar.txt > tmp && mv tmp ukbb_EOPD_case_covar_eval.txt # EOPD cases
#awk 'BEGIN{FS=OFS=","} FNR==1{print $0, "Status"} FNR!=1{print $0, 2}' ukbb_PD_proxy_covar.txt > tmp && mv tmp ukbb_PD_proxy_covar_eval.txt # PD proxies
awk 'BEGIN{FS=OFS=","} FNR==1{print $0, "Status"} FNR!=1{print $0, 1}' ukbb_PD_control_covar.txt > tmp && mv tmp ukbb_PD_control_covar_eval.txt # controls

# Merge selected covariates into one file
awk 'BEGIN{FS=OFS=","} NR==1{print $0} FNR!=1{print $0}' *_eval.txt > covar_UKBB.txt
tr ',' '\t' < covar_UKBB.txt > tmp && mv tmp covar_UKBB.txt

# Add FID & IID columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $0}' covar_UKBB.txt > tmp
awk 'BEGIN{FS=OFS="\t"} FNR==1 {$1 = "FID" ; $2 = "IID"; print} FNR!=1{print}' tmp > covar_UKBB.txt ; rm tmp

# Create Sex file for plink
awk 'BEGIN{FS=OFS="\t"} FNR==1{print "FID", "IID", "Sex"} FNR!=1{print $1, $2, $7}' covar_UKBB.txt > sex_UKBB.txt