# Initial half of this script can be written in xsv to save compute
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
ukb_path <- args[1]

PD <- readLines("PD.txt")
PD_proxy <- readLines("proxy.txt")
PD_control <- readLines("ctrl.txt")
EOPD <- readLines("EOPD.txt")

#-----------------------------------------------------------------------------------
#Read covar table
ukb <- as.data.frame(fread(ukb_path))

#Field list
field <- c("20002","41270","20111","20110","20107","22001","22006","22009","22000","34","21022","22021","22019","22027", "22189", "40001")

#Change into pattern recognisable by grep  
pattern <- paste0("^",field,"-",collapse = "|")

#Select fields from overall covar
ukb_filtered <- ukb[,c(1,grep(pattern,names(ukb)))]

#-----------------------------------------------------------------------------------


## Sample-level QC

#Perform filter for samples with known issue (aneupleudy, missingness, het outlier) and relatedness (0 = no closer than 3rd degree relative) & ancestry filter (1 = causacian)
unrelated <- readLines("~/runs/go_lab/GRCh37/ukbb/ukbb_raw_data_no_cousins.txt")

ukb_filtered_unrelated <- ukb_filtered[ukb_filtered$eid %in% unrelated,]
ukb_filtered_unrelated_euro <- ukb_filtered_unrelated[ukb_filtered_unrelated$"22006-0.0" %in% 1,]
ukb_filtered_unrelated_euro_aneu <- ukb_filtered_unrelated_euro[!(ukb_filtered_unrelated_euro$"22019-0.0" %in% 1),]
ukb_filtered_unrelated_euro_aneu_miss <- ukb_filtered_unrelated_euro_aneu[!(ukb_filtered_unrelated_euro_aneu$"22027-0.0" %in% 1),]


## Get PCs
PC <- as.data.frame(fread("~/runs/go_lab/GRCh37/ukbb/pc_euro.txt"))
PC$IID <- NULL
names(PC)[1] <- "eid"


## Select required columns and add PCs

#Select covariates used in GWAS
covar_field <- c("34","22189","21022","22000","22001")# PC is not included here
covar <- paste0("^",covar_field,"-",collapse = "|")
ukb_covar <- ukb_filtered_unrelated_euro_aneu_miss[,c(1,grep(covar,names(ukb_filtered_unrelated_euro_aneu_miss)))]

ukb_covar_pc10 <- merge(ukb_covar[,c("eid", "34-0.0", "22189-0.0", "21022-0.0", "22000-0.0", "22001-0.0")], PC, by="eid") #merge based on eID
#Rename covariates; fixed order
names(ukb_covar_pc10) <- c("ID", "YearAtBirth", "Townsend", "Age", "Batch", "Sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10")
ukb_covar_pc10$Sex = ifelse(ukb_covar_pc10$Sex == 0, 2, ukb_covar_pc10$Sex) # to keep it consistent with plink1.9 format

## Save groups
write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% PD,], "ukbb_PD_case_covar.txt", quote = F, row.names = F)
write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% EOPD,], "ukbb_EOPD_case_covar.txt", quote = F, row.names = F)
write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% PD_proxy,], "ukbb_PD_proxy_covar.txt", quote = F, row.names = F)
write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% PD_control,], "ukbb_PD_control_covar.txt", quote = F, row.names = F)
