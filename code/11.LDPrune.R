library(data.table)
library(stringr)
library(parallel)
chr_list<-1:22
source('/home2/ardo/Code/Calling_UpQC_DownQC_Code/3.DownQC_source.R')
output_fol<-"/home2/kjkim/analysis/Dong/20200608_GWAS"
DATA_fol<-c("PH","MB","MG")
C3_fol<-"3.Result"
FOL_snp<-"SNP"
DATA_ind<-1
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])


MR_table_cbd<-read.csv(file.path(file.path(output_fol,C3_fol,"MB",FOL_snp,"TTT2_GWA_midoutput"),paste0("MR_table_ALLALL.csv")),check.names = F,stringsAsFactors = F)