library(MendelianRandomization)
output_fol<-"/home2/kjkim/analysis/Dong/20200608_GWAS"
C1_fol<-"/KARE_QC"
C2_fol<-"2.DataReady"
C3_fol<-"3.Result"
FOL_snp<-"SNP"
log_binary<-F
DATA_fol<-c("PH","MB","MG")
DATA_ind<-1
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])

if(log_binary){
  MR_table_cbd_bf<-read.csv(file.path(file.path(output_fol,C3_fol,"MG",FOL_snp,"TTT2_GWA_midoutput"),paste0("logMR_table_ALL2MG.csv")),check.names = F,stringsAsFactors = F)  
}else{
  MR_table_cbd_bf<-read.csv(file.path(file.path(output_fol,C3_fol,"MG",FOL_snp,"TTT2_GWA_midoutput"),paste0("MR_table_ALL2MG.csv")),check.names = F,stringsAsFactors = F)  
}





.libPaths("~/tmp")
#install.packages("simex")
library(simex)

library(ieugwasr)
MR_table_cbd_sp<-split(MR_table_cbd_bf,MR_table_cbd_bf$'X-SNP')
TablePR<-MR_table_cbd_bf[,c("RS","X-SNP_pval","X-SNP")]
names(TablePR)<-c("rsid","pval","X-SNP")
TablePR_sp<-split(TablePR,TablePR$'X-SNP')
TablePRD_list<-lapply(TablePR_sp,function(df){
  output<-ld_clump(df,pop="EAS",clump_r2=0.01)
})

id_list<-lapply(TablePRD_list,"[[","rsid")

CBD_MR1_bf<-lapply(1:length(MR_table_cbd_sp),function(ind_data){
  MR_table_cbd_sp[[ind_data]][MR_table_cbd_sp[[ind_data]]$RS %in%id_list[[ind_data]],]
})
CBD_MR1<-do.call("rbind",CBD_MR1_bf)

DM_result<-read.csv(file.path(C3_PATH_OF,FOL_snp,"TTT_DM_WD_GWA_midoutput",paste0("DMWDresult_losenThre.csv")),check.names = F)
DM_res_readyORI<-DM_result[match(CBD_MR1$RS,DM_result$V2),c('Estimate','Std..Error','Pr...t..',paste0("PQL_",c('Estimate','Std..Error','Pr...t..')))]
names(DM_res_readyORI)<-c(paste0("DM1",c("_beta","_se","_pval")),paste0("DMPQL1",c("_beta","_se","_pval")))
MR_table_cbd_ORI<-cbind(CBD_MR1,DM_res_readyORI)

ind_PH_na<-which(is.na(MR_table_cbd_ORI$PH1_beta))
ind_MB_na<-which(is.na(MR_table_cbd_ORI$MB1_beta))
ind_MG_na<-which(is.na(MR_table_cbd_ORI$MG1_beta))

list_ind_NA<-list(ind_PH_na,ind_MB_na,ind_MG_na)


RS_find_list<-MR_table_cbd_ORI$RS[unique(do.call("c",list_ind_NA)[!is.na(do.call("c",list_ind_NA))])]
library(LDlinkR)
library(parallel)
tokens_bf<-c("002deef240a0","a691b83ea5ff","f79ecb07fdcf","fe549291c3a3","384269c403a4","2a6e728acf5d","e2442c6fc576","90eab294ac43","9051cb7a0903")
tokens<-c(tokens_bf,tokens_bf)
list_tmp.proxy1<-mclapply(1:9,function(ind_RS){
  print(ind_RS)
  cbind(LDproxy(RS_find_list[ind_RS],pop = "EAS", r2d = "r2", token = tokens[ind_RS], file = FALSE),ind_RS)
},mc.cores=9)

# ind_RS<-4
# list_tmp.proxy1[[ind_RS]]<-cbind(LDproxy(RS_find_list[ind_RS],pop = "EAS", r2d = "r2", token = tokens[ind_RS], file = FALSE),ind_RS)


list_tmp.proxy2<-mclapply(10:length(RS_find_list),function(ind_RS){
  print(ind_RS)
  cbind(LDproxy(RS_find_list[ind_RS],pop = "EAS", r2d = "r2", token = tokens[ind_RS], file = FALSE),ind_RS)
},mc.cores=9)
list_tmp.proxy<-c(list_tmp.proxy1,list_tmp.proxy2)



DATA_ind<-1
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
if(log_binary&DATA_ind==2){
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_log_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}else{
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}

RS_find_list_indiv<-MR_table_cbd_ORI$RS[list_ind_NA[[DATA_ind]]]
for (z in 1:length(RS_find_list_indiv)){
  ind_fill<-which(MR_table_cbd_ORI$RS==RS_find_list_indiv[z])
  z_proxy<-which(RS_find_list==RS_find_list_indiv[z])[1]
  
  # for(z_indi in z_proxy){
  tttt_RS<-list_tmp.proxy[[z_proxy]]$RS_Number[min(which(as.character(list_tmp.proxy[[z_proxy]]$RS_Number)%in%Result_cbd3sp[[1]]$V2))]
  # }
  
  
  for( u in 1:length(Result_cbd3sp)){
    tmp_res_tt <-Result_cbd3sp[[u]][which(Result_cbd3sp[[u]]$V2==tttt_RS),]
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_beta")] <- unlist(tmp_res_tt$Estimate)
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_se")] <- unlist(tmp_res_tt$'Std. Error')
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_pval")] <- unlist(tmp_res_tt$'Pr(>|t|)')  
  }
  print(z)
}

rm(Result_cbd3sp)



DATA_ind<-2
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
if(log_binary&DATA_ind==2){
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_log_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}else{
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}

RS_find_list_indiv<-MR_table_cbd_ORI$RS[list_ind_NA[[DATA_ind]]]
for (z in 1:length(RS_find_list_indiv)){
  ind_fill<-which(MR_table_cbd_ORI$RS==RS_find_list_indiv[z])
  z_proxy<-which(RS_find_list==RS_find_list_indiv[z])[1]
  tttt_RS<-list_tmp.proxy[[z_proxy]]$RS_Number[min(which(as.character(list_tmp.proxy[[z_proxy]]$RS_Number)%in%Result_cbd3sp[[1]]$V2))]
  
  for( u in 1:length(Result_cbd3sp)){
    tmp_res_tt <-Result_cbd3sp[[u]][which(Result_cbd3sp[[u]]$V2==tttt_RS),]
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_beta")] <- unlist(tmp_res_tt$Estimate)
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_se")] <- unlist(tmp_res_tt$'Std. Error')
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_pval")] <- unlist(tmp_res_tt$'Pr(>|t|)')  
  }
  print(z)
}

rm(Result_cbd3sp)
DATA_ind<-3
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
if(log_binary&DATA_ind==2){
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_log_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}else{
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}

RS_find_list_indiv<-MR_table_cbd_ORI$RS[list_ind_NA[[DATA_ind]]]
for (z in 1:length(RS_find_list_indiv)){
  ind_fill<-which(MR_table_cbd_ORI$RS==RS_find_list_indiv[z])
  z_proxy<-which(RS_find_list==RS_find_list_indiv[z])[1]
  tttt_RS<-list_tmp.proxy[[z_proxy]]$RS_Number[min(which(as.character(list_tmp.proxy[[z_proxy]]$RS_Number)%in%Result_cbd3sp[[1]]$V2))]
  
  for( u in 1:length(Result_cbd3sp)){
    tmp_res_tt <-Result_cbd3sp[[u]][which(Result_cbd3sp[[u]]$V2==tttt_RS),]
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_beta")] <- unlist(tmp_res_tt$Estimate)
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_se")] <- unlist(tmp_res_tt$'Std. Error')
    MR_table_cbd_ORI[ind_fill,paste0(DATA_fol[DATA_ind],u,"_pval")] <- unlist(tmp_res_tt$'Pr(>|t|)')  
  }
  print(z)
}

rm(Result_cbd3sp)

# sapply(list_tmp.proxy,function(data){
#   data$ind_RS
# })


# TablePRD_cbd_all<-TablePRD_cbd









DM_res_ready<-DM_result[match(MR_table_cbd_bf$RS,DM_result$V2),c('Estimate','Std..Error','Pr...t..',paste0("PQL_",c('Estimate','Std..Error','Pr...t..')))]
names(DM_res_ready)<-c(paste0("DM1",c("_beta","_se","_pval")),paste0("DMPQL1",c("_beta","_se","_pval")))

MR_table_cbd<-cbind(MR_table_cbd_bf,DM_res_ready)
head(MR_table_cbd)
MR_table_cbd2<-MR_table_cbd[!is.na(MR_table_cbd$DMPQL1_pval),]



MR_table_cbd2_sp<-split(MR_table_cbd2,MR_table_cbd2$'X-SNP')
TablePR2<-MR_table_cbd2[,c("RS","X-SNP_pval","X-SNP")]
names(TablePR2)<-c("rsid","pval","X-SNP")
TablePR2_sp<-split(TablePR2,TablePR2$'X-SNP')
TablePRD_list2<-lapply(TablePR2_sp,function(df){
  output<-ld_clump(df,pop="EAS",clump_r2=0.01)
})

TablePRD_cbd2<-do.call("rbind",TablePRD_list2)


id_list2<-lapply(TablePRD_list2,"[[","rsid")

CBD_MR1_bf2<-lapply(1:length(MR_table_cbd2_sp),function(ind_data){
  MR_table_cbd2_sp[[ind_data]][MR_table_cbd2_sp[[ind_data]]$RS %in%id_list2[[ind_data]],]
})
CBD_MR2<-do.call("rbind",CBD_MR1_bf2)

MR_table_cbd_DM<-CBD_MR2






ind_PH_na<-which(is.na(MR_table_cbd_DM$PH1_beta))
ind_MB_na<-which(is.na(MR_table_cbd_DM$MB1_beta))
ind_MG_na<-which(is.na(MR_table_cbd_DM$MG1_beta))

list_ind_NA<-list(ind_PH_na,ind_MB_na,ind_MG_na)


RS_find_list<-MR_table_cbd_DM$RS[unique(do.call("c",list_ind_NA)[!is.na(do.call("c",list_ind_NA))])]
library(LDlinkR)
library(parallel)
tokens_bf<-c("002deef240a0","a691b83ea5ff","f79ecb07fdcf","fe549291c3a3","384269c403a4","2a6e728acf5d","e2442c6fc576","90eab294ac43","9051cb7a0903")
tokens<-c(tokens_bf,tokens_bf)
list_tmp.proxy1<-mclapply(1:9,function(ind_RS){
  print(ind_RS)
  cbind(LDproxy(RS_find_list[ind_RS],pop = "EAS", r2d = "r2", token = tokens[ind_RS], file = FALSE),ind_RS)
},mc.cores=9)

# ind_RS<-4
# list_tmp.proxy1[[ind_RS]]<-cbind(LDproxy(RS_find_list[ind_RS],pop = "EAS", r2d = "r2", token = tokens[ind_RS], file = FALSE),ind_RS)


list_tmp.proxy2<-mclapply(10:length(RS_find_list),function(ind_RS){
  print(ind_RS)
  cbind(LDproxy(RS_find_list[ind_RS],pop = "EAS", r2d = "r2", token = tokens[ind_RS], file = FALSE),ind_RS)
},mc.cores=9)
list_tmp.proxy<-c(list_tmp.proxy1,list_tmp.proxy2)[1:length(RS_find_list)]



DATA_ind<-1
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
if(log_binary&DATA_ind==2){
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_log_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}else{
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}

RS_find_list_indiv<-MR_table_cbd_DM$RS[list_ind_NA[[DATA_ind]]]
for (z in 1:length(RS_find_list_indiv)){
  ind_fill<-which(MR_table_cbd_DM$RS==RS_find_list_indiv[z])
  z_proxy<-which(RS_find_list==RS_find_list_indiv[z])[1]
  
  # for(z_indi in z_proxy){
  tttt_RS<-list_tmp.proxy[[z_proxy]]$RS_Number[min(which(as.character(list_tmp.proxy[[z_proxy]]$RS_Number)%in%Result_cbd3sp[[1]]$V2))]
  # }
  
  
  for( u in 1:length(Result_cbd3sp)){
    tmp_res_tt <-Result_cbd3sp[[u]][which(Result_cbd3sp[[u]]$V2==tttt_RS),]
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_beta")] <- unlist(tmp_res_tt$Estimate)
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_se")] <- unlist(tmp_res_tt$'Std. Error')
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_pval")] <- unlist(tmp_res_tt$'Pr(>|t|)')  
  }
  print(z)
}


rm(RS_find_list_indiv)
rm(ind_fill)
rm(z_proxy)
rm(tttt_RS)
rm(Result_cbd3sp)



DATA_ind<-2
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
if(log_binary&DATA_ind==2){
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_log_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}else{
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}

RS_find_list_indiv<-MR_table_cbd_DM$RS[list_ind_NA[[DATA_ind]]]
for (z in 1:length(RS_find_list_indiv)){
  ind_fill<-which(MR_table_cbd_DM$RS==RS_find_list_indiv[z])
  z_proxy<-which(RS_find_list==RS_find_list_indiv[z])[1]
  tttt_RS<-list_tmp.proxy[[z_proxy]]$RS_Number[min(which(as.character(list_tmp.proxy[[z_proxy]]$RS_Number)%in%Result_cbd3sp[[1]]$V2))]
  
  for( u in 1:length(Result_cbd3sp)){
    tmp_res_tt <-Result_cbd3sp[[u]][which(Result_cbd3sp[[u]]$V2==tttt_RS),]
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_beta")] <- unlist(tmp_res_tt$Estimate)
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_se")] <- unlist(tmp_res_tt$'Std. Error')
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_pval")] <- unlist(tmp_res_tt$'Pr(>|t|)')  
  }
  print(z)
}
rm(RS_find_list_indiv)
rm(ind_fill)
rm(z_proxy)
rm(tttt_RS)
rm(Result_cbd3sp)


DATA_ind<-3
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
if(log_binary&DATA_ind==2){
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_log_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}else{
  Result_cbd3sp<-readRDS(file.path(C3_PATH_OF,"SNP/TTT2_GWA_midoutput",paste0("Result_cbd3sp_",DATA_fol[DATA_ind],".csv")))
}

RS_find_list_indiv<-MR_table_cbd_DM$RS[list_ind_NA[[DATA_ind]]]
for (z in 1:length(RS_find_list_indiv)){
  ind_fill<-which(MR_table_cbd_DM$RS==RS_find_list_indiv[z])
  z_proxy<-which(RS_find_list==RS_find_list_indiv[z])[1]
  tttt_RS<-list_tmp.proxy[[z_proxy]]$RS_Number[min(which(as.character(list_tmp.proxy[[z_proxy]]$RS_Number)%in%Result_cbd3sp[[1]]$V2))]
  
  for( u in 1:length(Result_cbd3sp)){
    tmp_res_tt <-Result_cbd3sp[[u]][which(Result_cbd3sp[[u]]$V2==tttt_RS),]
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_beta")] <- unlist(tmp_res_tt$Estimate)
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_se")] <- unlist(tmp_res_tt$'Std. Error')
    MR_table_cbd_DM[ind_fill,paste0(DATA_fol[DATA_ind],u,"_pval")] <- unlist(tmp_res_tt$'Pr(>|t|)')  
  }
  print(z)
}
rm(RS_find_list_indiv)
rm(ind_fill)
rm(z_proxy)
rm(tttt_RS)
rm(Result_cbd3sp)


  

  
mb_label<-c("X3.Hydroxybutyrate", "Acetoacetate")
mg_label<-c("C34", "R34")
ph_label_bf<-c("HBA1C_log","GLU0_TR_log","GLU60_TR_log","GLU120_TR_log","INS0_log","INS60_log","INS120_log","BMI_log","WHR")
# ind_pheno_go<-setdiff(1:length(ph_label_bf),c(3,5,6,8,9))
ind_pheno_go<-setdiff(1:length(ph_label_bf),c())
ph_label<-c("HBA1C_log","GLU0_TR_log","GLU60_TR_log","GLU120_TR_log","INS0_log","INS60_log","INS120_log","BMI_log","WHR")[ind_pheno_go]
dm_label<-c("DM_WD")

MB_list<-paste0("MB",1:2)
MG_list<-paste0("MG",1:length(mg_label))
PH_list<-paste0("PH",1:length(ph_label_bf))[ind_pheno_go]
DM_list<-paste0("DM",1:length(dm_label))

result_cbd_MBMG<-data.frame(matrix(nrow=0,ncol=7))

list1<-MB_list
list1_nm<-mb_label
list2<-PH_list
list2_nm<-ph_label

rm(list=c("list1","list2","list1_nm","list2_nm"))
add_CIs<-function(seq){
  c(seq[1:2],seq[1]-qnorm(0.975)*seq[2],seq[1]+qnorm(0.975)*seq[2],seq[3])
}
n_perm<-50000

n_cores<-32
# rm(list=c("list1","list2","list1_nm","list2_nm"))
library(parallel)
conduct_MR_WM<-function(clist1,clist1_nm,clist2,clist2_nm,cdata,crev=T){
  t_table<-matrix(nrow=length(clist1)*length(clist2),ncol=3)
  cnt<-0
  t_table[,1]<-1:dim(t_table)[1]
  for ( i in 1:length(clist1)){
    for ( j in 1:length(clist2)){
      cnt<-cnt+1
      
      t_table[cnt,2]<-i
      t_table[cnt,3]<-j
      
    }
  }

  result_cbd_MB_PH_bf<-mclapply(1:dim(t_table)[1],function(ind_aly,list1,list1_nm,list2,list2_nm,data,rev){
    rm(res)
    rm(res_rev)
    i<-t_table[ind_aly,2]
    j<-t_table[ind_aly,3]
    print(paste0(i,"_",j))
    Target<-c(list1[i],list2[j])
    Target_nm<-c(list1_nm[i],list2_nm[j])
    SNP_X<-Target[1]
    SNP_Y<-Target[2]
    MR_table_cbd2<-data
    ind_go<-MR_table_cbd2$'X-SNP'==SNP_X
    MRDATA_bf<-MR_table_cbd2[ind_go,]
    MRDATA<-MRDATA_bf[!(is.na(MRDATA_bf[,paste0(SNP_Y,"_beta")])|is.na(MRDATA_bf[,paste0(SNP_X,"_beta")])),]
    n_result<-5+7*5
    if(dim(MRDATA)[1]>3){
      bo<-paste0(SNP_Y,"_beta")
      be<-paste0(SNP_X,"_beta")
      so<-paste0(SNP_Y,"_se")
      se<-paste0(SNP_X,"_se")
      set.seed(1)
      MR_median<-mr_median(mr_input(bx = MRDATA[,be], bxse = MRDATA[,se], by = MRDATA[,bo], byse = MRDATA[,so]),weighting = "weighted", iterations = n_perm,seed=1)
      Egger <- mr_egger(mr_input(bx=MRDATA[,be],bxse=MRDATA[,se],by=MRDATA[,bo],byse=MRDATA[,so]))
      MR_IVW<-mr_ivw(mr_input(bx=MRDATA[,be],bxse=MRDATA[,se],by=MRDATA[,bo],byse=MRDATA[,so]))
      PRESSO <- mr_presso(BetaOutcome = bo, BetaExposure = be, SdOutcome = so, SdExposure = se, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = MRDATA, NbDistribution = n_perm,  SignifThreshold = 0.05)
      
      yyy <<- MRDATA[,paste0(SNP_Y,"_beta")]
      xxx <<- MRDATA[,paste0(SNP_X,"_beta")]
      se_yyy <<- MRDATA[,paste0(SNP_Y,"_se")]
      se_xxx <<- se<-MRDATA[,paste0(SNP_X,"_se")]
      
      dattt <- cbind(yyy,xxx,se_yyy,se_xxx)
      dattt <- as.data.frame(dattt)
      
      (Fit2 <- lm(yyy ~ xxx, weights=1/se_yyy^2,x=T, y=T))
      
    
      
      (mod.sim <- simex(Fit2,B=n_perm,
                            measurement.error = se_xxx,
                            SIMEXvariable="xxx", asymptotic = "FALSE"))
      
      simex_re <- summary(mod.sim)
      # res_rev
      # MR_IVW_rev@Estimate
      # MR_IVW_rev@StdError
      # MR_IVW_rev@Pvalue

      
      
      
      
      q_s<-c(Egger@Heter.Stat[2],MR_IVW@Heter.Stat[2],PRESSO$'MR-PRESSO r'$'Global Test'$Pvalue,Egger@I.sq)
      

      Median_s<-c(MR_median@Estimate,MR_median@StdError, MR_median@Pvalue)
      Egger_Int<-c(Egger@Intercept,Egger@StdError.Int,Egger@Pvalue.Int)
      Egger_Slope<-c(Egger@Estimate,Egger@StdError.Est,Egger@Causal.pval)
      SIMEX_Int<-unlist(simex_re$coefficients$jackknife[1,c(1,2,4)])
      SIMEX_Slope<-unlist(simex_re$coefficients$jackknife[2,c(1,2,4)])
      IVW_s<- c(MR_IVW@Estimate,MR_IVW@StdError, MR_IVW@Pvalue)
      
      
      
      if(is.na(PRESSO$'Main MR results'[2,3])){
        PRESSO_s<-PRESSO$'Main MR results'[1,c(3,4,6)]
        N_presso<-NA
      }else{
        PRESSO_s<-PRESSO$'Main MR results'[2,c(3,4,6)]
        N_presso<-dim(MRDATA)[1]-length(PRESSO$'MR-PRESSO results'$'Distortion Test'$'Outliers Indices')
      }
      
      list_results<-list(Median_s,Egger_Int,Egger_Slope,SIMEX_Int,SIMEX_Slope,IVW_s,PRESSO_s)
      list_results2<-lapply(list_results,add_CIs)
      
      res<-c(ind_aly,dim(MRDATA)[1],N_presso,q_s,do.call("c",list_results2))
      
    }else{
      print("Inst <= 2 skip")
      res<-c(ind_aly,dim(MRDATA)[1],rep(NA,n_result))
    }
    
    if(rev){
      SNP_X<-Target[2]
      SNP_Y<-Target[1]
      
      ind_go<-MR_table_cbd2$'X-SNP'==SNP_X
      MRDATA_bf<-MR_table_cbd2[ind_go,]
     
      
      MRDATA<-MRDATA_bf[!(is.na(MRDATA_bf[,paste0(SNP_Y,"_beta")])|is.na(MRDATA_bf[,paste0(SNP_X,"_beta")])),]
      if(dim(MRDATA)[1]>3){
        bo<-paste0(SNP_Y,"_beta")
        be<-paste0(SNP_X,"_beta")
        so<-paste0(SNP_Y,"_se")
        se<-paste0(SNP_X,"_se")
        set.seed(1)
        MR_median<-mr_median(mr_input(bx = MRDATA[,be], bxse = MRDATA[,se], by = MRDATA[,bo], byse = MRDATA[,so]),weighting = "weighted", iterations = n_perm,seed=1)
        Egger <- mr_egger(mr_input(bx=MRDATA[,be],bxse=MRDATA[,se],by=MRDATA[,bo],byse=MRDATA[,so]))
        
        MR_IVW<-mr_ivw(mr_input(bx=MRDATA[,be],bxse=MRDATA[,se],by=MRDATA[,bo],byse=MRDATA[,so]))
        
        PRESSO <- mr_presso(BetaOutcome = bo, BetaExposure = be, SdOutcome = so, SdExposure = se, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = MRDATA, NbDistribution = n_perm,  SignifThreshold = 0.05)
        
        yyy <<- MRDATA[,paste0(SNP_Y,"_beta")]
        xxx <<- MRDATA[,paste0(SNP_X,"_beta")]
        se_yyy <<- MRDATA[,paste0(SNP_Y,"_se")]
        se_xxx <<- se<-MRDATA[,paste0(SNP_X,"_se")]
        
        dattt <- cbind(yyy,xxx,se_yyy,se_xxx)
        dattt <- as.data.frame(dattt)
        
        (Fit2 <- lm(yyy ~ xxx, weights=1/se_yyy^2,x=T, y=T))
        
        (mod.sim <- simex(Fit2,B=n_perm,
                          measurement.error = se_xxx,
                          SIMEXvariable="xxx", asymptotic = "FALSE"))
        
        simex_re <- summary(mod.sim)
        # res_rev
        # MR_IVW_rev@Estimate
        # MR_IVW_rev@StdError
        # MR_IVW_rev@Pvalue
        
        
        
        
        
        q_s<-c(Egger@Heter.Stat[2],MR_IVW@Heter.Stat[2],PRESSO$'MR-PRESSO r'$'Global Test'$Pvalue,Egger@I.sq)
        
        
        Median_s<-c(MR_median@Estimate,MR_median@StdError, MR_median@Pvalue)
        Egger_Int<-c(Egger@Intercept,Egger@StdError.Int,Egger@Pvalue.Int)
        Egger_Slope<-c(Egger@Estimate,Egger@StdError.Est,Egger@Causal.pval)
        SIMEX_Int<-unlist(simex_re$coefficients$jackknife[1,c(1,2,4)])
        SIMEX_Slope<-unlist(simex_re$coefficients$jackknife[2,c(1,2,4)])
        IVW_s<- c(MR_IVW@Estimate,MR_IVW@StdError, MR_IVW@Pvalue)
        
        
        
        if(is.na(PRESSO$'Main MR results'[2,3])){
          PRESSO_s<-PRESSO$'Main MR results'[1,c(3,4,6)]
          N_presso<-NA
        }else{
          PRESSO_s<-PRESSO$'Main MR results'[2,c(3,4,6)]
          N_presso<-dim(MRDATA)[1]-length(PRESSO$'MR-PRESSO results'$'Distortion Test'$'Outliers Indices')
        }
        
        list_results<-list(Median_s,Egger_Int,Egger_Slope,SIMEX_Int,SIMEX_Slope,IVW_s,PRESSO_s)
        list_results2<-lapply(list_results,add_CIs)
        
        res_rev<-c(ind_aly,dim(MRDATA)[1],N_presso,q_s,do.call("c",list_results2))
        
      }else{
        print("Inst <= 2 skip")
        res_rev<-c(ind_aly,dim(MRDATA)[1],rep(NA,n_result))
      }
    }else{
      res_rev<-c(ind_aly,dim(MRDATA)[1],rep(NA,n_result))
    }
    
    ####Low -> Causal exist
    
    
    nms_togo<-lapply(list("Median","Egger_Int","Egger_Slope","SIMEX_Int","SIMEX_Slope","IVW","PRESSO"),function(strs){
      paste0(strs,"_",c("Est","SE","LCI","UCI","pvalue"))
    })
    
    nm<-c("Index","N_Exposure","N_presso","Q","Qp","RSS","Isq",do.call("c",nms_togo))
    names(res)<-nm
    names(res_rev)<-nm
    
    
    
    result<-cbind(WHAT=c(paste(Target_nm,collapse="#"),paste(Target_nm[2:1],collapse="#")),data.frame(rbind(res,res_rev)))
    
    return(result)
  },mc.cores=n_cores,list1=clist1,list1_nm=clist1_nm,list2=clist2,list2_nm=clist2_nm,data=cdata,rev=crev)
  
  result_cbd_MB_PH<-do.call("rbind",result_cbd_MB_PH_bf)
  # names(result_cbd_MB_PH)<-names(result)
  return(result_cbd_MB_PH)
}

source("~/analysis/Dong/20200608_GWAS/code/toKjkim/mr_presso.R")

# Target<-c("MB2","MG1")
dim(MR_table_cbd_ORI)
dim(MR_table_cbd_DM)


result_cbd_MBMG<-conduct_MR_WM(MB_list,mb_label,MG_list,mg_label,cdata=MR_table_cbd_ORI)
result_cbd_MBPH<-conduct_MR_WM(MB_list,mb_label,PH_list,ph_label,cdata=MR_table_cbd_ORI)
result_cbd_MGPH<-conduct_MR_WM(MG_list,mg_label,PH_list,ph_label,cdata=MR_table_cbd_ORI)

result_cbd_MGDM<-conduct_MR_WM(MG_list,mg_label,DM_list,dm_label,cdata=MR_table_cbd_DM,crev=F)
result_cbd_MBDM<-conduct_MR_WM(MB_list,mb_label,DM_list,dm_label,cdata=MR_table_cbd_DM,crev=F)
result_cbd_PHDM<-conduct_MR_WM(PH_list,ph_label,DM_list,dm_label,cdata=MR_table_cbd_DM,crev=F)

# MR_table_cbd_ORI[MR_table_cbd_ORI$'X-SNP'=="MG1",]
# MR_table_cbd_ORI[MR_table_cbd_ORI$'X-SNP'=="MB1",]




# be<-"MG1_beta"
# se<-"MG1_se"
# bo<-"PH6_beta"
# so<-"PH6_se"
# MR_IVW_rev<-mr_ivw(mr_input(bx=MRDATA[,be],bxse=MRDATA[,se],by=MRDATA[,bo],byse=MRDATA[,so]))


# MRRESULT<-rbind(result_cbd_MBMG,result_cbd_MBPH,result_cbd_MGPH)


MRRESULT_bf<-lapply(list(result_cbd_MBMG,result_cbd_MBPH,result_cbd_MGPH,result_cbd_MGDM,result_cbd_MBDM,result_cbd_PHDM),function(data_in){
  
  data_in2<-lapply(data_in[,-1],function(data){
    as.numeric(data)
  })
  data_in3<-cbind(as.character(data_in[,1]),data.frame(do.call("cbind",data_in2),stringsAsFactors=F))
  names(data_in3)<-names(data_in)
  return(data_in3)
})

MRRESULT<-do.call("rbind",MRRESULT_bf)

# MRRESULT<-rbind(result_cbd_MBMG,result_cbd_MBPH,result_cbd_MGPH)
library(stringr)
r_label<-str_split(MRRESULT$WHAT,"#")
r_label2<-t(sapply(r_label,function(data){
  data
}))
MRRESULT$X<-r_label2[,1]
MRRESULT$Y<-r_label2[,2]
MRRESULT2<-MRRESULT[,c("WHAT","X","Y",setdiff(names(MRRESULT),c("WHAT","X","Y")))]

write.csv(MRRESULT2,"LDproxMRRESULT_ALL1000.csv")
system("scp LDproxMRRESULT_ALL1000.csv ng:~")



MRRESULT2$Exposure
