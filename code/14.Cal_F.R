
library(data.table)
library(stringr)
library(parallel)


n_core<-22

output_fol<-"/home2/kjkim/analysis/Dong/20200608_GWAS"
C1_fol<-"/KARE_QC"
C2_fol<-"2.DataReady"
C3_fol<-"3.Result"
FOL_snp<-"SNP"

DATA_fol<-c("PH","MB","MG")
DATA_ind<-2
log_cri<-F
C2_PATH_OF<-file.path(output_fol,C2_fol,DATA_fol[DATA_ind])
C3_PATH_OF<-file.path(output_fol,C3_fol,DATA_fol[DATA_ind])
# chunk_undone<-readRDS(file.path(C3_PATH_OF,FOL_snp,"chunk_undone_NEW.rds"))

# do_list<-chunk_undone[601:length(chunk_undone)]
do_list<-2901:3800
# chr_picked<-"chr21"

output_files<-paste0(output_fol,C1_fol,c("/M_QC/","/M_QC/","/M_QC/snp_qc","/I_QC/1.sex.out","/I_QC/2.mind.out","/I_QC/3.hetero.out","/I_QC/4.Clean_IND_SNP","/3.PC/1.Prune","/3.PC/2.Pruned_data","/3.PC/3.Clean_IND_SNP_PC","/4.IBS_PCoutlier/1.IBS","/4.IBS_PCoutlier/2.PC_outlier","/3.PC/PRUNED_chr_split/2.Pruned_data_chr"))
names(output_files)<-c("M_QC_fol","I_QC_fol","M_QC","I_QC1","I_QC2","I_QC3","Clean_IND_SNP","PRUNNING","PRUNNED","Clean_IND_SNP_PC","IBS","PC_outlier","PRUNNED_chr")
path_plink<-"/bin2/plink"
path_wisard<-"/bin2/wisard"

SNPsInfo_bf<-read.csv("~/SNPsInfo.csv")

SNPsInfo<-SNPsInfo_bf[str_sub(SNPsInfo_bf$'X.SNP',1,-2)==DATA_fol[DATA_ind],]

chr_list2<-1:22
SNP_fol_list<-file.path(C2_PATH_OF,paste0(FOL_snp,"/chr",chr_list2))

OUTCOME_fol<-file.path(C2_PATH_OF,"OUTCOME")
names(output_files)
list.files("/home2/kjkim/analysis/Dong/20200608_GWAS/KARE_QC/3.PC/")
PCs<-fread(paste0(output_files["Clean_IND_SNP_PC"],".eigenvec"))



OC<-read.csv(file.path(OUTCOME_fol,"Pheno.csv"),stringsAsFactors = F)
head(OC)
OC$AGE_p6<-OC$AGE[OC$TIME==6]
OC$AGE_sq<-OC$AGE_p6^2

OC$AGE_p6_sc<-scale(OC$AGE_p6)
OC$AGE_sq_sc<-scale(OC$AGE_sq)



OC$SEX<-as.factor(OC$SEX)
# OC$Y<-
if(DATA_ind==1){
  var_togo<-names(OC)[c(6:12,17,18)]  
}else if(DATA_ind==2){
  var_togo<-c("X3.Hydroxybutyrate","Acetoacetate")  
}else if(DATA_ind==3){
  var_togo<-c("C34","R34")  
}


head(OC)

data_ready<-OC[,c("IID","SEX","AGE_p6_sc","AGE_sq_sc","TIME")]
i<-1
df_logT<-matrix(nrow=length(var_togo),ncol=2)
for ( i in 1:length(var_togo)){
  Y<-OC[,var_togo[i]]
  l_Y<-log(Y)
  shT<-shapiro.test(Y[which(OC$TIME==min(OC$TIME))])
  shT2<-shapiro.test(l_Y[which(OC$TIME==min(OC$TIME))])
  logTr<-shT$p.value < shT2$p.value
  df_logT[i,1]<-var_togo[i]
  df_logT[i,2]<-logTr
}

paste0(var_togo[as.logical(df_logT[,2])],"_log")
head(OC)

ind_logtr<-which(as.logical(df_logT[,2]))
tmp_logOC<-OC[,var_togo[ind_logtr]]
names(tmp_logOC)<-paste0(var_togo[which(as.logical(df_logT[,2]))],"_log")
for (i in 1:dim(tmp_logOC)[2]){
  tmp_logOC[,i]<-log(OC[,var_togo[i]])
}

OC2<-cbind(OC,tmp_logOC)
# OC2<-OC
# head(OC2)

if(DATA_ind==1){
  var_togo_bf<-c(var_togo[-ind_logtr],names(tmp_logOC))
  var_togo2<-c(setdiff(var_togo_bf,"WHR"),"WHR")
}else if(DATA_ind==2&log_cri){
  print(DATA_ind)
  var_togo_bf<-c(var_togo[-ind_logtr],names(tmp_logOC))
  var_togo2<-var_togo_bf
}else{
  print(DATA_ind)
  var_togo2<-var_togo
}





# PATH_target<-file.path(C2_PATH_OF,FOL_snp,chr_picked)
# chr_no_goto<-chr_list2
# PATH_target_list_chr<-file.path(C2_PATH_OF,FOL_snp,paste0("chr",chr_no_goto))


# listlist<-readRDS(file.path(C3_PATH_OF,FOL_snp,paste0("dolist_info_",DATA_fol[DATA_ind],".rds")))
# ttable<-listlist[[1]]
# table_chunks2<-listlist[[2]]
# ttable_BIG<-listlist[[3]]
# 
# 
# chosen_chunks<-table_chunks2[do_list,]
# 
# 
# 
# 
# chosen_works<-min(chosen_chunks[,2]):max(chosen_chunks[,3])
# chosen_ttable<-ttable[chosen_works,]





chosen_chr_list<-sort(unique(SNPsInfo$CHR))



PATH_target_list<-file.path(C2_PATH_OF,FOL_snp,paste0("chr",chosen_chr_list))

SNP_path<-file.path(PATH_target_list,"Imputed_cleanIND_QC1_mafhwe_recodedA.raw")
SNPsdf_list<-lapply(SNP_path,function(path){
  data.frame(fread(path))
})
names(SNPsdf_list)<-chosen_chr_list
names(PATH_target_list)<-chosen_chr_list
# SNPsdf_list[[as.character(21)]][1:6,1:6]
start_SNP<-7
names(PCs)[-c(1:2)]<-paste0("PC",1:10)

library(lmerTest)


# input<-1:(dim(SNPsdf)[2]-6)

# make_chenk_input<-function(input,chunk1=1000,n_core=24){
#   len_test<-length(input)
#   chunk2<-ceiling(len_test/n_core)
#   len_test%/%chunk2
#   chunk<-min(chunk1,chunk2)
#   J<-ceiling(len_test/chunk)
#   list_input<-vector("list",J)
#   print("J")
#   print(J)
#   
#   for (j in 1:J){
#     list_input[[j]]<-input[(1+chunk*(j-1)):min(chunk*j,len_test)]
#   }
#   return(list(list_input,chunk))
# }
# n_core<-96
# chuckoutput<-make_chenk_input(input,n_core = n_core)
# list_input<-chuckoutput[[1]]
chuck_selected<-1000

# resultGWA_list<-vector("list",length(var_togo2))
library(parallel)
library(nlme)
# SNPsdf_list[[as.character(21)]][1:6,1:6]
# do_list
OC_PCs<-PCs[match(data_ready$IID,PCs$V2),-c(1:2)]
sum(OC2$IID!=data_ready$IID)
# listlist<-list(ttable,table_chunks2,ttable_BIG)
# names(listlist)<-c("ttable","table_chunks2","SNPINFO_ttable")
# saveRDS(listlist,file.path(C3_PATH_OF,FOL_snp,"dolist_info.rds"))


SNPsInfo_splitedM<-split(SNPsInfo,SNPsInfo$CHR)

data_ready2<-cbind(data_ready,OC_PCs,R=OC2)
names(data_ready2)<-gsub("R.RID","RID",names(data_ready2))

anno_files<-"/data/sharedcode/kjkim/Ardo/rsID_snp138"
annofiles_path<-file.path(anno_files,list.files(anno_files))[order(as.numeric(str_extract(list.files(anno_files),"(?<=chr).+(?=\\.txt$)")))]
tmp_SNP<-data.frame(matrix(nrow=dim(data_ready2)[1],ncol=0))
for ( w in 1:length(SNPsInfo_splitedM)){
  CHR_ch<-as.character(SNPsInfo_splitedM[[w]]$CHR[1])
  SNPsdf_list2<-SNPsdf_list
  SNPSNPdfdf<-SNPsdf_list2[[CHR_ch]]
  SNPSNPdfdf_ready<-SNPSNPdfdf[match(data_ready2$IID,SNPSNPdfdf$IID),]
  # for ( q in 1:dim(SNPSNPdfdf_ready)[2]){
      BIMFILE<-fread(file.path(PATH_target_list,"Imputed_cleanIND_QC1_mafhwe_recodedA.bim")[w])
      BIMFILE$V2<-gsub("-","\\.",BIMFILE$V2)
      # BIMFILE$V2
      anno<-fread(annofiles_path[as.numeric(CHR_ch)])
      position_left<-unlist(BIMFILE[BIMFILE$V2!="rs","V4"])
      NM_right<-unlist(anno[match(position_left,anno$bp),"rsID"])
      names(SNPSNPdfdf_ready)[-c(1:6)]<-NM_right
    # }    
  # }
  
  tmp_SNP<-cbind(tmp_SNP,SNPSNPdfdf_ready[which(names(SNPSNPdfdf_ready)%in%as.character(SNPsInfo_splitedM[[w]]$RS))])
  print(w)
}



data_ready3<-cbind(data_ready2,tmp_SNP[1:dim(SNPsInfo)[1]])
RIDT_right<-paste0(data_ready3$RID,"_",data_ready3$TIME)

data_ready3_2<-read.csv("~/data_ready_SNP_MG.csv",stringsAsFactors = F)[,-1]
RIDT_left<-paste0(data_ready3_2$RID,"_",data_ready3_2$TIME)
length(setdiff(RIDT_right,RIDT_left))/3


data_ready4<-data_ready3[complete.cases(tmp_SNP),]
snps_include<-setdiff(names(tmp_SNP),"rs2259835")
fmla<-formula(paste0("R.X3.Hydroxybutyrate~",paste(snps_include,collapse="+"),"+factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"))
fmla2<-formula(paste0("R.X3.Hydroxybutyrate~factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"))
m2<-lme(fmla,random=~1|IID,data=data_ready4,method="ML",control=lmeControl(opt="optim"))
m1<-lme(fmla2,random=~1|IID,data=data_ready4,method="ML",control=lmeControl(opt="optim"))
rersme<-summary(m2)
coef_res<-rersme$tTable
ano<-anova(m2,adjustSigma = T)
anova(m2,m1)

data_ready3_PH<-read.csv("~/Data_target2_PH.csv",stringsAsFactors = F)[-1]

RIDT_right2<-paste0(data_ready3_PH$RID,"_",data_ready3_PH$TIME)


data_cbd_MBMGPH<-cbind(data_ready3_2,data_ready3[match(RIDT_left,RIDT_right),],P=data_ready3_PH[match(RIDT_left,RIDT_right2),])

data_cbd_MBMGPH$P.INS60_log<-log(data_cbd_MBMGPH$P.INS60)
data_cbd_MBMGPH$P.HBA1C_log<-log(data_cbd_MBMGPH$P.HBA1C)
corM<-cor(data_cbd_MBMGPH[,c("P.INS60_log","P.HBA1C_log","C34","R.X3.Hydroxybutyrate","R.Acetoacetate")],use="pairwise.complete.obs")
MM<-as.matrix(data_cbd_MBMGPH[,c("P.INS60_log","P.HBA1C_log","C34","R.X3.Hydroxybutyrate","R.Acetoacetate")])

png("~/corplot2.png")
plot(data_cbd_MBMGPH[,c("P.INS60_log","P.HBA1C_log","C34","R.X3.Hydroxybutyrate","R.Acetoacetate")])
dev.off()
system("scp ~/corplot2.png ng:~")

.libPaths("~/tmp")
library(mediation)
# install.packages("stargazer")
library(stargazer)
data_cbd_MBMGPH<-data_cbd_MBMGPH[,-c(52,53)]
MG_RS_label<-names(data_cbd_MBMGPH)[20:25]
MB_RS_label<-names(data_cbd_MBMGPH)[52:59]
# i<-2

label_Y<-"R.X3.Hydroxybutyrate"
label_X<-"C34"
targetLabel<-setdiff(MB_RS_label,"rs2259835")
library(mediation)

res_cb<-data.frame(matrix(nrow=length(targetLabel),ncol=9))

res_cb$Y<-label_Y
res_cb$X<-label_X
for ( i in 1:length(targetLabel)){
  X <- data_cbd_MBMGPH[,targetLabel[i]]
  M <- data_cbd_MBMGPH[,label_X]
  #M <- combined5$M_Acetoacetate
  Y <- data_cbd_MBMGPH[,label_Y]
  Meddata <- data.frame(X, M, Y)
  Meddata2<-Meddata[complete.cases(Meddata),]
  library(mediation)
  fitM <- lm(M ~ X,     data=Meddata2) #IV on M; Hours since dawn predicting coffee consumption
  fitY <- lm(Y ~ X + M, data=Meddata2) #IV and M on DV; Hours since dawn and coffee predicting wakefulness
  # gvlma(fitM) #data is positively skewed; could log transform (see Chap. 10 on assumptions)
  fitMed <- mediate(fitM, fitY, treat="X", mediator="M")
  summary(fitMed)

  
  #1. Total Effect
  fit <- lm(Y ~ X, data=Meddata)
  # summary(fit)
  
  Y_X_p<-summary(fit)$coef[2,4]
  Y_X_b<-summary(fit)$coef[2,1]
  
  #2. Path A (X on M)
  fita <- lm(M ~ X, data=Meddata)
  # summary(fita)
  
  M_X_p<-summary(fita)$coef[2,4]
  M_X_b<-summary(fita)$coef[2,1]
  
  #3. Path B (M on Y, controlling for X)
  fitb <- lm(Y ~ M + X, data=Meddata)
  # summary(fitb)
  
  McX_p<-summary(fitb)$coef[2,4]
  McX_b<-summary(fitb)$coef[2,1]
  
  XcM_p<-summary(fitb)$coef[3,4]
  XcM_b<-summary(fitb)$coef[3,1]
  
  #4. Reversed Path C (Y on X, controlling for M)
  fitc <- lm(X ~ Y + M, data=Meddata)
  
  Rev_p<-summary(fitc)$coef[2,4]
  
  
  res<-c(targetLabel[i],Y_X_p,Y_X_b,M_X_p,McX_b,McX_p,XcM_b,XcM_p,Rev_p)
  #install.packages("stargazer")
  #library(stargazer)
  #Summary Table
  # MG_RS_label[i]
  # stargazer(fit, fita, fitb, fitc, type = "text", title = "Baron and Kenny Method")
  res_cb[i,1:9]<-res
}


names(res_cb)[1:9]<-c("RS","Y_X_p","Y_X_b","M_X_p","McX_b",'McX_p','XcM_b',"XcM_p",'Rev_p')

write.csv(res_cb,paste0("~/",label_X,"_",label_Y,"_MED.csv"))
system(paste0("scp ",paste0("~/",label_X,"_",label_Y,"_MED.csv ng:~")))

# write.csv(data_ready3,"~/data_ready_SNP_MG.csv")

# for ( i in 1:length(var_togo2)){
rm(Y)
resultGWA<-mclapply(do_list,function(do_k){
  tc<-table_chunks2[do_k,]
  ctb<-ttable[tc[2]:tc[3],]
  chrchr_list<-unique(ctb$CHR)
  SNPsdf_list2<-SNPsdf_list[chrchr_list]
  chunkGroup<-paste0("G",tc["Group"])
  
  if(DATA_ind==3){
    df_SNP_res<-data.frame(matrix(nrow=dim(ctb)[1],ncol=7))
  }else{
    df_SNP_res<-data.frame(matrix(nrow=dim(ctb)[1],ncol=6))
  }
  
  # kaka1<-data_ready2$IID[data_ready2$TIME==6]
  # kaka2<-SNPSNPdfdf[match(data_ready2$IID[data_ready2$TIME==6],SNPSNPdfdf$IID),2]
  # sum(kaka1!=kaka2)
  # lala1<-rep(SNPSNPdfdf[match(data_ready2$IID[data_ready2$TIME==6],SNPSNPdfdf$IID),SNP_ind+6],3)
  # lala2<-SNPSNPdfdf[match(data_ready2$IID,SNPSNPdfdf$IID),SNP_ind+6]
  # sum(lala1!=lala2,na.rm=T)
  # sum(which(is.na(lala1))!=which(is.na(lala2)))
  # length(data_ready$IID)
  # length(SNPSNPdfdf[match(data_ready2$IID,SNPSNPdfdf$IID),2])
  # sum(data_ready$IID!=SNPSNPdfdf[match(data_ready2$IID,SNPSNPdfdf$IID),2])
  
  for ( j in 1:dim(ctb)[1]){
    if(log_cri){
      V_target<-paste0(ctb[j,"Y"],"_log")
    }
    
    SNP_ind<-ctb[j,"SNP_index"]
    CHR_ch<-ctb[j,"CHR"]
    head(data_ready)
    
    data_ready2<-cbind(data_ready,OC_PCs,Y=OC2[,V_target])
    SNPSNPdfdf<-SNPsdf_list2[[CHR_ch]]
    data_ready3<-cbind(data_ready2,X_SNP=SNPSNPdfdf[match(data_ready2$IID,SNPSNPdfdf$IID),SNP_ind+6])
    data_ready4<-data_ready3[!(is.na(data_ready3$Y)|is.na(data_ready3$X_SNP)),]
    # relme_H1<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4,REML=F)
    
    if(DATA_ind==3){
      m2<-lme(Y~X_SNP+factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,random=~1|IID,data=data_ready4,method="REML",control=lmeControl(opt="optim"))
      rersme<-summary(m2)
      coef_res<-rersme$tTable
      ano<-anova(m2,adjustSigma = T)
    }else{
      relme_H1_re<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4,REML=T)
      sm_re<-summary(relme_H1_re)
      coef_res<-coef(sm_re)
    }

    # coef_res
    # anova(m2,adjustSigma=F)
    # 
        # relme_H1_re<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4,REML=)
    # AIC(relme_H1_re)
    
    # relme_H1_re0<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+AGE_sq_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4)
    # AIC(relme_H1_re0)
    # relme_H1_re1<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+AGE_sq_sc+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4)
    # AIC(relme_H1_re1)
    # relme_H1_re2<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+TIME+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4)
    # AIC(relme_H1_re2)
    # relme_H1_re3<-lmer(Y~X_SNP+factor(SEX)+AGE_p6_sc+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+(1|(IID)), data=data_ready4)
    # AIC(relme_H1_re3)
    

    # re_LR<-drop1(relme_H1,test="Chisq")
    if(DATA_ind==3){
      df_SNP_res[j,]<-c(ctb[j,"Index"],unlist(coef_res["X_SNP",]),ano["X_SNP",4])  
    }else{
      df_SNP_res[j,]<-c(ctb[j,"Index"],unlist(coef_res["X_SNP",]))
    }
    
    
    # if(j%%10==0){
    #   print(paste0(CHR_ch,"_",V_target,"_",SNP_ind))
    # }
  }
  print(paste0("Chunk",do_k,"_done"))
  # nms<-c("SNP_ind","Estimate","Std. Error","df","t value","Pr(>|t|)","Sum Sq","Mean Sq","NumDF","DenDF","F value","PR(>F)")
  if(DATA_ind==3){
    nms<-c("Ind_ttable","Estimate","Std. Error","df","t value","Pr(>|t|)","LRP")
  }else{
    nms<-c("Ind_ttable","Estimate","Std. Error","df","t value","Pr(>|t|)")
  }
  
  
  names(df_SNP_res)<-nms
  if(log_cri){
    dir.create(file.path(C3_PATH_OF,FOL_snp,"TTT2_log_GWA_midoutput",chunkGroup),recursive = T)
    saveRDS(df_SNP_res,file.path(C3_PATH_OF,FOL_snp,"TTT2_log_GWA_midoutput",chunkGroup,paste0("D",do_k,"_",chuck_selected,".rds")))    
  }else{
    dir.create(file.path(C3_PATH_OF,FOL_snp,"TTT2_GWA_midoutput",chunkGroup),recursive = T)
    saveRDS(df_SNP_res,file.path(C3_PATH_OF,FOL_snp,"TTT2_GWA_midoutput",chunkGroup,paste0("D",do_k,"_",chuck_selected,".rds")))    
  }

  return(df_SNP_res)
},mc.cores = n_core)

# resultGWA_list[[i]]<-resultGWA
# print(paste0(var_togo2[i],"_done"))

# }
?lmer
# readRDS(df_SNP_res,file.path(PATH_target_toy,"GWA_midoutput",paste0(var_togo2[i],"_",1:20,".rds")))

