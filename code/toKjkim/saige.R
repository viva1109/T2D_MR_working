library(SAIGE)

args = commandArgs(trailingOnly=TRUE)
chr <- args[1]

setwd("~/NTM_revision/1.DiscoverySet/5.Imputation/")

#if (chr == 1) {
#	covar <- read.table("../covar.txt",head=T)
#
#	fam <- read.csv(paste0("4.QC/chr",chr,"/chr",chr,".fam"),sep=" ", comment.char = "",head=F)
#	colnames(fam) <-c ("FID","IID","PID","MID","SEX","PHENO")
#
#	covarTable <- merge(fam,covar,by="IID",sort=F)
#	covarTable <- covarTable[,c(6,2,1,8:ncol(covarTable))]
#	colnames(covarTable)[2] <- "FID"
#	covarTable$PHENO <- covarTable$PHENO - 1
#	sID <- paste0(covarTable$FID,"_",covarTable$IID)
#
#	write.table(covarTable,"7.SAIGE/pheno.txt",quote=F,row.names=F,sep="\t")
#	write.table(covarTable$IID,"7.SAIGE/sample_id.txt",quote=F,col.names=F,row.names=F,sep="\t")
#}

res1 <- fitNULLGLMM(plinkFile = paste0("4.QC/chr",chr,"/chr",chr), phenoFile = "7.SAIGE/pheno.txt", phenoCol = "PHENO",traitType="binary",covarColList=c("sex","bmi","age","age2"),sampleIDColinphenoFile="IID",nThreads=5,outputPrefix=paste0("7.SAIGE/chr",chr,"/chr",chr))

res2 <- SPAGMMATtest(dosageFile = paste0("6.PrediXcan/chr",chr,"/chr",chr,".txt"), dosageFileNrowSkip=0, dosageFileNcolSkip=6,minMAF=0.0001,sampleFile="7.SAIGE/sample_id.txt",GMMATmodelFile = paste0("7.SAIGE/chr",chr,"/chr",chr,".rda"),varianceRatioFile = paste0("7.SAIGE/chr",chr,"/chr",chr,".varianceRatio.txt"), SAIGEOutputFile = paste0("7.SAIGE/chr",chr,"/chr",chr,".saige.res"),numLinesOutput=100,IsOutputAFinCaseCtrl=TRUE)


