library(MRInstruments)
data(gwas_catalog)
head(gwas_catalog)
data(metab_qtls)
head(metab_qtls)

MET2<-metab_qtls[metab_qtls$phenotype=="AcAce",]
COME<-cbind(MET2$chromosome,MET2$position,data.frame(MET2$SNP,stringsAsFactors=F),-MET2$beta,MET2$se,MET2$pval)
names(COME)<-c("CHr","Position","RSno","beta","se","pval")
beta_outcome<-c(-0.007582,0.004415)
se_outcome<-c(0.002382,0.002103)
pval_outcome<-c(0.001457,0.03573)
chr_outcome<-c(8 ,11)
phypo_outcome<-c(9181395,116648917 )
COME2<-cbind(COME,beta_outcome,se_outcome,pval_outcome)

 
TargetX<-"MB2"
TargetY<-"PH1"
SNP_Y<-TargetY
SNP_X<-TargetX
MRMR_toplot<-MR_table_cbd_ORI[MR_table_cbd_ORI$'X-SNP'==TargetX,c(paste0(c("X-SNP"),c("_beta","_se","_pval")),paste0(c(TargetX),c("_beta","_se","_pval")),paste0(c(TargetY),c("_beta","_se","_pval")))]

# MRMR_toplot
names(MRMR_toplot)<-gsub("X-SNP","XSNP",names(MRMR_toplot))
MRMR_toplot1<-MRMR_toplot[complete.cases(MRMR_toplot),]
MRMR_toplot2<-MRMR_toplot[complete.cases(MRMR_toplot),c(paste0(c("XSNP"),c("_beta")),paste0(c(TargetY),c("_beta")))]
names(MRMR_toplot1)<-gsub(TargetY,"Y",names(MRMR_toplot1))
names(MRMR_toplot2)<-gsub(TargetY,"Y",names(MRMR_toplot2))


library(ggplot2)
# pd 
labelX<-"Acetoacetate"
labelY<-"HbA1c"


library(simex)

library(ieugwasr)


yyy <<- MRMR_toplot[,paste0(SNP_Y,"_beta")]
xxx <<- MRMR_toplot[,paste0(SNP_X,"_beta")]
se_yyy <<- MRMR_toplot[,paste0(SNP_Y,"_se")]
se_xxx <<- se<-MRMR_toplot[,paste0(SNP_X,"_se")]


#yyy <<- mer_scz2$Effect
#xxx <<- log(mer_scz2$OR)
#se_yyy <<- mer_scz2$StdErr
#se_xxx <<- mer_scz2$SE


dattt <- cbind(yyy,xxx,se_yyy,se_xxx)
dattt <- as.data.frame(dattt)

(Fit2 <- lm(yyy ~ xxx, weights=1/se_yyy^2,x=T, y=T))

set.seed(1)
(mod.sim <- simex(Fit2,B=10000,
                  measurement.error = se_xxx,
                  SIMEXvariable="xxx", asymptotic = "FALSE"))

simex_re <- summary(mod.sim); simex_re
simex_inter_lower <- simex_re$coefficients$jackknife[1,1]-1.96*simex_re$coefficients$jackknife[1,2]; simex_inter_lower
simex_inter_upper <- simex_re$coefficients$jackknife[1,1]+1.96*simex_re$coefficients$jackknife[1,2]; simex_inter_upper
simex_slope_lower <- simex_re$coefficients$jackknife[2,1]-1.96*simex_re$coefficients$jackknife[2,2]; simex_slope_lower
simex_slope_upper <- simex_re$coefficients$jackknife[2,1]+1.96*simex_re$coefficients$jackknife[2,2]; simex_slope_upper

MRDATA<-MRMR_toplot
bo<-paste0(SNP_Y,"_beta")
be<-paste0(SNP_X,"_beta")
so<-paste0(SNP_Y,"_se")
se<-paste0(SNP_X,"_se")
MR_median<-mr_median(mr_input(bx = MRDATA[,be], bxse = MRDATA[,se], by = MRDATA[,bo], byse = MRDATA[,so]),weighting = "weighted", iterations = 10000,seed=1)

res<-c(1,dim(MRDATA)[1],c(MR_median@Estimate, MR_median@StdError, MR_median@Pvalue))



res_rev <- mr_presso(BetaOutcome = bo, BetaExposure = be, SdOutcome = so, SdExposure = se, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = MRDATA, NbDistribution = 1000,  SignifThreshold = 0.05)
res_rev
Egger_rev <- mr_egger(mr_input(bx=MRDATA[,be],bxse=MRDATA[,se],by=MRDATA[,bo],byse=MRDATA[,so]))

####low -> pleio exist.
Egger_rev@Pleio.pval
####Low -> Causal exist
Egger_rev@Causal.pval
####No. SNPs
Egger_rev@SNPs
##low -> Weak
Egger_rev@I.sq

Egger_rev@StdError.Est



egger_slope<-Egger_rev@Estimate
egger_int<-Egger_rev@Intercept

Wmedian_slope<-MR_median@Estimate

betaJ<-MRMR_toplot1$Y_beta/MRMR_toplot1$XSNP_beta
weights<-(MRMR_toplot1$XSNP_beta^2)/(MRMR_toplot1$Y_se^2)
betaCausal<-(weights%*%betaJ)/sum(weights)
simex_inter<-simex_re$coefficients$jackknife[1,1]
simex_slop<-simex_re$coefficients$jackknife[2,1]
max(MRMR_toplot1$XSNP_beta)


# COME2$beta,COME2$beta_outcome

COME3<-COME2[,-c(1:3)]
COME4<-cbind(COME3[,1:3],COME3,Data="UKB")

names(COME4)<-c(names(MRMR_toplot1),"Data")
MRMR_toplot2<-rbind(cbind(MRMR_toplot1,Data="KARE"),COME4)
DFintSlop<-data.frame(intercept=c(0,simex_inter,0,egger_int),slope=c(betaCausal,simex_slop,Wmedian_slope,egger_slope),group=c("IVW","MR-Egger(SIMMEX)","Weighted median","MR-Egger"))
DFintSlop$group<-factor(DFintSlop$group)

  
  ggplot(MRMR_toplot1,aes(x=XSNP_beta,y= Y_beta)) +
    geom_errorbar(aes(ymin=Y_beta-Y_se, ymax=Y_beta+Y_se), width=15,col="black",alpha=0.5)+
    geom_errorbarh(aes(xmin=XSNP_beta-XSNP_se,xmax=XSNP_beta+XSNP_se), height=0.001,col="black",alpha=0.5)+
    geom_point(aes())+
    stat_summary(fun.data=mean_cl_normal) + geom_point(color='black') +
    geom_abline(data=DFintSlop,aes(intercept = intercept, slope = slope,colour = group,linetype=group))+
    # xlim(c(min(min(MRMR_toplot1$XSNP_beta)*1.1,0),max(max(MRMR_toplot1$XSNP_beta)*1.1,0))) +
    theme_bw()+xlab(paste0("Genetic assocition with exposure: ",labelX))+ylab(paste0("Genetic assocition with outcome: ",labelY))+
    theme(legend.position="bottom")
  
  ggsave(paste0(labelX,"_",labelY,"LDproxyMB2.tiff"),device="tiff",width=6,height=6,dpi=300)
  
  system(paste0("scp ",paste0(labelX,"_",labelY,"LDproxyMB2.tiff ng:~")))
  


betaJCOME<-COME4$Y_beta/COME4$XSNP_beta
weightsCOME<-(COME4$XSNP_beta^2)/(COME4$Y_se^2)
betaCausalCOME<-(weightsCOME%*%betaJCOME)/sum(weightsCOME)




ggplot(COME4,aes(x=XSNP_beta,y= Y_beta)) +
  geom_errorbar(aes(ymin=Y_beta-Y_se, ymax=Y_beta+Y_se), width=0.01,col="black",alpha=0.5)+
  geom_errorbarh(aes(xmin=XSNP_beta-XSNP_se,xmax=XSNP_beta+XSNP_se), height=0.001,col="black",alpha=0.5)+
  geom_point(aes())+
  stat_summary(fun.data=mean_cl_normal) + geom_point(color='black') +
  geom_abline(aes(intercept = 0, slope = betaCausalCOME,colour = "red"))+
  # xlim(c(min(min(MRMR_toplot1$XSNP_beta)*1.1,0),max(max(MRMR_toplot1$XSNP_beta)*1.1,0))) +
  theme_bw()+xlab(paste0("Genetic assocition with exposure: ",labelX))+ylab(paste0("Genetic assocition with outcome: ",labelY))

ggsave(paste0(labelX,"_",labelY,"LDproxyMB2COME.tiff"),device="tiff",width=6,height=6,dpi=300)

system(paste0("scp ",paste0(labelX,"_",labelY,"LDproxyMB2COME.tiff ng:~")))

COME4

MR_median<-mr_median(mr_input(bx = COME4[,be], bxse = COME4[,se], by = COME4[,"Y_beta"], byse = COME4[,"Y_se"]),weighting = "weighted", iterations = 10000,seed=1)

MR_IVW_rev<-mr_ivw(mr_input(bx=COME4[,be],bxse=COME4[,se],by=COME4[,"Y_beta"],byse=COME4[,"Y_se"]))
MR_IVW_rev@Estimate
MR_IVW_rev@StdError
MR_IVW_rev@Pvalue