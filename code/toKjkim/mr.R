library(MendelianRandomization)
library(ggplot2)
source("~/tools/mr_presso.R")
library(parallel)

setwd("/data4/ktpark/NTM_revision/1.DiscoverySet/5.Imputation")

## Load Data

logit <- read.table("5.Logistic/NTMxTWIN_PC10.add.pvalTable.txt",head=T)
eqtl <- read.table("16.MR/Lung.v7.signif_variant_STK17A_pairs.txt",head=T,stringsAsFactors=F)
chr <- sapply(eqtl$variant_id,function(x)strsplit(x,"_")[[1]][1])
bp <- sapply(eqtl$variant_id,function(x)strsplit(x,"_")[[1]][2])
snp <- paste0(chr,":",bp)

eqtl <- data.frame(eqtl,SNP=snp)

dat <- merge(eqtl,logit,by="SNP",sort=F)
dim(dat) 

## MR-PRESSO adjusting pairwise correlation of SNPs by boostraping
# Estimate beta 
system("plink --bfile 4.QC/chr7/chr7 --extract 16.MR/mr.snp.list --indep-pairwise 50 5 0.2 --out 16.MR/pruning")
pruned_snp <- unlist(read.table("16.MR/pruning.prune.in",head=F,stringsAsFactors=F))

tmp <- dat[dat$SNP %in% pruned_snp,]; tmp

alt <- sapply(tmp$variant_id,function(x)strsplit(x,"_")[[1]][4])
tmp_loc <- which(alt!=tmp$A1)
tmp[tmp_loc,]$slope <- -tmp[tmp_loc,]$slope
tmp <- data.frame(tmp,logOR=log(tmp$OR))

res <- mr_presso(BetaOutcome = "logOR", BetaExposure = "slope", SdOutcome = "SE", SdExposure = "slope_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = tmp, NbDistribution = 1000,  SignifThreshold = 0.05)
beta <- res$"Main MR results"$"Causal Estimate"

# Generate empirical distribution of beta by boostraping

GenDis <- function(iter) {

	if(iter %% 100 == 0) print(iter)
	set.seed(100*iter)
	foo <- dat[sample(c(1:nrow(dat)),5,replace=T),]; foo
	alt <- sapply(foo$variant_id,function(x)strsplit(x,"_")[[1]][4])

	foo_loc <- which(alt!=foo$A1)
	foo[foo_loc,]$slope <- -foo[foo_loc,]$slope
	foo <- data.frame(foo,logOR=log(foo$OR))

	foo_res <- mr_presso(BetaOutcome = "logOR", BetaExposure = "slope", SdOutcome = "SE", SdExposure = "slope_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = foo, NbDistribution = 1000,  SignifThreshold = 0.05)
	foo_beta <- foo_res$"Main MR results"$"Causal Estimate"

	return(foo_beta)

}

beta_dist <- do.call(rbind,mclapply(1:100000,GenDis,mc.cores=40))
beta_res <- apply(beta_dist,1,function(x)ifelse(is.na(x[2]),x[1],x[2]))

z <- (beta[1] - mean(beta_res))/sd(beta_res)
pval <- pnorm(-abs(z))*2; pval

sd(beta_res) 
quantile(beta_res,probs=c(0.025,0.975))

## other MR methods
# mr_egger
res <- mr_egger(mr_input(bx=tmp$slope,bxse=tmp$slope_se,by=log(tmp$OR),byse=tmp$SE))
res <- mr_egger(mr_input(bx=tmp$slope,bxse=tmp$slope_se,by=log(tmp$OR),byse=tmp$SE,corr=cor_mat))

# others
res_all <- mr_allmethods(mr_input(bx=tmp$slope,bxse=tmp$slope_se,by=log(tmp$OR),byse=tmp$SE,corr=cor_mat))
res_median <- mr_median(mr_input(bx=tmp$slope,bxse=tmp$slope_se,by=log(tmp$OR),byse=tmp$SE,corr=cor_mat),weighting="penalized")
res_ivw <- mr_ivw(mr_input(bx=tmp$slope,bxse=tmp$slope_se,by=log(tmp$OR),byse=tmp$SE,corr=cor_mat))
res_conmix <- mr_conmix(mr_input(bx=tmp$slope,bxse=tmp$slope_se,by=log(tmp$OR),byse=tmp$SE,corr=cor_mat))


