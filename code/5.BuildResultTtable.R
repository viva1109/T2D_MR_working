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
dolist_info<-readRDS(file.path(C3_PATH_OF,FOL_snp,"dolist_info.rds"))
ttable<-dolist_info[[1]]
table_chunks2<-dolist_info[[2]]
ttable_BIG<-dolist_info[[3]]
head(ttable)
head(table_chunks2)

# table_chunks2[1:80,]

Result_fol_list1<-file.path(C3_PATH_OF,FOL_snp,"TTT2_GWA_midoutput")
Result_fol_list2<-file.path(Result_fol_list1,list.files(Result_fol_list1))
Result_file_lists<-lapply(Result_fol_list2,function(t_path){
  file.path(t_path,list.files(t_path))
})

Result_file_path<-do.call("c",Result_file_lists)

Chunk_done<-as.numeric(str_extract(Result_file_path,"(?<=D)\\d+(?=_1000.rds)"))
length(Chunk_done)
Sys.time()
Chunk_undone<-setdiff(1:dim(table_chunks2)[1],Chunk_done)

length(Chunk_undone)
# saveRDS(Chunk_undone,file.path(C3_PATH_OF,FOL_snp,"chunk_undone_NEW.rds"))

CuttedChunks<-table(cut(Chunk_undone,
    breaks=seq(1,dim(table_chunks2)[1],by=2000),
    include.lowest=F))

library(ggplot2)
dfdf<-data.frame(Range=names(CuttedChunks),CuttedChunks)
ggplot(data=dfdf,aes(x=factor(Range),y=Freq))+geom_bar(stat="identity")

hist(Chunk_undone)
sort(Chunk_undone)
tail(Chunk_undone)

Result_files<-lapply(Result_file_path,readRDS)
Result_cbd<-do.call("rbind",Result_files)
head(Result_cbd)

Result_cbd2<-cbind(ttable_BIG[match(Result_cbd$Ind_ttable,ttable_BIG$Index),],Result_cbd)

dim(Result_cbd2)

# Result_cbd2<-cbind(Result_cbd,ttable[match(Result_cbd$Ind_ttable,ttable$Index),c("CHR","SNP_index","Y")])


setdiff(Result_cbd$Ind_ttable,ttable_BIG$Index)



# SNP_raw<-lapply(SNP_list,fread)
# SNP_name_ext<-str_extract(names(SNP_raw[[1]])[-(1:6)],".+(?=_[ATCG]$)")

midresult_fol<-file.path(C3_PATH_OF,"SNP","TTT2_GWA_midoutput")
plot_fol<-file.path(C3_PATH_OF,"SNP","TTT2_GWA_plot")
dir.create(midresult_fol,recursive = T)
dir.create(plot_fol,recursive = T)
lapply(midresult_fol,dir.create)
lapply(plot_fol,dir.create)



var_candi<-unique(ttable$Y)
prefix<-"newOG"
for ( j in 1:length(var_candi)){
  var_target<-var_candi[j]
  
  
  CBD_cbd<-Result_cbd2[Result_cbd2$Y==var_target,]
  names(CBD_cbd)[c(2,3,4)]<-c("CHR","SNP","BP")
  CBD_cbd<-data.frame(CBD_cbd,check.names = F)
  sort(unique(CBD_cbd$CHR))
  # CBD_cbd[is.na(CBD_cbd$CHR),]
  
  label_chrs<-paste(c(min(chr_list),max(chr_list)),collapse = "-")
  # make_qq(CBD_cbd,outdir=file.path(plot_fol,paste0(prefix,"QQ_LR_",var_target,label_chrs,".png")),p="PR(>F)",width=3000,height=3000,res=600,units='px')
  make_qq(CBD_cbd,outdir=file.path(plot_fol,paste0(prefix,"QQ_Wald_",var_target,label_chrs,".png")),p="Pr(>|t|)",width=3000,height=3000,res=600,units='px')
  make_manhattan(CBD_cbd,file.path(plot_fol,paste0(prefix,"MAN_Wald_",var_target,label_chrs,".png")), 
                 chr_coln = 'CHR', bp_coln='BP', snp_coln='SNP', p_coln='Pr(>|t|)', width=4937, height=2598, res=660, units='px')
  # make_manhattan(CBD_cbd,file.path(plot_fol,paste0(prefix,"MAN_LR_",var_target,label_chrs,".png")), 
  #                chr_coln = 'CHR', bp_coln='BP', snp_coln='SNP', p_coln='PR(>F)', width=4937, height=2598, res=660, units='px')
  
  system(paste0("scp ",file.path(plot_fol,paste0(prefix,"QQ_*",var_target,label_chrs,".png"))," ng:~"))
  system(paste0("scp ",file.path(plot_fol,paste0(prefix,"MAN_*",var_target,label_chrs,".png"))," ng:~"))
  print(var_target)
}
