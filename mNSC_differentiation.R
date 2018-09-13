########################################################################
# PHN Gene Expression Analysis 
########################################################################
library(ggstance)
library(sleuth)
library(reshape2)
library(ggplot2)
library(ggpubr)

library(data.table)
library(limma)

SO<-sleuth_load("~/GitHub/mNSC-Analysis/SleuthApp/so.rds")

#LR test
Diff_sleuth_table <- sleuth_results(SO, 'age:full', 'lrt', show_all = FALSE)
Age_sleuth_table <- sleuth_results(SO, 'diff:full', 'lrt', show_all = FALSE)
AgeDiff_sleuth_table <- sleuth_results(SO, 'null:full', 'lrt', show_all = FALSE)

Diff_sleuth_table<-merge(Diff_sleuth_table,annotables::grcm38,by.x = "target_id",by.y="ensgene")
Age_sleuth_table<-merge(Age_sleuth_table,  annotables::grcm38,by.x = "target_id",by.y="ensgene")

Diff_sleuth_significant <- dplyr::filter(Diff_sleuth_table, qval <= 0.01)
Age_sleuth_significant <- dplyr::filter(Age_sleuth_table, qval <= 0.01)

write.csv(Diff_sleuth_significant,file = "mNSC_Diff_Sig-q01.csv")
write.csv(Age_sleuth_significant,file = "mNSC_AgeDiff_Sig-q01.csv")

#Wald test
Wald_Diff_sleuth_table <- sleuth_results(SO,test = "diffD",test_type ='wt', show_all = F)
Wald_Age_sleuth_table <- sleuth_results(SO,test = "age3mo",test_type ='wt', show_all = F)

WDiff_sleuth_table<-merge(Wald_Diff_sleuth_table,annotables::grcm38,by.x = "target_id",by.y="ensgene")
WAge_sleuth_table<-merge(Wald_Age_sleuth_table  ,annotables::grcm38,by.x = "target_id",by.y="ensgene")

WDiff_sleuth_significant<-dplyr::filter(WDiff_sleuth_table, qval <= 0.01)
WAge_sleuth_significant<-dplyr::filter(WAge_sleuth_table, qval <= 0.01)

write.csv(WDiff_sleuth_significant,file = "~/GitHub/mNSC-Analysis/DataTables/mNSC_Diff_Sig-Wald-q01.csv")
write.csv(WAge_sleuth_significant,file = "~/GitHub/mNSC-Analysis/DataTables/mNSC_Age_Sig-Wald-q01.csv")

ggplot(SO$tests$wt$full$diffD)+geom_point(aes(b,qval,col=qval<0.01))+scale_y_log10()
ggplot(SO$tests$wt$full$age3mo)+geom_point(aes(b,qval,col=qval<0.01))+scale_y_log10()

ggplot(SO$tests$lrt$`age:full`)+geom_point(aes(pval,qval,col=qval<0.01))+scale_y_log10()
ggplot(SO$tests$lrt$`diff:full`)+geom_point(aes(pval,qval,col=qval<0.01))+scale_y_log10()




########################################################################
# PHN analysis 
########################################################################

Chaps<-read.delim("~/GitHub/mNSC-Analysis/DataTables/Chaperones.txt")
Chaps$Names<-strsplit2(Chaps$Gene.names,split = " ")
Chaps<-reshape2::melt(Chaps)
Proteasome<-read.delim("~/GitHub/mNSC-Analysis/DataTables/proteasome.txt")
PSMlist<-data.frame(stringsAsFactors = F,symbol=unique(c(strsplit2(Proteasome$Gene.names,split = " "))),subset="",set="Proteasome")
ChapList<-unique(data.frame(stringsAsFactors = F,symbol=Chaps$Names,subset=Chaps$Class,set="Chaperone"))

UBI<-read.delim("~/GitHub/mNSC-Analysis/DataTables/Ubl.txt")
UBIlist<-data.frame(stringsAsFactors = F,symbol=unique(c(strsplit2(UBI$Gene.names,split = " "))),subset="",set="Ubiquitin")

allNames<-rbind.data.frame(ChapList, PSMlist, UBIlist)

modelAnnot<-function(model){
  model$ensgene<-model$target_id
  genesDF<-merge.data.frame(model,data.frame(stringsAsFactors = F,annotables::grcm38_tx2gene),by="ensgene")
  annotDF<-merge(genesDF,annotables::grcm38,by="ensgene")
  annot.DF<-merge(annotDF,allNames,by="symbol",all.x = T)
  annot.DF$qval[annot.DF$qval<1E-20]<-1E-20
  return(annot.DF)
  print("values limited to 1E-20")
}

model=SO$tests$lrt$`age:full`
LRTannot<-modelAnnot(model)

model=SO$tests$wt$full$diffD
WTannot<-modelAnnot(model)

ggsave(useDingbats=F,file="~/GitHub/mNSC-Analysis/Figures/all.pdf",height=4,width=5,
       ggplot(WTannot[!is.na(WTannot$set),])+
         geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=set))+
         geom_vline(xintercept=0)+
         #geom_hline(yintercept=1e-25,lty=2)+
         geom_text(data=WTannot[!is.na(WTannot$set)&WTannot$qval<0.01,],cex=2,aes(b,-log10(qval),label=symbol))
)

ggsave(useDingbats=F,file="~/GitHub/mNSC-Analysis/Figures/PSM.pdf",height=4,width=5,
       ggplot(WTannot[WTannot$set=="Proteasome"&(!is.na(WTannot$set)),])+
         geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=set))+
         geom_vline(xintercept=0)+
         #geom_hline(yintercept=1e-25,lty=2)+
         geom_text(data=WTannot[WTannot$set=="Proteasome"&WTannot$qval<0.01,],cex=2,aes(b,-log10(qval),label=symbol))
         
)

ggsave(useDingbats=F,file="~/GitHub/mNSC-Analysis/Figures/Chaps_grid.pdf",height=1.7,width=7,
       ggplot(WTannot[(!is.na(WTannot$set))&!(is.na(WTannot$subset)),])+
         geom_pointrangeh(size=.2,aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=subset))+
         geom_vline(xintercept=0)+
         #geom_hline(yintercept=1e-25,lty=2)+
         geom_text(nudge_x = .3,cex=2,data=unique(WTannot[(WTannot$qval<0.01)&(WTannot$set=="Chaperone")&(!is.na(WTannot$set))&(WTannot$subset!=""),-19]),aes(b,-log10(qval),label=symbol))+
         facet_grid(~subset,scales = 'free_y')+scale_color_brewer(palette = "Dark2")+coord_cartesian(xlim = c(-5.5,5.5))
)

ggsave(useDingbats=F,file="~/GitHub/mNSC-Analysis/Figures/Chaps.pdf",height=4,width=5,
       ggplot(WTannot[WTannot$set=="Chaperone"&(!is.na(WTannot$set))&(WTannot$subset!=""),])+
         geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=subset))+
         geom_vline(xintercept=0)+
         geom_text(nudge_x = .3,cex=2,data=unique(WTannot[WTannot$set=="Chaperone"&(!is.na(WTannot$set))&(WTannot$subset!=""),-19]),aes(b,-log10(qval),label=symbol))+
         scale_color_brewer(palette = "Dark2")+xlim(-4,NA)
)


ggsave(useDingbats=F,file="~/GitHub/mNSC-Analysis/Figures/Ubiq.pdf",height=4,width=5,
       ggplot(WTannot[WTannot$set=="Ubiquitin"&(!is.na(WTannot$set)),])+
         geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=b+se_b,alpha=qval<0.01,col=set))+
         geom_vline(xintercept=0)+
         #geom_hline(yintercept=1e-25)+
         geom_text(data=unique(WTannot[!is.na(WTannot$set)&WTannot$set=="Ubiquitin"&WTannot$qval<0.01,-19]),aes(b,-log10(qval),cex=1,label=symbol))
)



################


genesDF<-data.frame(SO$obs_norm_filt)
annotDF<-merge(genesDF,annotables::grcm38,by.x = "target_id",by.y="ensgene")
annot.DF<-merge(annotDF,allNames,by="symbol",all.x = T)

annot.DF$age=strsplit2(annot.DF$sample,split = "-")[,1]
annot.DF$diff=strsplit2(annot.DF$sample,split = "-")[,3]
annot.DF$sample=strsplit2(annot.DF$sample,split = "-")[,2]

mergedDF<-merge.data.frame(WTannot,annot.DF)

ggplot(annot.DF[!is.na(annot.DF$set),])+geom_boxplot(aes(set,tpm,col=diff))+scale_y_log10()

ggplot(annot.DF[!is.na(annot.DF$subset),])+geom_boxplot(aes(subset,tpm,col=diff))+scale_y_log10()

#ggsave(useDingbats=F,width= 6,height=6,filename = "SmoothedRPB.pdf",
ggsave(useDingbats=F,"~/GitHub/mNSC-Analysis/Figures/RPB_ChangePlot.pdf",width=7, height=5,
       ggplot(mergedDF[!is.na(mergedDF$subset),])+
         geom_line(aes(diff,scaled_reads_per_base,group=paste(age,sample,target_id)),alpha=0.3)+
         geom_smooth(data = mergedDF[!is.na(mergedDF$subset)&mergedDF$qval<0.01,],aes(diff,scaled_reads_per_base,group=paste(target_id),col=subset))+
         geom_boxplot(cex=.5,mapping = aes(x = factor(diff),y=scaled_reads_per_base,fill=subset,group=diff),width=0.3)+
         facet_grid(set~subset,scales = "free_y")+
         #geom_smooth(mapping = aes(diff,scaled_reads_per_base,color=subset,group=subset),method = "lm")+
         #coord_cartesian(ylim = c(1,100000))+
         scale_y_log10()+
         theme_bw()+scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")
)

ggsave(useDingbats=F,"~/GitHub/mNSC-Analysis/Figures/TPM_ChangePlot.pdf",width=7, height=5,
       ggplot(mergedDF[!is.na(mergedDF$subset),])+
         geom_line(aes(diff,tpm,group=paste(age,sample,target_id)),alpha=0.3)+
         geom_smooth(data = mergedDF[!is.na(mergedDF$subset)&mergedDF$qval<0.01,],aes(diff,tpm,group=paste(target_id),col=subset))+
         #geom_violin(mapping = aes(x = factor(diff),y=scaled_reads_per_base,color=subset,group=diff))+
         geom_boxplot(cex=.5,mapping = aes(x = factor(diff),y=tpm,fill=subset,group=diff),width=0.3)+
         facet_grid(set~subset,scales = "free_y")+
         #geom_smooth(mapping = aes(diff,scaled_reads_per_base,color=subset,group=subset),method = "lm")+
         #coord_cartesian(ylim = c(1,100000))+
         scale_y_log10()+
         theme_bw()+scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")
)

means<-dcast(annot.DF,symbol~diff,value.var=c("tpm"),fun.aggregate=mean)
sds<-dcast(annot.DF,symbol~diff,value.var=c("tpm"),fun.aggregate=sd)
annot<-dcast(annot.DF,symbol~diff,value.var=c("subset"),fun.aggregate=function(X){as.character(X[1])})
qvals<-dcast(mergedDF,symbol~diff,value.var=c("qval"),fun.aggregate=mean)

ggsave(useDingbats=F,"~/GitHub/mNSC-Analysis/Figures/Diff_meanTPM_Scatter.pdf",width=6,height=5,
       ggplot()+geom_abline(slope=1)+theme_bw()+
         geom_point(data=means,aes(Var.2,D),col="black",alpha=0.2)+
         geom_abline(slope=1)+
         geom_point(pch=21,data=means[annot$Var.2!=""&!is.na(as.factor(annot$Var.2)),],cex=2.5,aes(Var.2,D,fill=as.factor(annot$Var.2[annot$Var.2!=""&!is.na(as.factor(annot$Var.2))])))+scale_x_log10()+scale_y_log10()+
         guides(fill=guide_legend(title= "Chaperone Class"))+
         ylab("mean Transcripts/Million, Differentiated mNSCs")+
         xlab("mean Transcripts/Million, Undifferentiated mNSCs")
)

means<-dcast(annot.DF,symbol~diff,value.var=c("scaled_reads_per_base"),fun.aggregate=mean)
sds<-dcast(annot.DF,symbol~diff,value.var=c("scaled_reads_per_base"),fun.aggregate=sd)
annot<-dcast(annot.DF,symbol~diff,value.var=c("subset"),fun.aggregate=function(X){as.character(X[1])})
qvals<-dcast(mergedDF,symbol~diff,value.var=c("qval"),fun.aggregate=mean)

ggsave(useDingbats=F,"~/GitHub/mNSC-Analysis/Figures/Diff_meanRPB_Scatter.pdf",width=6,height=5,
       ggplot()+geom_abline(slope=1)+theme_bw()+
         geom_point(data=means,aes(Var.2,D),col="black",alpha=0.2)+
         geom_abline(slope=1)+
         geom_point(pch=21,data=means[annot$Var.2!=""&!is.na(as.factor(annot$Var.2)),],cex=2.5,aes(Var.2,D,fill=as.factor(annot$Var.2[annot$Var.2!=""&!is.na(as.factor(annot$Var.2))])))+scale_x_log10()+scale_y_log10()+
         guides(fill=guide_legend(title= "Chaperone Class"))+
         ylab("mean scaled reads/base, Differentiated mNSCs")+
         xlab("mean scaled reads/base, Undifferentiated mNSCs")
)

ggsave(useDingbats=F,filename = "~/GitHub/mNSC-Analysis/Figures/SmoothedCount.pdf",width=6,height=5,
       ggplot(annot.DF[!is.na(annot.DF$subset),])+
         geom_line(aes(diff,scaled_reads_per_base,group=target_id,col=subset),alpha=0.2)+
         facet_wrap(~subset)+
         stat_summary(data=annot.DF[annot.DF$diff=="D",],aes(x=diff,group=target_id,y=scaled_reads_per_base,label=target_id),cex =2,geom = "text",fun.y = function(X){mean(X,na.rm=T)})+
         geom_smooth(mapping = aes(diff,scaled_reads_per_base,color=subset,group=subset),method = "lm", lwd=2)+scale_y_log10()+theme_pubr()+scale_color_brewer(palette="Spectral")
)


############################
MarkerGenes<-read.delim("~/GitHub/mNSC-Analysis/DataTables/MarkerGenes.txt",header = T,stringsAsFactors = F)

AnnotGenes<-merge(MarkerGenes,annot.DF,by.x = "ensgene", by.y="target_id")

AnnotStats<-merge(AnnotGenes,Wald_Diff_sleuth_table,by.x = "ensgene",by.y = "target_id")
AnnotStats$qval[AnnotStats$qval<1E-50]<-1E-50

ggsave(useDingbats=F,"~/GitHub/mNSC-Analysis/Figures/OtherAnnotation_ChangePlot.pdf",width = 7, height = 5,
ggplot(data=AnnotStats)+theme_bw()+
  #geom_smooth(method = "loess",aes(color=celltype,group=paste(celltype,age),diff,scaled_reads_per_base),alpha=0.2)+
  geom_smooth(lwd=1,aes(color=celltype,group=paste(ensgene),diff,scaled_reads_per_base),alpha=0.2)+
  stat_summary(data=AnnotStats[AnnotStats$diff=="D",],aes(x=diff,group=gene,y=scaled_reads_per_base,label=gene),cex =2,geom = "text",fun.y = function(X){mean(X,na.rm=T)})+
  guides(fill=guide_legend(title= "Class"))+
  scale_y_log10()+facet_wrap(~celltype)+
  ylab("mean scaled reads/base")+
  xlab("Differentiation State")
)

ggsave(useDingbats=F,"~/GitHub/mNSC-Analysis/Figures/OtherAnnotation_VolcanoPlot.pdf",width = 7, height = 5,
ggplot(AnnotStats)+
  geom_pointrangeh(aes(b,-log10(qval),xmin=(b-se_b),xmax=(b+se_b),color=celltype))+
  geom_text(data=unique(AnnotStats[AnnotStats$qval<0.01,-c(6:19)]),cex =2,aes(b,-log10(qval),label=symbol))
)

