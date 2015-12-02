#install.packages("xlsx")
#library(xlsx)


if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

##reading in plate files
setwd("~/Dropbox/BH3/Final/")
filenames<-system("find . -name '*alt*.txt'", intern=TRUE)
filenames

Batch1=read.table(filenames[1],header = T,sep='\t', row.names = T)
Batch2=read.table(filenames[2],header = T,sep='\t', row.names = T)
Batch3=read.table(filenames[3],header = T,sep='\t', row.names = T)

AllBatches=merge(Batch1,Batch2)
AllBatches=merge(AllBatches,Batch3)
colnames(AllBatches)

colnames(AllBatches)[2:26]=paste(colnames(AllBatches)[2:26],"B1",sep='_')
colnames(AllBatches)[27:49]=paste(colnames(AllBatches)[27:49],"B2",sep='_')
colnames(AllBatches)[50:76]=paste(colnames(AllBatches)[50:76],"B3",sep='_')

Subtypes<-as.matrix(c("Subtypes", rep("luminal/HER2 positive",5), rep("basal",5), rep("basal/HER2 positive",5),rep("basal",5), rep("luminal",5), rep("claudin-low",5),rep("luminal/HER2 positive",5),rep("basal",5),rep("basal/HER2 positive",5),rep("luminal",5),rep("basal",5),rep("basal/HER2 positive",5),rep("claudin-low",5),rep("luminal/HER2 positive",5),rep("luminal/HER2 positive",5)))

AllBatches=t(AllBatches)
AllBatches=cbind(Subtypes,AllBatches)
colnames(AllBatches)=AllBatches[1,]
AllBatches=AllBatches[2:nrow(AllBatches),]

Controls=AllBatches[,8:9]
AllBatches=AllBatches[,-8]
AllBatches=AllBatches[,-8]

AllBatches=cbind(AllBatches,Controls)
colnames(AllBatches)[13]="DMSO_No_Peptide"

write.table(AllBatches,"~/Dropbox/BH3/Final/Alt_Normalized_All_Batches.txt",sep='\t',quote=F,col.names = NA)


pdf("~/Dropbox/BH3/Batch2/Batch2_0918.pdf")
heatmap.2(as.matrix(bad[,1:4]),col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Response with \nBAD 10 uM")
heatmap.2(as.matrix(bim),col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Response with \nBIM 0.1 uM")
heatmap.2(as.matrix(noxa),col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Response with\nNOXA 100 uM")
tmp_drug=(cor(t(bad),method="spearman"))
tmp_drug[c(1,4:12),c(1,4:12)]
heatmap.2(tmp_drug[c(1,4:12),c(1,4:12)],col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Drug correlations \nin BAD 10 uM")
tmp_drug=(cor(t(bim),method="spearman"))
tmp_drug[c(1:6,8:11),c(1:6,8:11)]
heatmap.2(tmp_drug[c(1,4:12),c(1,4:12)],col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Drug correlations \nin BIM 0.1 uM")
tmp_drug=(cor(t(noxa),method="spearman"))
tmp_drug[c(1:6,8:11),c(1:6,8:11)]
heatmap.2(tmp_drug[c(1,4:12),c(1,4:12)],col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Drug correlations \nin NOXA 100 uM")




cell_lines_subtypes<-read.table("~/Dropbox/BH3/cellLine_subtypes.txt", row.names = 1, header = 1,sep='\t')
subtype_bad<-merge_drop(cell_lines_subtypes,t(bad))
for(i in (ncol(cell_lines_subtypes)+1):ncol(subtype_bad)){
  boxplot2(subtype_bad[,i]~subtype_bad[,5],main=paste("Response across subtypes:\n BAD 10 uM with",colnames(subtype_bad)[i],sep=" "))
}
subtype_bim<-merge_drop(cell_lines_subtypes,t(bim))
for(i in (ncol(cell_lines_subtypes)+1):ncol(subtype_bim)){
  boxplot2(subtype_bim[,i]~subtype_bim[,5],main=paste("Response across subtypes:\n BIM 0.1 uM with",colnames(subtype_bim)[i],sep=" "))
}
subtype_noxa<-merge_drop(cell_lines_subtypes,t(noxa))
for(i in (ncol(cell_lines_subtypes)+1):ncol(subtype_noxa)){
  boxplot2(subtype_noxa[,i]~subtype_noxa[,5],main=paste("Response across subtypes:\n NOXA 100 uM with",colnames(subtype_noxa)[i],sep=" "))
}
dev.off()





for(i in 1:ncol(mean_normalized)){
  tmp_col=mean_normalized[,i]
  #colnames(tmp_col)=colnames(mean_normalized)[i]
  #print(names(tmp_col))
  #names(tmp_col)=colnames(mean_normalized)[i]
  if(strsplit(colnames(mean_normalized)[i],split = "_")[[1]][1]=="BAD"){
    bad=cbind(bad,tmp_col)
    colnames(bad)[ncol(bad)]=strsplit(colnames(mean_normalized)[i],split = "_")[[1]][3]
    print(colnames(bad))
    #colnames(bad)[-1]=colnames(mean_normalized)[i]
    #print(colnames(bad)[-1])
  }
  else if(strsplit(colnames(mean_normalized)[i],split = "_")[[1]][1]=="BIM"){
    bim=cbind(bim,tmp_col)
    colnames(bim)[ncol(bim)]=strsplit(colnames(mean_normalized)[i],split = "_")[[1]][3]
  }
  else if(strsplit(colnames(mean_normalized)[i],split = "_")[[1]][1]=="NOXA"){
    noxa=cbind(noxa,tmp_col)
    colnames(noxa)[ncol(noxa)]=strsplit(colnames(mean_normalized)[i],split = "_")[[1]][3]
  }
  else{
    controls=cbind(controls,tmp_col)
    colnames(controls)[ncol(controls)]=strsplit(colnames(mean_normalized)[i],split = "_")[[1]][3]
  }
  tmp_col=NULL
}
pdf("~/Dropbox/BH3/Batch1/Batch1.pdf")
heatmap.2(as.matrix(bad),col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Response with BAD 10 uM")
heatmap.2(as.matrix(bim),col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Response with BIM 0.1 uM")
heatmap.2(as.matrix(noxa),col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Response with NOXA 100 uM")
tmp_drug=(cor(t(bad),method="spearman"))
tmp_drug[c(1:6,8:11),c(1:6,8:11)]
heatmap.2(tmp_drug[c(1:6,8:11),c(1:6,8:11)],col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Drug correlations in BAD 10 uM")
tmp_drug=(cor(t(bim),method="spearman"))
tmp_drug[c(1:6,8:11),c(1:6,8:11)]
heatmap.2(tmp_drug[c(1:6,8:11),c(1:6,8:11)],col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Drug correlations in BIM 0.1 uM")
tmp_drug=(cor(t(noxa),method="spearman"))
tmp_drug[c(1:6,8:11),c(1:6,8:11)]
heatmap.2(tmp_drug[c(1:6,8:11),c(1:6,8:11)],col = bluered,density.info = "none",trace="none",margins = c(10,8),main="Drug correlations \nin NOXA 100 uM")





cell_lines_subtypes<-read.table("~/Dropbox/BH3/cellLine_subtypes.txt", row.names = 1, header = 1,sep='\t')
subtype_bad<-merge_drop(cell_lines_subtypes,t(bad))
for(i in (ncol(cell_lines_subtypes)+1):ncol(subtype_bad)){
  boxplot2(subtype_bad[,i]~subtype_bad[,5],main=paste("Response across subtypes:\n BAD 10 uM with",colnames(subtype_bad)[i],sep=" "))
}
subtype_bim<-merge_drop(cell_lines_subtypes,t(bim))
for(i in (ncol(cell_lines_subtypes)+1):ncol(subtype_bim)){
  boxplot2(subtype_bim[,i]~subtype_bim[,5],main=paste("Response across subtypes:\n BIM 0.1 uM with",colnames(subtype_bim)[i],sep=" "))
}
subtype_noxa<-merge_drop(cell_lines_subtypes,t(noxa))
for(i in (ncol(cell_lines_subtypes)+1):ncol(subtype_noxa)){
  boxplot2(subtype_noxa[,i]~subtype_noxa[,5],main=paste("Response across subtypes:\n NOXA 100 uM with",colnames(subtype_noxa)[i],sep=" "))
}
dev.off()