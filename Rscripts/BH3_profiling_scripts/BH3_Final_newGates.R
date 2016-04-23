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
setwd("~/Dropbox/BH3/Final/Batch3/ForMoom/")
alive=CytoCMed=NULL
filenames<-system("ls B3*.txt", intern=TRUE)
cell_lines=NULL
for (i in 1:length(filenames)){
  cell_lines[i]=strsplit(filenames[i],split = "_")[[1]][2]
}

for(fn in 1:length(filenames)){
  f=read.table(filenames[fn],header = T,sep='\t')[1:220,]
  print(dim(f))
  ##creating a live cell matrix
  tmp_alive=matrix(data = NA,nrow = 10,ncol = 22)##instantiating new matrix that will have 10 rows and 22 columns
  i=1
  for(k in 1:220)
  {
    if(k%%44==1){
      print(k)
      tmp_alive[i:(i+1),]<-rbind(t(f[k:(k+21),2]),t(f[(k+22):(k+43),2]))
      i=i+2
    }
    
  }
  rownames(tmp_alive)<-c(paste(c("Alive_BIM_0.1","Alive_BIM_0.1_","Alive_BAD_10","Alive_BAD_10_", "Alive_NOXA_100", "Alive_NOXA_100_", "Alive_Pos_control","Alive_Pos_control_","Alive_Neg_control", "Alive_Neg_control_"),cell_lines[fn],sep="_"))
  colnames(tmp_alive)<-c("Doxo_0.1","Doxo_0.1_", "Doxo_0.3", "Doxo_0.3_", "Carbo_100", "Carbo_100_", "Carbo_300", "Carbo_300_","AKT_10","AKt_10_", "AKT_30","AKT_30_", "MEK_3","MEK_3_","MEK_10","MEK_10_", "Bafilo_0.001","Bafilo_0.001_","Bafilo_0.01","Bafilo_0.01_","DMSO","DMSO_")
  
  alive=cbind(alive,t(tmp_alive))
  ###converting the file read to a 16X24 matrix; This is to match the plate lay out
  tmp_CytoCMed=matrix(data = NA,nrow = 10,ncol = 22)##instantiating new matrix that will have 16 rows and 24 columns
  i=1
  for(k in 1:220)
  {
    if(k%%44==1){
      print(k)
      tmp_CytoCMed[i:(i+1),]<-rbind(t(f[k:(k+21),3]),t(f[(k+22):(k+43),3]))
      i=i+2
    } 
  }
  rownames(tmp_CytoCMed)<-c(paste(c("CytoCMed_BIM_0.1","CytoCMed_BIM_0.1_", "CytoCMed_BAD_10","CytoCMed_BAD_10_", "CytoCMed_NOXA_100", "CytoCMed_NOXA_100_", "CytoCMed_Pos_control","CytoCMed_Pos_control_","CytoCMed_Neg_control", "CytoCMed_Neg_control_"),cell_lines[fn],sep="_"))
  
  colnames(tmp_CytoCMed)<-c("Doxo_0.1","Doxo_0.1_", "Doxo_0.3", "Doxo_0.3_", "Carbo_100", "Carbo_100_", "Carbo_300", "Carbo_300_","AKT_10","AKt_10_", "AKT_30","AKT_30_", "MEK_3","MEK_3_","MEK_10","MEK_10_", "Bafilo_0.001","Bafilo_0.001_","Bafilo_0.01","Bafilo_0.01_","DMSO","DMSO_")
  
  CytoCMed=cbind(CytoCMed,t(tmp_CytoCMed))
  
}
write.table(alive,"~/Dropbox/BH3/Final/Batch3/BH3_alive_cells.txt",sep='\t',quote=F,col.names = NA)
write.table(CytoCMed,"~/Dropbox/BH3/Final/Batch3/BH3_CytoC.txt",sep='\t',quote=F,col.names = NA)
Alive_CytoC=rbind(alive,colnames(CytoCMed),CytoCMed)
View(Alive_CytoC)
write.table(Alive_CytoC,"~/Dropbox/BH3/Final/Batch3/BH3_AllPlates_Alive_CytoC.txt",sep='\t',quote=F,col.names = NA)

####Peptide normalization###
cyto_c_release<-CytoCMed
cyto_c_release[21,49]<-NA#this replicate was way off than the rest. Therefore,removing it.
#colnames(cyto_c_release)<-Alive_CytoC[25,]
View(cyto_c_release)
dim(cyto_c_release)
peptide_norm<-function(x,pos,neg){
  return (1-((x-pos)/(neg-pos)))
}
###transforming the matrix so that the replicates are side by side
new_cyto_c_release=matrix(data = NA,nrow = 11,ncol = 100)
k=1
for(i in 1:11){
  for(j in 1:100){
    if(j%%4==1){
      new_cyto_c_release[i,j]=cyto_c_release[(2*i-1),k]
      new_cyto_c_release[i,(j+1)]=cyto_c_release[2*i,k]
      new_cyto_c_release[i,(j+2)]=cyto_c_release[(2*i-1),(k+1)]
      new_cyto_c_release[i,(j+3)]=cyto_c_release[2*i,(k+1)]
      k=(k+2)%%50#
    } 
  }
}

View(new_cyto_c_release)
dim(new_cyto_c_release)
mean_cyto=matrix(data = 0,nrow = 11,ncol = 25)
k=1
for(i in 1:11){
  for(j in 1:25){
    #print(mean(new_cyto_c_release[i,(k:(k+3))],na.rm = T))
    mean_cyto[i,j]<-mean(new_cyto_c_release[i,(k:(k+3))],na.rm = T)
    k=(k+4)%%100
  }
}

View(mean_cyto)
dim(mean_cyto)
rownames(mean_cyto)<-c("Doxo_0.1", "Doxo_0.3", "Carbo_100", "Carbo_300", "AKT_10", "AKT_30","MEK_3","MEK_10", "Bafilo_0.001","Bafilo_0.01","DMSO")
colnames_mean_cyto=NULL
j = 1
for (i in 1:ncol(cyto_c_release)){
  print(colnames(cyto_c_release)[i])
  print(i)
  if (i%%2 == 1) {
    tmp=colnames(cyto_c_release)[i]
    colnames_mean_cyto[j] = sub(pattern = "CytoCMed_",replacement =  "",tmp)
    j = j + 1
  }
}
colnames(mean_cyto)<-colnames_mean_cyto
View(mean_cyto)
write.table(mean_cyto,"~/Dropbox/BH3/Final/Batch3/BH3_mean_CytoC_batch3.txt",sep='\t',quote=F,col.names = NA)


##########normalizing the peptide data using the positive and negative peptides#######
normalized_cyto_c_release=matrix(data = NA,ncol = ncol(mean_cyto),nrow = nrow(mean_cyto))
dimnames(normalized_cyto_c_release)=dimnames(mean_cyto)
j=k=1
dmso_neg_peptide=NULL
for(i in 1:(ncol(mean_cyto)/length(filenames))){
  tmp<-mean_cyto[,j:(j+4)]
  pos<-tmp[,4]
  neg<-tmp[,5]
  for(k in 1:nrow(tmp)){
    for(l in 1:5){

      normalized_cyto_c_release[k,(5*(i-1)+l)]=peptide_norm(x =tmp[k,l],pos = pos[k],neg=neg[k] ) 
    }
  }
  dmso_neg_peptide=c(dmso_neg_peptide,rep(peptide_norm(x =tmp[11,5],pos = min(pos,na.rm = T),neg=max(neg,na.rm=T)) ,5))
  j=j+5
}
normalized_cyto_c_release<-rbind(normalized_cyto_c_release,dmso_neg_peptide)


View(normalized_cyto_c_release)
write.table(normalized_cyto_c_release,"~/Dropbox/BH3/Final/Batch3/normalized_BH3_CytoC_batch3.txt",sep='\t',quote=F,col.names = NA)

###alternate normalization where pos=max(tmp) and neg=min(tmp)
alt_normalized_cyto_c_release=matrix(data = NA,ncol = ncol(mean_cyto),nrow = nrow(mean_cyto))
dimnames(alt_normalized_cyto_c_release)=dimnames(mean_cyto)
dmso_neg_peptide=NULL
j=k=1
for(i in 1:(ncol(mean_cyto)/length(filenames))){
  tmp<-mean_cyto[,j:(j+4)]
  pos<-min(tmp,na.rm = T)
  neg<-max(tmp, na.rm = T)
  for(k in 1:nrow(tmp)){
    for(l in 1:5){
      
      alt_normalized_cyto_c_release[k,(5*(i-1)+l)]=peptide_norm(x =tmp[k,l],pos = pos,neg=neg ) 
    }
  }
  dmso_neg_peptide=c(dmso_neg_peptide,rep(peptide_norm(x =tmp[11,5],pos = pos,neg=neg) ,5))
  j=j+5
}
alt_normalized_cyto_c_release=rbind(alt_normalized_cyto_c_release,dmso_neg_peptide)
View(alt_normalized_cyto_c_release)
write.table(alt_normalized_cyto_c_release,"~/Dropbox/BH3/Final/Batch3/alternative_normalized_BH3_CytoC_batch3.txt",sep='\t',quote=F,col.names = NA)
###
