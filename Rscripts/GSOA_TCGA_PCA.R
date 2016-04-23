---
  title: " "
output: html_document
---
install.packages("ROCR")
library(ROCR)
library(GSOA)
  
getwd()
# name the file paths
dataFilePath = "~/Downloads/PANCAN24_BRCA_1119_TPMlog2.txt"
classFilePath = "~/Dropbox/TCGA_PCAClasses_PC1.txt"
gmtFilePath = "/Users/shelley/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/GMT_Files/c2.cp.v4.0.symbols.gmt"
  
#Read in the files

# BRCA RNA-seq data
data = as.matrix(read.table(dataFilePath, sep="\t", stringsAsFactors=F, header=T, row.names=1, check.names=F, quote="\""))
dim(data)

# Read class data from file
classData = as.matrix(read.table(classFilePath, sep="\t", stringsAsFactors=F, header=T, row.names=1, check.names=F, quote="\""))

dim(classData)
# Read the gene-set data from the GMT file
geneSetDatabase <- getGmt(gmtFilePath)

# Run GSOA
results_PC1 = GSOA(data, classData, geneSetDatabase, classificationAlgorithm="svm", numCrossValidationFolds=1, numRandomIterations=1, numCores=4)
View(results_PC1)

View(results_PC1)

