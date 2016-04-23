
# "December 16, 2014"
#Download and load packages

source("http://bioconductor.org/biocLite.R")
biocLite("GSVA")
biocLite("GSVAdata")
biocLite("limma")
biocLite("RColorBrewer")
biocLite("edgeR")
biocLite("gage")
library(GSVA)
library(gage)
library(limma)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(RColorBrewer)

# name file paths
BRCAdataFilePath = "~/Downloads/PANCAN24_BRCA_1119_TPMlog2.txt"
classFilePath = "~/Dropbox/TCGA_PCAClasses_PC1.txt"
gmtFilePath = "/Users/shelley/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/GMT_Files/c2.cp.v4.0.symbols.gmt"

# Read in the BRCA RNA-seq data, Class File, and GmtFiles 
TCGA_Breast_RNAseq = as.matrix(read.table(BRCAdataFilePath, sep="\t", stringsAsFactors=F, header=T, row.names=1, check.names=F, quote="\""))
TCGA_Breast_RNAseq
classDataPC1 = as.matrix(read.table(classFilePath, sep="\t", stringsAsFactors=F, header=T, row.names=1, check.names=F, quote="\""))
View(classDataPC1)
dim(classDataPC1)

geneSetDatabase <- readList(gmtFilePath)
length(geneSetDatabase)


# run GSVA on PC1
PC1_gsva <- gsva(TCGA_Breast_RNAseq, geneSetDatabase , min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=TRUE, rnaseq=TRUE)$es.obs
View(PC1_gsva)
dim(PC1_gsva)

IDs= row.names(classDataPC1)
View(IDs)

PC1_gsva_filtered=merge(IDs,PC1_gsva, by=1)
head(PC1_gsva_filtered)

PC1_gsva_t=t(PC1_gsva)
head(PC1_gsva_t)
dim(PC1_gsva_t)
PC_gsva_classfile=merge_drop(PC1_gsva_t,classDataPC1)
PC_gsva_classfile=PC_gsva_classfile[,-1321]
PC1_gsva_new=t(PC_gsva_classfile)
View(PC1_gsva_new)
View()

head(colnames(PC1_gsva_new))
head(colnames(PC1_gsva))
head(classDataPC1)


colnames(PC1_gsva_t)

dim(PC_gsva_classfile)
colnames(PC_gsva_classfile)

# make the model matrix with the class file
Breast_design_PC1<- model.matrix(~ factor(classDataPC1[,1]))

dim(Breast_design_PC1)
colnames(Breast_design_PC1) <- c("ALL", "HighPCvsLowPC") # 0 is high and 1 is low 





Breast_fitPC1 <- lmFit(PC1_gsva, Breast_design_PC1)

head(Breast_fitC2)
head(Breast_fitC6)
Breast_fitC2 <- eBayes(Breast_fitC2)
Breast_fitC6 <- eBayes(Breast_fitC6)
View(Breast_fitC2)
View(Breast_fitC6)

Breast_allGeneSetsC2 <- topTable(Breast_fitC2, coef="WTvsM", number=Inf)
View(Breast_allGeneSetsC2)

Breast_allGeneSetsC6 <- topTable(Breast_fitC6, coef="WTvsM", number=Inf)
View(Breast_allGeneSetsC6)

write.table(Breast_allGeneSetsC2, file = "RAS_C2_GSVA.txt", sep = "\t")
write.table(Breast_allGeneSetsC6, file = "RAS_C6_GSVA.txt", sep = "\t")
```



GSOA comparisons RAS (RNA-seq)
```{r}
gmtFileName_cp_C2="~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/GMT_Files/c2.cp.v4.0.symbols.gmt"
gmtFileC2_cp=readList(gmtFileName_cp_C2)
length(gmtFileC2_cp) #1320

gmtFileName_C6_Bild="~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/GMT_Files/c6.all.v4.0.symbols_Bild.gmt"
gmtFile_C6_Bild=readList(gmtFileName_C6_Bild)
length(gmtFile_C6_Bild) #194

gmtFile_CGP="~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/GMT_Files/c2.cgp.v4.0.symbols.gmt"
gmtFile_CGP=readList(gmtFile_CGP)
length(gmtFile_CGP) #3402

gmtFile_C2_Bild="~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/GMT_Files/c2.cp.v4.0.symbols_Bild.gmt"
gmtFile_C2_Bild=readList(gmtFile_C2_Bild)
length(gmtFile_C2_Bild) #1325


# Read in the TCGA RNA-seq data for breast only 
TCGA_Lung_RNAseq <- read.delim("~/Documents/ThesisWork/GSOA/RNA_Seq_Files/PANCAN12.IlluminaHiSeq_RNASeqV2.geneExp.tumor_whitelist_Lung_nodup", row.names=1, check.names=FALSE, stringsAsFactors=F)
TCGA_Lung_RNAseq=as.matrix(TCGA_Lung_RNAseq)

head(TCGA_Lung_RNAseq)
dim(TCGA_Lung_RNAseq) # 20499 x 169 RNA-seq samples
row.names(TCGA_Lung_RNAseq)
colnames(TCGA_Lung_RNAseq)

head(gmtFileC2_cp)
head(gmtFile_C6_Bild)
head(gmtFile_C2_Bild)
head(gmtFile_CGP)

#removed values with constant expression
RAS_C2_gsva <- gsva(TCGA_Lung_RNAseq, gmtFileC2_cp, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=TRUE, rnaseq=TRUE)$es.obs
View(RAS_C2_gsva)

RAS_C6_gsva <- gsva(TCGA_Lung_RNAseq, gmtFile_C6_Bild, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=TRUE, rnaseq=TRUE)$es.obs
View(RAS_C6_gsva)

RAS_CGP_gsva <- gsva(TCGA_Lung_RNAseq, gmtFile_CGP, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=TRUE, rnaseq=TRUE)$es.obs
View(RAS_CGP_gsva)

RAS_C2_bild_gsva <- gsva(TCGA_Lung_RNAseq, gmtFile_C2_Bild, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=TRUE, rnaseq=TRUE)$es.obs
View(RAS_CGP_gsva)

IDs=data.frame(colnames(TCGA_Lung_RNAseq))
length(IDs)
class(IDs)
IDs

classFile_TCGA_Lung = read.table("~/Dropbox/GenomeMedicine/Revisons/RAS_HER2/LUAD_Ras_MuationStatus.txt", sep="\t", header=F, stringsAsFactors=FALSE, row.names=NULL, quote="\"", check.names=FALSE)
dim(classFile_TCGA_Lung) #227

# remove duplicates
classFile_TCGA_Lung = classFile_TCGA_Lung[!duplicated(classFile_TCGA_Lung),]
dim(classFile_TCGA_Lung) 
View(classFile_TCGA_Lung)

#sort by Name
classFile_TCGA_Lung = classFile_TCGA_Lung[order(classFile_TCGA_Lung[,1]),]
classFile_TCGA_Lung

#delate
test <- factor(classFile_TCGA_Lung [,2])
test

#filter only RNA-seq samples
classFile_TCGA_Lung_RNA_seq=merge(IDs,classFile_TCGA_Lung, by=1)
dim(classFile_TCGA_Lung_RNA_seq) #169
classFile_TCGA_Lung_RNA_seq

lung_design2 <- model.matrix(~ factor(classFile_TCGA_Lung_RNA_seq [,2]))
lung_design2
colnames(lung_design2) <- c("ALL", "WTvsM")

#lung_fitC2 <- lmFit(RAS_C2_gsva, lung_design2)
#lung_fitC6 <- lmFit(RAS_C6_gsva, lung_design2)
lung_fit_CGP <- lmFit(RAS_CGP_gsva, lung_design2)
lung_fit_C2_bild <- lmFit(RAS_C2_bild_gsva, lung_design2)
#head(lung_fitC2)
#head(lung_fitC6)
head(lung_fit_CPG)
head(lung_fit_C2_bild)

#lung_fitC2 <- eBayes(lung_fitC2)
#lung_fitC6 <- eBayes(lung_fitC6)
lung_fit_CGP <- eBayes(lung_fit_CGP)
View(lung_fit_CGP)
lung_fit_C2_bild <- eBayes(lung_fit_C2_bild)

View(lung_fit_CGP)
View(lung_fit_C2_bild)

#lung_allGeneSetsC2 <- topTable(lung_fitC2, coef="WTvsM", number=Inf)
View(lung_allGeneSetsC2)

#lung_allGeneSetsC6 <- topTable(lung_fitC6, coef="WTvsM", number=Inf)
View(lung_allGeneSetsC6)

lung_allGeneSets_CGP <- topTable(lung_fit_CGP, coef="WTvsM", number=Inf)
View(lung_allGeneSets_CGP)

lung_allGeneSetsC2_bild <- topTable(lung_fit_C2_bild, coef="WTvsM", number=Inf)
View(lung_allGeneSetsC2_blid)


#write.table(lung_allGeneSetsC2, file = "RAS_C2_GSVA.txt", sep = "\t")
#write.table(lung_allGeneSetsC6, file = "RAS_C6_GSVA.txt", sep = "\t")
write.table(lung_allGeneSets_CGP, file = "RAS_CGP_GSVA.txt", sep = "\t")
write.table(lung_allGeneSetsC2_bild, file = "RAS_C2_bild_GSVA.txt", sep = "\t")
```

Endometrial
```{r}
TCGA_Endo_RNAseq <- read.delim("~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/RNA_Seq_Files/PANCAN12.IlluminaHiSeq_RNASeqV2.geneExp.tumor_whitelist_Endometrial_nodup", row.names=1, check.names=FALSE, stringsAsFactors=F)

TCGA_Endo_RNAseq=as.matrix(TCGA_Endo_RNAseq)

head(TCGA_Endo_RNAseq)
dim(TCGA_Endo_RNAseq) # 20499 x 323 RNA-seq samples
row.names(TCGA_Endo_RNAseq)
colnames(TCGA_Endo_RNAseq)

head(gmtFileC2_cp)
length(gmtFileC2_cp)

#removed values with constant expression
Endo_C2_gsva <- gsva(TCGA_Endo_RNAseq, gmtFileC2_cp, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=TRUE, rnaseq=TRUE)$es.obs
View(Endo_C2_gsva)

IDs=data.frame(colnames(TCGA_Endo_RNAseq))
length(IDs)
class(IDs)
IDs

classFile_TCGA_Endo = read.table("~/Documents/ThesisWork/GSOA_Manuscript/GSOA_Files/RNA_Seq_Files/TCGA_UCEC_Serous_new.txt", sep="\t", header=F, stringsAsFactors=FALSE, row.names=NULL, quote="\"", check.names=FALSE)
dim(classFile_TCGA_Endo) #360

# remove duplicates
classFile_TCGA_Endo = classFile_TCGA_Endo[!duplicated(classFile_TCGA_Endo),]
dim(classFile_TCGA_Endo) 
View(classFile_TCGA_Endo)

#sort by Name
classFile_TCGA_Endo = classFile_TCGA_Endo[order(classFile_TCGA_Endo[,1]),]
classFile_TCGA_Endo

#filter only RNA-seq samples
classFile_TCGA_Endo_RNA_seq=merge(IDs,classFile_TCGA_Endo, by=1)
dim(classFile_TCGA_Endo_RNA_seq) #323
classFile_TCGA_Endo_RNA_seq

Endo_design2 <- model.matrix(~ factor(classFile_TCGA_Endo_RNA_seq [,2]))
Endo_design2
colnames(Endo_design2) <- c("ALL", "SvsNS")

Endo_fitC2 <- lmFit(Endo_C2_gsva, Endo_design2)
head(Endo_fitC2)
Endo_fitC2 <- eBayes(Endo_fitC2)
Endo_allGeneSetsC2 <- topTable(Endo_fitC2, coef="SvsNS", number=Inf)
View(Endo_allGeneSetsC2)

write.table(Endo_allGeneSetsC2, file = "Endo_C2_GSVA.txt", sep = "\t", col.names = NA)
```











Leukemia example from GSVA vinettge (microarray data)
```{r}
data(c2BroadSets)
length(c2BroadSets)
c2BroadSets

cacheDir <- system.file("extdata", package="GSVA")
cachePrefix <- "cache4vignette_"
file.remove(paste(cacheDir, list.files(cacheDir, pattern=cachePrefix), sep="/"))

data(leukemia)
#RMA-processed expresson values, background corrected 
leukemia_eset
experimentData(leukemia_eset) #empty
featureData(leukemia_eset) # empty
sampleNames(leukemia_eset)
exprs(leukemia_eset)
varLabels(leukemia_eset)
table(leukemia_eset$subtype)

filtered_eset <- nsFilter(leukemia_eset, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE, feature.exclude="^AFFX")
head(filtered_eset)
leukemia_filtered_eset <- filtered_eset$eset

cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets, min.sz=10, max.sz=500, verbose=TRUE)$es.obs, dir=cacheDir, prefix=cachePrefix)
View(leukemia_es)

adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

DEgeneSetsdesign <- model.matrix(~ factor(leukemia_es$subtype))
DEgeneSetsdesign

colnames(DEgeneSetsdesign) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, DEgeneSetsdesign)
head(fit)
fit <- eBayes(fit)
View(fit)
allGeneSets <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgeneSets <- topTable(fit, coef="MLLvsALL", number=Inf, p.value=adjPvalueCutoff, adjust="BH")
View(DEgeneSets)
res <- decideTests(fit, p.value=adjPvalueCutoff)
View(res)
View(res)
```

RNA-seq data example from the GSVA vinette (compares RNA-seq to microarray)
```{r}
data(commonPickrellHuang)
pickrellCountsYaleCQNcommon_eset

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset), featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset), sampleNames(pickrellCountsArgonneCQNcommon_eset)))

canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
+ grep("^REACTOME", names(c2BroadSets)),
+ grep("^BIOCARTA", names(c2BroadSets)))]
length(canonicalC2BroadSets)

data(genderGenesEntrez)
MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c2"), setName="MSY")
MSY

XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(), collectionType=BroadCollection(category="c2"), setName="XiE")
XiE

canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))
canonicalC2BroadSets

esmicro <- gsva(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, rnaseq=FALSE, parallel.sz=1)$es.obs
dim(esmicro)
View(esmicro)


esrnaseq <- gsva(pickrellCountsArgonneCQNcommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, rnaseq=TRUE, parallel.sz=1)$es.obs
dim(esrnaseq)
class(esrnaseq)
View(esrnaseq)
esrnaseq

library(edgeR)
data(annotEntrez220212)
head(annotEntrez220212)

cpm <- cpm(exprs(pickrellCountsArgonneCQNcommon_eset))
dim(cpm)

common <- intersect(rownames(cpm), rownames(annotEntrez220212))
length(common)
head(common)

rpkm <- sweep(cpm[common, ], 1, annotEntrez220212[common, "Length"] / 10^3, FUN="/")
dim(rpkm)

dim(huangArrayRMAnoBatchCommon_eset[rownames(rpkm), ])

corsrowsgene <- sapply(1:nrow(huangArrayRMAnoBatchCommon_eset[rownames(rpkm), ]), function(i, expmicro, exprnaseq) cor(expmicro[i, ], exprnaseq[i, ], method="pearson"), exprs(huangArrayRMAnoBatchCommon_eset[rownames(rpkm), ]), log2(rpkm+0.1))
names(corsrowsgene) <- rownames(rpkm)
corsrowsgs <- sapply(1:nrow(esmicro), function(i, esmicro, esrnaseq) cor(esmicro[i, ], esrnaseq[i, ], method="spearman"), exprs(esmicro), exprs(esrnaseq))
names(corsrowsgs) <- rownames(esmicro)

View(corsrowsgs
```

