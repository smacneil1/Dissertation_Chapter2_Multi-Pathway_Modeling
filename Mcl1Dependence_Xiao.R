# This scripts analyzes the difference bettwen MCL-1 dependent and independent cell lines 

setwd("~/Documents/Multipathway_Modeling/")
source("http://www.r-statistics.com/wp-content/uploads/2010/02/Barnard.R.txt")
source("http://www.r-statistics.com/wp-content/uploads/2012/01/source_https.r.txt")
source_https("https://raw.github.com/talgalili/R-code-snippets/master/Barnard.R")
MCLdata=read.csv("BH3_Profiling/Mcl1_Xiao_Data.txt",header = T,sep='\t')[1:23,1:5]
MCLdata_ccle=read.csv("BH3_Profiling/Mcl1_Xiao_Data_ccle.txt",header = T,sep='\t')[1:30,1:5]



View(MCLdata)
rownames(MCLdata)

MCLdata[1:13,4]

pdf("~/Documents/Multipathway_Modeling/BH3_Profiling/Mcl1_Plots.pdf")

boxplot(MCLdata_ccle[1:19,4],MCLdata_ccle[20:30,4], ylab="% Viability upon treatment (MCL-1 siRNA) CCLE",names = c("AKT Phenotype \n n=19", " EGFR Phenotype \n n=11"), col = c("red", "grey"), main="Mcl-1 Dependence Across Growth Phenotypes", ylim=c(0,100))
text(1,2,"p-value = 0.32")
t.test(MCLdata_ccle[1:19,4],MCLdata_ccle[20:30,4])


boxplot(MCLdata[1:12,4],MCLdata[13:23,4], ylab="% Viability upon treatment (MCL-1 siRNA)",names = c("AKT Phenotype \n n=12", " EGFR Phenotype \n n=10"), col = c("red", "grey"), main="Mcl-1 Dependence Across Growth Phenotypes", ylim=c(0,100))
text(1,2,"p-value = 0.43")
t.test(MCLdata[1:14,4],MCLdata[15:23,4])
dev.off()



MclDependence_before<- matrix(c(9,2,5,7), nrow=2, dimnames=list(c("EGFR", "AKT"), c("Mcl Dependent", "Non-Mcl Dependent")))
MclDependence_after<- matrix(c(4,7,6,5), nrow=2, dimnames=list(c("EGFR", "AKT"), c("Mcl Dependent", "Non-Mcl Dependent")))
MclDependence_ccle<- matrix(c(4,6,11,8), nrow=2, dimnames=list(c("EGFR", "AKT"), c("Mcl Dependent", "Non-Mcl Dependent")))

View(Mcl)
fisher.test(MclDependence_after)
fisher.test(MclDependence_ccle)

fisher.test(MclDependence_after, alternative = "greater")
View(MclDependence)
Barnard(MclDependence_after)
Barnard(MclDependence_ccle)


TeaTasting <-
  matrix(c(3, 1, 1, 3),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))
View(TeaTasting)
fisher.test(TeaTasting, alternative = "greater")

Convictions <-
  matrix(c(2, 10, 15, 3),
         nrow = 2,
         dimnames =
           list(c("Dizygotic", "Monozygotic"),
                c("Convicted", "Not convicted")))
View(Convictions)



