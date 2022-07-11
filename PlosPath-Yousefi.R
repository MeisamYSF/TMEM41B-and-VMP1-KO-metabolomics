library(ggplot2)
setwd("../Results/Scripts/DENV-TMEM-VMP")

###Read data:

Metabolome <- read.csv(file ="MetabolitesRaw.csv",header = T,quote="")
Metabolome <- Metabolome[,1:10]
Metabolome_mat <- as.matrix(Metabolome[,-1])
class(Metabolome_mat) <- "numeric"
row.names(Metabolome) <- Metabolome[,1]
Metabolome <- Metabolome[,-1]

###PCA:

PCAmet <- prcomp(Metabolome_mat)
PCArot <- data.frame(PCAmet$rotation)
PCArot$Group <- c(rep("WT",3),rep("TMEM41B KO",3),rep("VMP1 KO",3))
Summary(PCAmet)
ggplot(PCArot, aes(PC1,PC2,color=Group)) + geom_point(size=4) + theme_bw()

###LFC Plots:

Metabolome <- read.csv(file ="MetabolitesRaw.csv",header = T,quote="")
Metabolome <- Metabolome[,c(1,13,14)]
row.names(Metabolome) <- Metabolome[,1]
Metabolome <- Metabolome[,-1]
Metabolome$Group <- c(rep("FA",12),rep("TG",27),rep("Sphingolipids",7),rep("Glycerophospholipids",22),rep("Diacylglycerol",6),rep("Acylcarnitine",13))
plotting for LFC TMEM vs WT:
ggplot(data = Metabolome, aes(Group,FC_TMEM41BKO_vs_WT, color=Group)) + geom_boxplot() + geom_point() + theme_classic()
plotting for LFC VMP vs WT:
ggplot(data = Metabolome, aes(Group,FC_VMP1KO_vs_WT, color=Group)) + geom_boxplot() + geom_point() + theme_classic()
