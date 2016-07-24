#                                                                              #
#	Long-term fert and mowing at WRC:  Community Charactorization                #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Ariane Peralta (2016/05/13)                                      #
# Modified by:                                                                 #
#                                                                              #
#                                                                              #
################################################################################

#Code for multivariate analyses - AP work in progress - used file WRC_Importance_final.csv b/c had to hand enter 4 samples

# Setup Work Environment
rm(list=ls())
setwd("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology")
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}
library(vegan)
PCC <- read.csv ("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology/WRC_Importance_final.csv", header=TRUE)
labels(PCC)
soil <- read.csv("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology/data_original/WRC_Soil_Data.csv", header=TRUE)
labels(soil)


treatments <- PCC$treatment
levels(treatments) <- c("UM/UF", "UM/F", "M/UF", "M/F")

###########################
# Simple Hypothesis Testing
###########################

adonis = adonis(PCC[,-c(1:9)] ~ Fertilizer*Mowing*Year, method = "bray", data = PCC, perm=1000)
adonis

adonis2 = adonis(PCC[,-c(1:9)] ~ Fertilizer*Mowing*Year+(1|BLOCK/QUADRAT..), method = "bray", data = PCC, perm=1000)
adonis2
#adonis1 and adonis2 models were the same

WL.simper <- simper(PCC[,-c(1:9)], group = PCC$treatment)
simper <- summary(WL.simper)
simper

#adonis(formula = PCC[, -c(1:9)] ~ Fertilizer * Mowing * Year,      data = PCC, permutations = 1000, method = "bray") 

#Permutation: free
#Number of permutations: 1000

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
#Fertilizer               1    10.988 10.9883  56.034 0.04454 0.000999 ***
#Mowing                   1    22.378 22.3785 114.118 0.09070 0.000999 ***
#Year                     1    16.649 16.6492  84.902 0.06748 0.000999 ***
#Fertilizer:Mowing        1     2.778  2.7775  14.164 0.01126 0.000999 ***
#Fertilizer:Year          1     1.345  1.3446   6.857 0.00545 0.000999 ***
#Mowing:Year              1     5.316  5.3162  27.110 0.02155 0.000999 ***
#Fertilizer:Mowing:Year   1     0.589  0.5894   3.006 0.00239 0.003996 ** 
#Residuals              952   186.687  0.1961         0.75664             
#Total                  959   246.731                 1.00000             
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#visualize plant community composition based on PCoA by year x treatment ([centroid of 8 replicate block x 3 quadrats] for each year x treatment)

# Setup Work Environment
rm(list=ls())
setwd("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology")
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}
ci <- function(x, ...){1.96 * sd(x,na.rm = TRUE)}

# Code Dependencies
library(MASS)
library(nlme)
library(reshape2)
library(vegan)
library(reshape)
library(lme4)
library(ggplot2)


se <- function(x){sd(x)/sqrt(length(x))}
new.data <- read.csv("~/Dropbox/WRC_Project_Writing/Manuscript_WRC_Plant_Ecology/WRC_Importance_final.csv", header=TRUE)
sampleREL.dist <- vegdist(decostand((new.data[,-c(1:9)]),method="log"),method="bray")
WL_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE)
explainvar1 <- round(WL_pcoa$eig[1]/sum(WL_pcoa$eig)*100,2)
explainvar2 <- round(WL_pcoa$eig[2]/sum(WL_pcoa$eig)*100,2)
explainvar3 <- round(WL_pcoa$eig[3]/sum(WL_pcoa$eig)*100,2)
explainvar1
explainvar2
explainvar3
pcoap <- merge(as.data.frame(WL_pcoa$points),new.data$treatment, by=0,all.x=T)
rownames(pcoap) <- rownames(WL_pcoa$points)
treatments <- new.data$treatment
levels(treatments) <- c("UM/UF", "UM/F", "M/UF", "M/F")
points <- cbind(as.data.frame(WL_pcoa$points), treatments)
L.centroids <- melt(points, id="treatments", measure.vars = c("V1", "V2", "V3"))
centroids <- cast(L.centroids, treatments ~ variable, mean)
centroids <- cast(L.centroids, treatments ~ variable, fun.aggregate=c(mean,se))
write.csv(centroids, file="newdataClassified_genus_OTUsP.centroids.csv")

------------------------------------------------------------------------
  
  par(mar=c(6,6,1,1), oma=c(3,1,1,1)+0.1 )
x.dim <- c(min(pcoap$V1)-(max(pcoap$V1)*0.1),max(pcoap$V1)+(max(pcoap$V1)*0.1))
x.dim <- c(-1, 1)
y.dim <- c(min(pcoap$V2)-(max(pcoap$V2)*0.1),max(pcoap$V2)+(max(pcoap$V2)*0.1))
plot(pcoap$V1, pcoap$V2, xlab=paste("PCoA Axis 1 (",explainvar1, "%)", sep=""), ylab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""), xlim=x.dim, ylim=y.dim, pch=16, cex=2.0, type="n",xaxt="n",yaxt="n", cex.lab=1.5, cex.axis=1.2)
axis(side=1, las=1)
axis(side=2, las=1)
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")
box(lwd=2)
points(pcoap$V1, pcoap$V2, pch=19, cex=2.0, bg="gray", col="gray")
text(pcoap$V1, pcoap$V2, labels=pcoap$y, pos=3)
ordiellipse(cbind(pcoap$V1, pcoap$V2), pcoap$y, kind="se", conf=0.95, lwd=2, draw = "polygon", col="gray", border = "black", label=TRUE, cex=2)
levels(treatments) <- c("UM/UF", "UM/F", "m/UF", "MF")
myColors <- c("#FFF000", "#CCFF00", "#33CC33", "#339933")
names(myColors) <- levels(treatments)
