################################################################################
#                                                                              #
#	Henry County Mitigation Project 2013 : Microbial Community Charactorization               #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Jay Lennon (2014/02/03)                                        #
# Modified by: M. Muscarella (2014/02/03)                                      #
# Modified by: Ariane Peralta                                                  #
#	Last update: 2014/02/04                                                    #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/IL_Wetlands")
se <- function(x){sd(x)/sqrt(length(x))}

# Code Dependencies
source("DiversityFunctions.r")  
require(vegan)                               
require(reshape)
require(BiodiversityR)
require("ecodist")

write.table(OTUsREL, file = "WLbact_STDnew.csv", sep = ",", col.names = NA)
read.table("WLbact_STDnew.csv", header = TRUE, sep = ",", row.names = 1)

#I'm not sure why the OTUsREL didn't get standardized - so the WLbact_STDnew.csv isn't standardized. Please help :)

new.data <- read.csv("./WLbact_STDnew.csv", header=TRUE)
df.env <- read.csv("./data/WL_plant_soil.csv", header=TRUE)

new.data.env <- merge(df.env, new.data) #but we need to get rid of the ID column in the new.data dataframe

#rows,col
adonis = adonis(new.data.env[,-c(1:17)] ~ Treatment, method = "bray", data = new.data, perm=1000)
adonis

Simper <- simper(sampleREL.dist, group = Treatment)
summary(Simper)


library(ggplot2)
df <- new.data[,-c(1:6)]	# bacteria
df.env <- new.data[,c(1:6)]	# ENV
df.mds <- metaMDS(df, k=2, trymax=50, zerodist="add")

treatments <- new.data$Treatment
levels(treatments) <- c("BallBurlapped", "Bareroot", "Seedling", "Acorn", "Seedbank", "Reference")

points <- cbind(as.data.frame(df.mds$points), treatments)
L.centroids <- melt(points, id="treatments", measure.vars = c("MDS1", "MDS2")) 
centroids <- cast(L.centroids, variable ~ treatments, mean)


df <- data.frame(new.data[,3],new.data[,2], df.mds$points[,1], df.mds$points[,2])
str(df)

levels(new.data$Treatment) <- c("BallBurlapped", "Bareroot", "Seedling", "Acorn", "Seedbank", "Reference")

names(df) <- c("Treatment", "Plot", "Axis1", "Axis2")
df$Treatment <- factor(df$Treatment, levels=c("BallBurlapped", "Bareroot", "Seedling", "Acorn", "Seedbank", "Reference"))
myColors <- c("#FFF000", "#CCFF00", "#33CC33", "#339933", "#336633", "#FF9933") #pick new
names(myColors) <- levels(df$Treatment)
colScale <- scale_colour_manual(values = myColors)
p1<-ggplot(df,aes(x=Axis1,y=Axis2,label=TRUE))
p2 <- p1+geom_point(aes(colour=Treatment), size=5)
# p3 <-p2 + geom_text(aes(label=new.data$taxon), size=4) #added labels
# p3 + colScale
p2 + colScale

which(df$Axis1 > 0.05)
#[1] 11 42 49 56


library(ggplot2)
df <- new.data[,-c(1:8)]	# bacteria
df.env <- new.data[,c(1:8)]	# ENV
df.mds <- metaMDS(df, k=2, trymax=50, zerodist="add")
df <- data.frame(new.data[,5],new.data[,4], df.mds$points[,1], df.mds$points[,2])
str(df)

levels(new.data$Treatment) <- c("BallBurlapped", "Bareroot", "Seedling", "Acorn", "Seedbank", "Reference")

names(df) <- c("Treatment", "Plot", "Axis1", "Axis2")
df$Treatment <- factor(df$Treatment, levels=c("BallBurlapped", "Bareroot", "Seedling", "Acorn", "Seedbank", "Reference"))
myColors <- c("#FFF000", "#CCFF00", "#33CC33", "#339933", "#336633", "#FF9933") #pick new
names(myColors) <- levels(df$Treatment)
colScale <- scale_colour_manual(values = myColors)
p1<-ggplot(df,aes(x=Axis1,y=Axis2,label=TRUE))
p2 <- p1+geom_point(aes(colour=Treatment), size=5)
#p3 <-p2 + geom_text(aes(label=new.data$RunID), size=4) #added labels
#p3 + colScale
p2 + colScale

#to figure out outliers
which(df$Axis1 > 0.05)
#[1] 11 12 13 45 53 61

# Create Distance Matrix
samplePA.dist <- vegdist(t(dataPA),method="bray")
sampleREL.dist <- vegdist(t(dataREL),method="bray")

# Principal Coordinates Analysis
EC_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE) 
  # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
  # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

# Responder Analysis Based on PCoA 
pcoaS <- add.spec.scores(EC_pcoa,t(dataREL),method="cor.scores",Rscale=TRUE,
  scaling=1,multi=1) 
  # retrieves correlation coefficient for each taxon's relative 
  # abudnace with respect to PCoA coordinates (k = 3)
  
# PCoA Axis 1 Responders
cor_spp_a1 <- cbind(EC_tax[,3],pcoaS$cproj[,1]) 
  # creates matrix of taxonomy and correlations (r)
colnames(cor_spp_a1)[7] <- "Corr" 
  # renames correlation column name
cor_spp_a1 <- cor_spp_a1[order(cor_spp_a1[,7],decreasing=TRUE),] 
  # sorts based on r
responders_a1 <- cor_spp_a1[a1 <- abs(cor_spp_a1[,7])>0.7,] 
  # subset based on r > |0.7| 

# PCoA Axis 2 Responders
cor_spp_a2 <- cbind(EC_tax[,3],pcoaS$cproj[,2]) 
  # creates matrix of taxonomy and correlations(r) axis 2
colnames(cor_spp_a2)[7] <- "Corr" 
  # renames correlation column name
cor_spp_a2 < -cor_spp_a2[order(cor_spp_a1[,7],decreasing=TRUE),] 
  # sorts based on r
responders_a2 <- cor_spp_a2[a2 <- abs(cor_spp_a2[,7])>0.7,] 
  # subset based on r > |0.7| 

# Percent Variacne Explained Using PCoA 
Eig_sum <- sum(pcoaS$eig) 
  # sum of k eigenvalues
eig1 <- pcoaS$eig[1] 
  # eigenvalue for PCoA 1
explainvar1 <- pcoaS$eig[1]/Eig_sum 
  # eigenvalue for PCoA 1 divided by sum of k eigenvalues
explainvar1 
  # percent variance explained by PCoA 1
eig2 <- pcoaS$eig[2] 
  # eigenvalue for PCoA 2
explainvar2 <- pcoaS$eig[2]/Eig_sum 
  # eigenvalue for PCoA 2 divided by sum of k eigenvalues
explainvar2 
  # percent variance explained by PCoA 2
eig3 <- pcoaS$eig[3] 
  # eigenvalue for PCoA 3
explainvar3 <- pcoaS$eig[3]/Eig_sum 
  # eigenvalue for PCoA 3 divided by sum of k eigenvalues
explainvar3 
  # percent variance explained by PCoA 3

## PCoA Plots ##################################################################
################################################################################

# Organize Samples for Plotting
# Plot Codes:
# "St" = station; 1 and 2 = African (dry) slope; 5 and 6 = European (wet) slope 
# "Loc" = locations; 10 replicate plots on each station, identified (1, 2, 5, 6)

# Stripping PCoA Loadings out of cmdscale Output
pcoap <- as.data.frame(pcoaS$points)
site <- rep(NA,76)
site[c(1,3,5,7,9,11,13,15,17,19)] <- "AF1_DNA"
site[c(21,23,24,25,27,29,31,33,35)] <- "AF2_DNA"
site[c(37,39,41,43,45,47,49,51,53,55)] <- "EU1_DNA"
site[c(57,59,61,63,65,67,69,71,73,75)] <- "EU2_DNA"
site[c(2,4,6,8,10,12,14,16,18,20)] <- "AF1_RNA"
site[c(22,26,28,30,32,34,36)] <- "AF2_RNA"
site[c(38,40,42,44,46,48,50,52,54,56)] <- "EU1_RNA"
site[c(58,60,62,64,66,68,70,72,74,76)] <- "EU2_RNA"
pcoap <- cbind(pcoap, site)

#--> AF: Station 1 - DNA
AF1_DNA_1 <- pcoap[c(1,3,5,7,9,11,13,15,17,19),1]             #PCoA 1 scores
AF1_DNA_2 <- pcoap[c(1,3,5,7,9,11,13,15,17,19),2]             #PCoA 2 scores 
#--> AF: Station 2 - DNA
AF2_DNA_1 <- pcoap[c(21,23,24,25,27,29,31,33,35),1]           #PCoA 1 scores  
AF2_DNA_2 <- pcoap[c(21,23,24,25,27,29,31,33,35),2]           #PCoA 2 scores 
#--> EU: Station 5 - DNA
EU1_DNA_1 <- pcoap[c(37,39,41,43,45,47,49,51,53,55),1]           #PCoA 1 scores
EU1_DNA_2 <- pcoap[c(37,39,41,43,45,47,49,51,53,55),2]           #PCoA 2 scores 
#--> EU: Station 6 - DNA
EU2_DNA_1 <- pcoap[c(57,59,61,63,65,67,69,71,73,75),1]        #PCoA 1 scores
EU2_DNA_2 <- pcoap[c(57,59,61,63,65,67,69,71,73,75),2]        #PCoA 2 scores 
#--> AF: Station 1 - RNA
AF1_RNA_1 <- pcoap[c(2,4,6,8,10,12,14,16,18,20),1]            #PCoA 1 scores
AF1_RNA_2 <- pcoap[c(2,4,6,8,10,12,14,16,18,20),2]            #PCoA 2 scores 
#--> AF: Station 2 - RNA
AF2_RNA_1 <- pcoap[c(22,26,28,30,32,34,36),1]                 #PCoA 1 scores  
AF2_RNA_2 <- pcoap[c(22,26,28,30,32,34,36),2]                 #PCoA 2 scores  
#--> EU: Station 5 - RNA
EU1_RNA_1 <- pcoap[c(38,40,42,44,46,48,50,52,54,56),1]        #PCoA 1 scores 
EU1_RNA_2 <- pcoap[c(38,40,42,44,46,48,50,52,54,56),2]        #PCoA 2 scores   
#--> EU: Station 6 - RNA
EU2_RNA_1 <- pcoap[c(58,60,62,64,66,68,70,72,74,76),1]        #PCoA 1 scores 
EU2_RNA_2 <- pcoap[c(58,60,62,64,66,68,70,72,74,76),2]        #PCoA 2 scores  

# Plot Parameters
par(mfrow=c(1,1), mar=c(5,5,1,1)) 
x <- c(min(pcoap[,2])+min(pcoap[,2])*0.25,max(pcoap[,2])+max(pcoap[,2])*0.1)
y <- c(min(pcoap[,1])+min(pcoap[,1])*0.1,max(pcoap[,1])+max(pcoap[,1])*0.1)

# Initiate Plot
plot(pcoap[,2], pcoap[,1], xlab="PCoA Axis 2", ylab="PCoA Axis 1", 
  xlim=rev(range(x)),ylim= y, pch=16, cex=2.0, type='n',xaxt="n",yaxt="n", 
  cex.lab=1.5, cex.axis=1.2) 
    # Creates PCoA byplot for PCoA axis 1 and PCoA axis 2; adds axes labels, 
    # Adjusts axes lengths and size
    
  axis(side=1,at=c(-0.4,-0.2,0,0.2), las=1)    # adds x-axis ticks
  axis(side=2,at=c(-0.2,0,0.2,0.4), las=1)     # adds y-axis ticks
  segments(-1500, -0, 1500, 0, lty="dotted")   # adds horizontal reference line
  segments(0, -1500, 0, 1500, lty="dotted")    # adds vertical reference line
  
  points(AF1_DNA_2,AF1_DNA_1,pch=21,cex=2.0,col="black",bg="brown3",lwd=2)
  points(AF2_DNA_2,AF2_DNA_1,pch=21,cex=2.0,col="black",bg="brown3",lwd=2)
  points(EU1_DNA_2,EU1_DNA_1,pch=21,cex=2.0,col="black",bg="green3",lwd=2)
  points(EU2_DNA_2,EU2_DNA_1,pch=21,cex=2.0,col="black",bg="green3",lwd=2)
  points(AF1_RNA_2,AF1_RNA_1,pch=22,cex=2.0,col="black",bg="brown3",lwd=2)
  points(AF2_RNA_2,AF2_RNA_1,pch=22,cex=2.0,col="black",bg="brown3",lwd=2)
  points(EU1_RNA_2,EU1_RNA_1,pch=22,cex=2.0,col="black",bg="green3",lwd=2)
  points(EU2_RNA_2,EU2_RNA_1,pch=22,cex=2.0,col="black",bg="green3",lwd=2)
    # adds PCoA scores for each to existing graph
  
  ordiellipse(cbind(pcoap[,2], pcoap[,1]), pcoap$site, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=FALSE)  
  #legend(-0.25, -0.08,"AF-1-DNA", cex=2.0,col="red",pch=16,bty="n")  
  #legend(-0.25, -0.12,"AF-2-DNA", cex=2.0,col="pink",pch=16,bty="n")
  #legend(-0.25, -0.12,"EU-1-DNA", cex=1.5,col="darkgreen",pch=16,bty="n")
  #legend(-0.25, -0.12,"EU-2-DNA", cex=1,col="chartreuse4",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="red",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="pink",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="red",pch=16,bty="n")
  #legend(-0.25, -0.12,"AF-1-RNA", cex=1,col="pink",pch=16,bty="n")
    # adds legends to existing graph
  box(lwd=2)

## Hypothesis Testing ##########################################################
################################################################################

site <- rep(NA,76)
site[c(1,3,5,7,9,11,13,15,17,19)] <- "AF1_DNA"
site[c(21,23,24,25,27,29,31,33,35)] <- "AF2_DNA"
site[c(37,39,41,43,45,47,49,51,53,55)] <- "EU1_DNA"
site[c(57,59,61,63,65,67,69,71,73,75)] <- "EU2_DNA"
site[c(2,4,6,8,10,12,14,16,18,20)] <- "AF1_RNA"
site[c(22,26,28,30,32,34,36)] <- "AF2_RNA"
site[c(38,40,42,44,46,48,50,52,54,56)] <- "EU1_RNA"
site[c(58,60,62,64,66,68,70,72,74,76)] <- "EU2_RNA"

station <- rep(NA, 76)
station[1:20] <- "AF1"
station[21:36] <- "AF2"
station[37:56] <- "EU1"
station[57:76] <- "EU2"

molecule <- rep(NA, 76)
molecule[c(1,3,5,7,9,11,13,15,17,19,21,23,24,25,27,29,31,33,35,
  37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75)] <- "DNA"
molecule[c(2,4,6,8,10,12,14,16,18,20,22,26,28,30,32,34,36,
  38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76)] <- "RNA"
  
slope <- rep(NA, 76)
slope[1:36] <- "AF"
slope[37:76] <- "EU"

#conduct an Analysis of Similarities (ANOSIM)
#if consulting dissimilarity ranks, use values from 50% column
Anosim <- anosim(sampleREL.dist, grouping = site, permutations = 1000, distance = "bray")
summary(Anosim)

#conduct a Similarity Percentages analysis
#an overall average proportional dissimilarity between treatments is not reported, but it can be obtained by summing the values in the contribution column
#ratio is the average proportional contribution of the respective OTU divided by its SD
Simper <- simper(sampleREL.dist, group = site)
summary(Simper)


# anosim performs an analysis of simliarities on a distance matrix
# The default distance measure is bray-curtis, but other measures (Chao, Jaccard, Euclidean) can be used when specified
# Example specifying distance measure: anosim(lter3REL.dist, grouping=c(rep(1,3), rep(2,3)), distance="euclidean", permutations=1000)
# The grouping argument used here specifies what group (1 or 2) each sample belongs to (here the first 3 samples belong to group 1)
anosim(sampleREL.dist, grouping=site, permutations=1000)


# mrpp performs a multi-response permutation procedure to test for within vs. among group differences
# The default distance measure is euclidean, but other measures (Chao, Jaccard, Bray-Curtis) can be used when specified
mrpp(t(sampleREL.dist), grouping=site, distance="bray", permutations=1000)

# adonis runs a permanova (Created by Marti J. Anderson) this is very similar to ANOVA but for multivariate data. You can make very complex
# experimental designs with it.
# The default distance measure is bray-curtis, but other measures (Chao, Jaccard, Euclidean) can be used when specified
adonis(sampleREL.dist ~ molecule + site , method="bray", permutations=1000)


*******
data.bact<- read.csv("~/Documents/Projects_ACTIVE/Wetlands_Illinois/2013_HenryCountyMitigation/WLbact_only.csv", header=TRUE)
data.env<- read.csv("~/Documents/Projects_ACTIVE/Wetlands_Illinois/2013_HenryCountyMitigation/WL_soil.csv", header=TRUE)

data(data.bact)
data(data.env)
(sim <- with(data.env, simper(data.bact, Treatment)))
summary(sim, ordered=TRUE)

data(dune)
data(dune.env)
(sim <- with(dune.env, simper(dune, Management)))
summary(sim)

simper(comm, group,  ...)
## S3 method for class 'simper'
summary(object, ordered = TRUE, 
     digits = max(3, getOption("digits") - 3), ...)
     