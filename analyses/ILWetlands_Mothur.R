################################################################################
#                                                                              #
#  IL Wetlands Mothur Output Check                                             #
#                                                                              #
################################################################################
#                                                                              #
#  Written by: Mario Muscarella                                                #
#                                                                              #
#	Last update: 2015/05/18                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes:                                                                       #
#                                                                              #
# Dependencies: 2) R v3.1.3 (Smooth Sidewalk), or newer                        #
#               2) Vegan v2.2-1                                                #
#               2) MothurTools.R (2015/04/04)                                  #
#                                                                              #
# Issues:                                                                      #
#         1.                                                                   #
#                                                                              #
# Recent Changes:                                                              #
#         1.                                                                   #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/IL_Wetlands/")
source("./bin/MothurTools.R")
require("vegan")
se <- function(x, ...){sd(x, ...)/sqrt(length(na.omit(x)))}

# Define Inputs
# Design = general design file for experiment
# shared = OTU table from mothur with sequence similarity clustering
# design <- ""
shared <- "./data/WL.final.shared"

# Import Design
# design <- read.delim(design, header=T, row.names=1)

# Import Shared Files
OTUs <- read.otu(shared = shared, cutoff = "0.03")         # 97% Similarity

# Remove OTUs with less than two occurances across all sites
OTUs <- Pond97[, which(colSums(OTUs) >= 2)]

# Coverage Stats
coverage <- rowSums(OTUs)
summary(coverage)

# Make Presence Absence Matrices
OTUsPA <- (OTUs > 0) * 1

rich <- rowSums(OTUsPA)
summary(rich)

# Make Relative Abundence Matrices
OTUsREL <- OTUs
for(i in dim(OTUs)[1]){
  OTUs[i,]<-OTUs[i,]/sum(OTUs[i,])
}
