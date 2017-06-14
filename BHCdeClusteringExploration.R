# BHCdeClusteringExploration.R
# R version 3.3.1 (2016-06-21)
# June 13, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread

#-----------------------------------------------------------------------
library(BHC)
library(bnlearn)
#library(stringr)
#-----------------------------------------------------------------------

#source("https://bioconductor.org/biocLite.R")
#biocLite("BHC")

#require(BHC)



gdata <- read.csv(file = "largeDE.csv", row.names = 1)

RNA <- as.data.frame(t(gdata))

discRNA <- discretize(RNA, method = "quantile", breaks = 3)

for (i in 1:dim(discRNA)[2]){
  levels(discRNA[, i]) <- c(-1, 0, 1)
  discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
}

discRNA <- t(discRNA)



##DATA DIMENSIONS
nDataItems <- nrow(discRNA)
nFeatures  <- ncol(discRNA)

itemLabels <- rownames(discRNA)


hct <- bhc(discRNA, itemLabels, numReps=2)
WriteOutClusterLabels(hct, "labels.txt", verbose=TRUE)


x <- fread(file = "labels.txt")


#####################################################################

##BUILD SAMPLE DATA AND LABELS
data         <- matrix(0,15,10)
itemLabels   <- vector("character",15)
data[1:5,]   <- 1 ; itemLabels[1:5]   <- "a"
data[6:10,]  <- 2 ; itemLabels[6:10]  <- "b"
data[11:15,] <- 3 ; itemLabels[11:15] <- "c"
timePoints   <- 1:10 # for the time-course case

##DATA DIMENSIONS
nDataItems <- nrow(data)
nFeatures  <- ncol(data)

##RUN MULTINOMIAL CLUSTERING
hc1 <- bhc(data,itemLabels,verbose=TRUE)

##RUN TIME-COURSE CLUSTERING
hc2 <- bhc(data, itemLabels, 0, timePoints, "time-course",
           numReps=1, noiseMode=0, numThreads=1, verbose=TRUE)

##OUTPUT CLUSTER LABELS TO FILE
WriteOutClusterLabels(hc1, "labels.txt", verbose=TRUE)

##FOR THE MULTINOMIAL CASE, THE DATA CAN BE DISCRETISED
newData      <- data[] + rnorm(150, 0, 0.1);
percentiles  <- FindOptimalBinning(newData, itemLabels, transposeData=TRUE, verbose=TRUE)
discreteData <- DiscretiseData(t(newData), percentiles=percentiles)
discreteData <- t(discreteData)
hc3          <- bhc(discreteData, itemLabels, verbose=TRUE)



###################################################
### code chunk number 2: bhc.Rnw:50-53
###################################################
plot(hc1, axes=FALSE)
plot(hc2, axes=FALSE)
plot(hc3, axes=FALSE)
