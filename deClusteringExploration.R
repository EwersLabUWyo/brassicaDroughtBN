# deClusteringExploration.R
# R version 3.3.1 (2016-06-21)
# March 11, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread

#-----------------------------------------------------------------------
#library(maSigPro)
library(bnlearn)
library(stringr)
#-----------------------------------------------------------------------

degenes <- read.csv(file.choose(), row.names = 1)        
# Cluster by sample. 
hr <- hclust(dist(t(degenes)))
plot(hr, xlab = "de")
rect.hclust(hr, k = 6, border = "red")

d <- scale(degenes)
hr <- hclust(dist(t(d)))
plot(hr, xlab = "scaledDE")
rect.hclust(hr, k = 6, border = "red")

# Discretize scaled genes and then cluster.
 
RNA <- degenes

# Convert to data.frame. 
RNA <- as.data.frame(t(RNA))

# Discretize the data. 
discRNA <- discretize(RNA, method = "quantile", breaks = 3)

r <- discRNA

for (i in 1:25){
  levels(r[, i]) <- c(-1, 0, 1)
  r[, i] <- as.numeric(as.character(r[, i]))
}

r = t(r)

# Rewrite rownames of r.
#rownames(r) <- rownames(degenes)
colnames(r) <- colnames(degenes)

f <- as.data.frame(t(r))

# Split rownames into Treatment/Timepoint and Replicate. 
rnaNames <- str_split_fixed(rownames(RNA), '_', 2)

# Create a Timepoint column. 
f$Timepoint <- as.numeric(str_split_fixed(rnaNames[, 1], '', 2)[, 2])

# Create a treatment column named int. 
f$int <- as.factor(str_split_fixed(rnaNames[, 1], '', 2)[, 1])
f <- as.data.table(f)


d <- f[, lapply(.SD, mean), by= .(Timepoint, int) ]
n <- c(paste(rep("D", 12), d$Timepoint[1:12], sep = ""), 
       paste(rep("W", 12), d$Timepoint[1:12], sep = ""))
d <- d[, Timepoint := NULL]
d <- d[, int := NULL]
d <- t(d)
colnames(d) <- n

hr <- hclust(dist(t(d)))
plot(hr, xlab = "disc")
rect.hclust(hr, k = 6, border = "red")

d["Bra001534", ] 
d["Bra001379", ]


sum(d["Bra001534", ]  != d["Bra001379", ])
sum(d["Bra001534", ]  != d["Bra000159", ])


# Split rownames into Treatment/Timepoint and Replicate. 
rnaNames <- str_split_fixed(rownames(RNA), '_', 2)

# Create a Timepoint column. 
discRNA$Timepoint <- as.numeric(str_split_fixed(rnaNames[, 1], '', 2)[, 2])

# Create a treatment column named int. 
discRNA$int <- as.factor(str_split_fixed(rnaNames[, 1], '', 2)[, 1])



source("https://bioconductor.org/biocLite.R")
biocLite("BHC")

require(BHC)
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






itemLabels   <- rep(c("D", "W"), each = 24)
timePoints   <- 1:12 # for the time-course case

##DATA DIMENSIONS
nDataItems <- nrow(data)
nFeatures  <- ncol(data)


hct <- bhc(r, itemLabels, 0, timePoints, "time-course",
           numReps=2, noiseMode=0, numThreads=1, verbose=TRUE)


# Using rnaPheno from network with first 50 genes. 

for (i in 7:56){
  levels(rnaPheno[, i]) <- c(-1, 0, 1)
  rnaPheno[, i] <- as.numeric(as.character(rnaPheno[, i]))
}

for (i in 1:2){
  levels(rnaPheno[, i]) <- c(-2, -1, 0, 1, 2)
  rnaPheno[, i] <- as.numeric(as.character(rnaPheno[, i]))
}

for (i in 4:6){
  levels(rnaPheno[, i]) <- c(-2, -1, 0, 1, 2)
  rnaPheno[, i] <- as.numeric(as.character(rnaPheno[, i]))
}

  levels(rnaPheno[, 3]) <- c(-1, 0, 1)
  rnaPheno[, 3] <- as.numeric(as.character(rnaPheno[, 3]))

  hr <- hclust(dist(t(rnaPheno)))
  plot(hr)
  
  hr <- hclust(dist(rnaPheno))
  plot(hr)
  
  library(cluster)
  set.seed(125)
  gap <- clusGap(rnaPheno, FUN = kmeans, iter.max = 30, K.max = 20, B = 50, verbose=interactive())
  plot(gap, main = "Gap Statistic")
  with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
  
  plot(hr)
  rect.hclust(hr, k = 6, border = "red")
  
  plot(hr)
  rect.hclust(hr, k = 50, border = "red")
  
  cutree(hr, k = 50)
  rnaPheno[which(cutree(hr, k = 50) == 2),]
  
  hr <- hclust(dist(rnaPheno[, 7:56]))
  plot(hr)
  #cutree(hr, k = 50)
  rnaPheno[which(cutree(hr, k = 50) == 2),]
  
  hr <- hclust(dist(t(rnaPheno[, 7:56]))) 
  plot(hr)
  rect.hclust(hr, k = 14, border = "red")
  
  c <- as.data.frame(cutree(hr, k=10))
  rnaPheno[, rownames(c)[which(c==2)]]
  
  g <- rnaPheno
  
  
  
# Don't bind to pheno dataset. Use genes from last network. 
  # Read in RNAseq data. 
  RNA <- read.csv(file = "/Users/mblai/Documents/Thesis/PhenoRNAnetworkBrassica/DEgenesLog.csv", row.names = 1)
  
  # Limit RNA to first 10 genes. 
  #RNA <- RNA[1:50, ]
  
  ## Transpose data. 
  RNA <- t(RNA)
  
  # Convert to data.frame. 
  RNA <- as.data.frame(RNA)
  
  # Discretize the data. 
  discRNA <- discretize(RNA, method = "quantile", breaks = 3)
  
  # Split rownames into Treatment/Timepoint and Replicate. 
  rnaNames <- str_split_fixed(rownames(RNA), '_', 2)
  
  # Create a Timepoint column. 
  discRNA$Timepoint <- as.numeric(str_split_fixed(rnaNames[, 1], '', 2)[, 2])
  
  # Create a treatment column named int. 
  discRNA$int <- as.factor(str_split_fixed(rnaNames[, 1], '', 2)[, 1])
  
  # Order RNA dataframe by Timepoint and int. 
  discRNA <- discRNA[with(discRNA, order(Timepoint, int)), ]
  
  # Remove unneccesary dataframes. 
  #rm(RNA)
  rm(rnaNames)
  
  for (i in 1:43){
    levels(discRNA[, i]) <- c(-1, 0, 1)
    discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
  }
  # -1 is D
  
  hr <- hclust(dist(t(discRNA[, 1:41]))) 
  plot(hr)
  rect.hclust(hr, k = 14, border = "red")
  
  c <- as.data.frame(cutree(hr, k=10))
  d <- discRNA[, rownames(c)[which(c==2)]]
  t <- d[as.character(1:48),]
  
  
## Make sure that there aren't a large number of duplicates when averaged.   
  RNA <- read.csv(file = "largeDE.csv", row.names = 1)
  
  RNA <- as.data.frame(t(RNA))
  
  discRNA <- discretize(RNA, method = "quantile", breaks = 3)
  
  for (i in 1:dim(discRNA)[2]){
    levels(discRNA[, i]) <- c(-1, 0, 1)
    discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
  }
  
  discRNA <- t(discRNA)
  f <- t(discRNA)

  Replicates <- rep(c(1:24), each = 2)
  
  f <- cbind(f, Replicates)
  f <- as.data.table(f)
  
  g <- f[, lapply(.SD, mean), by = "Replicates"]
  g <- g[, -1]
  g <- t(g)
  
  c <- unique(g)
 
  hcg <- hclust(dist(g))
  plot(hcg)
  
  library(cluster)
  set.seed(125)
  gap <- clusGap(g, FUN = kmeans, iter.max = 30, K.max = 150, B = 50, verbose=interactive())
  plot(gap, main = "Gap Statistic")
  with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
  #10
  plot(hcg)
  rect.hclust(hcg, k = 10, border = "red")
  
  trg <- as.data.frame(cutree(hcg, k = 10))
  g[which(cutree(hcg, k = 10) == 7),]
  # So far, the most promising.
 
  
# Use hclust on RNA + genes.
  # Read in phenotype file. 
  Pheno <- read.csv(file = "PhenoBrassicaImp.csv", row.names = 1)
  
  # Rename SM... to get rid of periods. 
  colnames(Pheno)[8] <- "SM"
  
  # Add a column for Time of Day, named TOD.
  Pheno$TOD <- rep(c(7, 11, 15, 19, 23, 3), each = 24, 2)
  
  #### Discretize data. 
  
  # Discretize the phenotype data, excluding Fv.Fm. 
  phenoDisc <- discretize(Pheno[, 3:8], 
                          method = "interval", 
                          breaks = c(5, 5, 3, 5, 5, 5))
  
  # Add INT column to discretized data. 
  phenoDisc$INT <- as.factor(Pheno$Treatment)
  
  # Add Timepoint column to discretized data. 
  phenoDisc$TP <- as.factor(Pheno$Timepoint)
  
  # Add Time of Day to discretized data. 
  phenoDisc$TOD <- as.factor(Pheno$TOD)
  
  # Order Pheno dataframe by Timepoint and int. 
  phenoDisc <- phenoDisc[with(phenoDisc, order(INT, TP)), ]
  
  # Remove unnecesary dataframes.
  rm(Pheno)
  
  for (i in 1:2){
    levels(phenoDisc[, i]) <- c(-2, -1, 0, 1, 2)
    phenoDisc[, i] <- as.numeric(as.character(phenoDisc[, i]))
  }
  
  for (i in 4:6){
    levels(phenoDisc[, i]) <- c(-2, -1, 0, 1, 2)
    phenoDisc[, i] <- as.numeric(as.character(phenoDisc[, i]))
  }
  
  levels(phenoDisc[, 3]) <- c(-1, 0, 1)
  phenoDisc[, 3] <- as.numeric(as.character(phenoDisc[, 3]))
  
  Replicates <- rep(c(1:24), each = 12)
  phenoDisc$Replicates <- Replicates
  
  phenoDisc$INT <- NULL
  phenoDisc$TOD <- NULL
  phenoDisc$TP <- NULL
  
  p <- as.data.table(phenoDisc)
  
  p <- p[, lapply(.SD, mean), by = "Replicates"]
  p <- p[, -1]
  p <- t(p)

  pg <- rbind(p, g)
  
  hc <- hclust(dist(pg))
  plot(hc)
  
  #library(cluster)
  set.seed(125)
  gap <- clusGap(pg, FUN = kmeans, iter.max = 30, K.max = 150, B = 50, verbose=interactive())
  plot(gap, main = "Gap Statistic")
  with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
  
  plot(hc)
  rect.hclust(hc, k = 11, border = "red")
  
  tr <- as.data.frame(cutree(hc, k = 11))
  pg[which(cutree(hc, k = 11) == 7),]
  # So far, the most promising.
  
  # Looking at overlap between two clustering methods.
  rownames(g[which(cutree(hcg, k = 10) == 1),]) %in% rownames(pg[which(cutree(hc, k = 11) == 5),])
  # Fair amount of overlap.
  rownames(g[which(cutree(hcg, k = 10) == 2),]) %in% rownames(pg[which(cutree(hc, k = 11) == 6),])
  
  # Use clustering of pheno/RNA for network. 
  write.csv(tr, "pg.csv")
  write.csv(trg, "g.csv")
  
  
  
  
## More clustering techniques for time series data.  
  RNA <- read.csv(file = "largeDE.csv", row.names = 1)
  
  RNA <- as.data.frame(t(RNA))
  
  discRNA <- discretize(RNA, method = "quantile", breaks = 3)
  
  for (i in 1:dim(discRNA)[2]){
    levels(discRNA[, i]) <- c(-1, 0, 1)
    discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
  }
  
  discRNA <- t(discRNA)
  f <- t(discRNA)
  
  Replicates <- rep(c(1:24), each = 2)
  
  f <- cbind(f, Replicates)
  f <- as.data.table(f)
  
  g <- f[, lapply(.SD, mean), by = "Replicates"]
  g <- g[, -1]
  g <- t(g)
  
  library(TSclust)
  
  d <- diss(g, METHOD = "FRECHET")
  
  
  
  hcg <- hclust(dist(g, method = "DTW"))
  plot(hcg)  
  
  