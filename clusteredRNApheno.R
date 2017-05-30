# clusteredRNApheno.R
# R version 3.3.1 (2016-06-21)
# May 28, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating combined pheno & RNA seq BN for Brassica data   
# using bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
#-----------------------------------------------------------------------

#### Preprocessing: 

# Set working directory.
# setwd("E:/Mallory Lai/PhenoRNAnetworkBrassica")
#setwd("/Users/mblai/Documents/GitHub/PhenoRNAnetworkBrassica")

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
phenoDisc <- phenoDisc[with(phenoDisc, order(TP, INT)), ]

# Remove unnecesary dataframes.
rm(Pheno)







# Read in RNAseq data. 
RNA <- read.csv(file = "largeDE.csv", row.names = 1)

# Read in cluster classification. 
tr <- read.csv(file = "pg.csv")

# Rename column names to "variable" and "cluster."
colnames(tr) <- c("variable", "cluster")

# Remove phenotype variables. 
ph <- c("Photo", "gs", "Fv.Fm.", "Starch", "NSC", "SM")
tr <- tr[-which(tr$variable %in% ph), ]

# Create empty matrix for storing clusters.
cl <- matrix(data = NA, nrow = dim(tr)[1], ncol = (length(unique(tr$cluster))))
colnames(cl) <- paste("c", 1:(length(unique(tr$cluster))), sep = "")

# Separate by clusters.
for (i in 1:(length(unique(tr$cluster))))
{
  cl[, i] <- c(as.character(tr[tr$cluster == i, "variable"]), 
               rep(NA, length(cl[,i]) - dim(tr[tr$cluster == i, ])[1]))
}

# Look at dim of each cluster
for (i in 1:(length(unique(tr$cluster))))
{
  print(length(cl[,i]) - sum(is.na(cl[,i])))
}

# Note: 488 is the maximum cluster size. (x2 for replicates is 976)


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

# Convert data intervals to -1, 0, and 1 representing low, medium, and high. 
for (i in 1:(dim(discRNA)[2] - 2)){
  levels(discRNA[, i]) <- c(-1, 0, 1)
  discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
}

# Transform RNA data frame. 
discRNA <- t(discRNA)

# Rename column names to match timepoint and treatment. 
colnames(discRNA) <- paste(discRNA[dim(discRNA)[1] - 1,], 
                           discRNA[dim(discRNA)[1],], sep = "")

# Separate clusters. 
c1 <- discRNA[c(cl[1:(length(cl[,1]) - sum(is.na(cl[,1]))), 1]), ]









# Remove unneccesary dataframes. 
rm(RNA)
rm(rnaNames)

# Combine pheno and RNA data frames.
rnaPheno <- cbind(phenoDisc, discRNA)

# Remove unneccesary dataframes. 
rm(discRNA)
rm(phenoDisc)

# Remove extra columns.
rnaPheno$Timepoint <- NULL
rnaPheno$int <- NULL
rnaPheno$INT <- NULL
rnaPheno$TP <- NULL
rnaPheno$TOD <- NULL

# Subset training data and set aside test data. 
# To subset 10% of the data, we need to randomly select 30 samples
# to use as test data. 
set.seed(3)
testData <- sample(1:288, 30, replace = F)

# Subset test data. 
test <- rnaPheno[testData, ]

# Isolate test data from training data. 
training <- rnaPheno[-testData, ]

# Create blacklist for arcs to soil moisture. 
bl <- data.frame(From = colnames(rnaPheno)[- which(colnames(rnaPheno) == "SM")], 
                 To = rep("SM", (dim(rnaPheno)[2] - 1)))

# Learn network structure. 
bn <- suppressWarnings(tabu(training, score = "bde", blacklist = bl,
                            iss = 10, tabu = 50))

# Write csv of network arcs. 
write.csv(bn$arcs, file = "blLargeBNarcs.csv")    
