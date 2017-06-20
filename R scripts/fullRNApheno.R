# fullRNApheno.R
# R version 3.3.1 (2016-06-21)
# March 15, 2017. Mallory B. Lai.
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

# Repeat rows to match number of observations in pheno dataframe.
discRNA <- do.call(rbind, replicate(6, discRNA, simplify = F))

# Order RNA dataframe by Timepoint and int. 
discRNA <- discRNA[with(discRNA, order(Timepoint, int)), ]

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

# Create blacklist for arcs to soil moisture. 
bl <- data.frame(From = colnames(rnaPheno)[- which(colnames(rnaPheno) == "SM")], 
                 To = rep("SM", (dim(rnaPheno)[2] - 1)))

# Create a whitelist using expert knowledge of 
# physiological interactions.  
wh <- data.frame(from = c("SM", "Photo", "Starch"), 
                 to = c("gs", "Fv.Fm.", "NSC"))

# Learn network structure. 
bn <- rsmax2(rnaPheno, restrict = "aracne", 
             blacklist = bl, whitelist = wh, 
             maximize = "tabu", score = "bde", 
             maximize.args = list(iss = 5))
#plot(bn)

# Write csv of network arcs. 
write.csv(bn$arcs, file = "fullLargeBNarcs.csv")    
