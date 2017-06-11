# clusteredRNApheno3.R
# R version 3.3.1 (2016-06-21)
# June 8, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating combined pheno & RNA seq BN for clustered Brassica gene data   
# using bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
#-----------------------------------------------------------------------

#### Preprocessing: 

# Set working directory.
# setwd("E:/Mallory Lai/PhenoRNAnetworkBrassica")
#setwd("/Users/mblai/Documents/GitHub/PhenoRNAnetworkBrassica")

######## RNA

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

# Remove rnaNames vector and RNA dataframe. 
rm(rnaNames)
rm(RNA)

# Order RNA dataframe by Timepoint and int. 
discRNA <- discRNA[with(discRNA, order(Timepoint, int)), ]

# Convert data intervals to -1, 0, and 1 representing low, medium, and high. 
for (i in 1:(dim(discRNA)[2] - 2)){
  levels(discRNA[, i]) <- c(-1, 0, 1)
  discRNA[, i] <- as.factor(discRNA[, i])
}


# Transform RNA data frame. 
discRNA <- t(discRNA)

# Rename column names to match timepoint and treatment. 
colnames(discRNA) <- paste(discRNA[dim(discRNA)[1] - 1,], 
                           discRNA[dim(discRNA)[1],], sep = "")

# Read in cluster classification. 
tr <- read.csv(file = "pg.csv")

# Rename column names to "variable" and "cluster."
colnames(tr) <- c("variable", "cluster")

# Remove phenotype variables. 
ph <- c("Photo", "gs", "Fv.Fm.", "Starch", "NSC", "SM")
tr <- tr[-which(tr$variable %in% ph), ]

# Create empty matrix for storing clusters.
cl <- matrix(data = NA, nrow = dim(tr)[1], ncol = (max(unique(tr$cluster))))
colnames(cl) <- paste("c", 1:(max(unique(tr$cluster))), sep = "")

# Separate by clusters.
for (i in 1:(max(unique(tr$cluster))))
{
  cl[, i] <- c(as.character(tr[tr$cluster == i, "variable"]), 
               rep(NA, length(cl[,i]) - dim(tr[tr$cluster == i, ])[1]))
}


# Separate clusters to form gene modules. 
c1 <- discRNA[c(cl[1:(length(cl[,1]) - sum(is.na(cl[,1]))), 1]), ]

# Add timepoint and treatment row. 
c1 <- rbind(c1, TP = discRNA[dim(discRNA)[1] - 1,], Trmt = discRNA[dim(discRNA)[1],])

# Transpose. 
c1 <- as.data.frame(t(c1))

# Gather all genes into one gene module, M1. 
m1 <- gather(c1, Gene, M1, -TP, -Trmt)

# Count number of values per treatment and time point. 
mCounts <- as.data.frame(table(m1[, c(1, 2, 4)]))

# Order mCounts dataframe by Timepoint, treatment, and module value. 
mCounts <- mCounts[with(mCounts, order(TP, Trmt, M1)), ]

# Add column for proportion per timepoint, treatment, and module value. 
mCounts$Prop <- mCounts$Freq/((length(cl[,1]) - sum(is.na(cl[,1])))*2)

# Create a column with the number of counts proportional to 12.
mCounts$twelve <- mCounts$Prop * 12

# Round the counts proportional to 12. 
mCounts$round <- round(mCounts$twelve)

# Convert to data.table. 
mCounts <- as.data.table(mCounts)

# Add a column that counts the total in round.  
mCounts <- mCounts[, total := sum(round), by = list(TP, Trmt)]

# Convert to dataframe for easier updates. 
mCounts <- as.data.frame(mCounts)

# If total count is 13, round to nearest 1/2 number and then round that. 
mCounts[mCounts$total == 13, 'round'] <- round(round(mCounts[mCounts$total == 13, 'twelve']/0.5) *.5)

# Convert back to data table for easy updating.   
mCounts <- as.data.table(mCounts)

# Update the total column with new rounding. 
mCounts <- mCounts[, total := sum(round), by = list(TP, Trmt)]

# Add a column that identifies the max proportion per treatment and time point. 
mCounts <- mCounts[, max := max(Prop), by = list(TP, Trmt)]

# Convert back to dataframe for easy subsetting. 
mCounts <- as.data.frame(mCounts)

# Add 1 to the count with the max proportion. 
mCounts[mCounts$total == 11 & mCounts$Prop == mCounts$max, 'round'] <- 
       mCounts[mCounts$total == 11 
          & mCounts$Prop == mCounts$max, 'round'] + 1

# Convert back to data table for easy updating.   
mCounts <- as.data.table(mCounts)

# Update the total column with new rounding. 
mCounts <- mCounts[, total := sum(round), by = list(TP, Trmt)]

# Extract the round column as a dataframe. 
mod <- data.frame(TP = mCounts$TP, Trmt = mCounts$Trmt, 
                  Value = mCounts$M1, M1 = mCounts$round)

# Remove unnecessary dataframes. 
rm(m1)
rm(c1)
rm(mCounts)

############ TO DO: Update loop to perform a similar function as above. 


# Loop through remaining clusters. 
for (i in 3:(dim(cl)[2]))
{
  # Separate clusters to form gene modules. 
  c1 <- discRNA[c(cl[1:(length(cl[,i]) - sum(is.na(cl[,i]))), i]), ]
  
  # Add timepoint and treatment row. 
  c1 <- rbind(c1, TP = discRNA[dim(discRNA)[1] - 1,], Trmt = discRNA[dim(discRNA)[1],])
  
  # Transpose. 
  c1 <- as.data.frame(t(c1))
  
  # Gather all genes into one gene module 
  module <- gather(c1, Gene, M, -TP, -Trmt)
  
  # Order module dataframe by Timepoint and int. 
  module <- module[with(module, order(TP, Trmt)), ]
  
  # Count number of values per treatment and time point. 
  counts <- as.data.frame(table(module[, c(1, 2, 4)]))
  
  # Order mCounts dataframe by Timepoint, treatment, and module value. 
  counts <- counts[with(counts, order(TP, Trmt, M)), ]
  
  # Add column for proportion per timepoint, treatment, and module value. 
  counts$Prop <- counts$Freq/((length(cl[,i]) - sum(is.na(cl[,i])))*2)
  
  # Create a column with the number of counts proportional to 12.
  counts$twelve <- counts$Prop * 12
  
  # Round the counts proportional to 12. 
  counts$round <- round(counts$twelve)
  
  # Convert to data.table. 
  counts <- as.data.table(counts)
  
  # Add a column that counts the total in round.  
  counts <- counts[, total := sum(round), by = list(TP, Trmt)]
  
  # Convert to dataframe for easier updates. 
  counts <- as.data.frame(counts)
  
  # If total count is 13, round to nearest 1/2 number and then round that. 
  counts[counts$total == 13, 'round'] <- round(round(counts[counts$total == 13, 'twelve']/0.5) *.5)
  
  # Convert back to data table for easy updating.   
  counts <- as.data.table(counts)
  
  # Update the total column with new rounding. 
  counts <- counts[, total := sum(round), by = list(TP, Trmt)]
  
  # Add a column that identifies the max proportion per treatment and time point. 
  counts <- counts[, max := max(Prop), by = list(TP, Trmt)]
  
  # Add a column that identifies the second highest proportion per treatment 
  # and time point. 
  counts <- counts[, max2 := as.numeric(Prop)][, max2 := sort(Prop, T)[2], by = list(TP, Trmt)]
  
  # Add a column that identifies the min proportion per treatment and time point. 
  counts <- counts[, min := min(Prop), by = list(TP, Trmt)]
  
  # Convert back to dataframe for easy subsetting. 
  counts <- as.data.frame(counts)
  
  # Add 1 to the count with the max proportion. 
  counts[counts$total <= 11 & counts$Prop == counts$max, 'round'] <- 
    counts[counts$total <= 11 
            & counts$Prop == counts$max, 'round'] + 1
  
  # Convert back to data table for easy updating.   
  counts <- as.data.table(counts)
  
  # Update the total column with new rounding. 
  counts <- counts[, total := sum(round), by = list(TP, Trmt)]
  
  # Convert back to dataframe for easy subsetting. 
  counts <- as.data.frame(counts)
  
  # Add 1 to the count with the 2nd highest proportion if 
  # still less than 12. 
  counts[counts$total <= 11 & counts$Prop == counts$max2, 'round'] <- 
    counts[counts$total <= 11 
           & counts$Prop == counts$max2, 'round'] + 1
  
  # Convert back to data table for easy updating.   
  counts <- as.data.table(counts)
  
  # Update the total column with new rounding. 
  counts <- counts[, total := sum(round), by = list(TP, Trmt)]
  
  # If there are any column totals of 13, subtract one from the
  # value with the lowest proportion. 
  counts[counts$total == 13 & counts$Prop == counts$min, 'round'] <- 
    counts[counts$total == 13 
           & counts$Prop == counts$min, 'round'] - 1
  
  # Convert back to data table for easy updating.   
  counts <- as.data.table(counts)
  
  # Update the total column with new rounding. 
  counts <- counts[, total := sum(round), by = list(TP, Trmt)]
  
  # Bind round column to mod dataframe. 
  mod <- cbind(mod, counts$round)

}

# Remove unneccesary dataframes. 
rm(tr)
rm(c1)
rm(i)
rm(ph)
rm(cl)
rm(discRNA)
rm(module)

# Rename modules in mod dataframe. 
colnames(mod)[5:dim(mod)[2]] <- paste("M", 2:10, sep = "")

# Expand modules to match the number of -1, 0, and 1's that should be 
# in the RNA-Seq dataframe. 
RNAmod <- as.data.frame(sapply(mod[, 4:dim(mod)[2]], function(x) rep(mod$Value, times = x)))

#### Pheno.

# Read in phenotype file. 
Pheno <- read.csv(file = "PhenoBrassicaImp.csv", row.names = 1)

# Rename SM... to get rid of periods. 
colnames(Pheno)[8] <- "SM"

# Add a column for Time of Day, named TOD.
Pheno$TOD <- rep(c(7, 11, 15, 19, 23, 3), each = 24, 2)

#### Discretize data. 

# Discretize the phenotype data, excluding fluorescence. 
phenoDisc <- discretize(Pheno[, c(3, 4, 6, 7, 8)], 
                        method = "interval", 
                        breaks = c(5, 5, 5, 5, 5))

# Use arules package to discretize fluorescence and detach due
# to the overlap in package functions with bnlearn. 
library(arules)
fluor <- discretize(Pheno[, 5], method = "cluster")
detach("package:arules", unload=TRUE)

# Attach fluorescence data to phenoDisc dataframe. 
phenoDisc$fluor <- fluor

# Add INT column to discretized data. 
phenoDisc$INT <- as.factor(Pheno$Treatment)

# Add Timepoint column to discretized data. 
phenoDisc$TP <- as.factor(Pheno$Timepoint)

# Order Pheno dataframe by Timepoint and int. 
phenoDisc <- phenoDisc[with(phenoDisc, order(TP, INT)), ]

# Combine pheno and RNA data frames.
rnaPheno <- cbind(phenoDisc, RNAmod)

# Remove unneccesary dataframes. 
rm(phenoDisc)
rm(Pheno)

# Remove extra columns.
rnaPheno$Trmt <- NULL
rnaPheno$INT <- NULL
rnaPheno$TP <- NULL
rnaPheno$TP <- NULL

# Subset training data and set aside test data. 
# To subset 10% of the data, we need to randomly select 30 samples
# to use as test data. 
set.seed(3)
testData <- sample(1:dim(rnaPheno)[1], round(dim(rnaPheno)[1]*.1), replace = F)

# Subset test data. 
test <- rnaPheno[testData, ]

# Isolate test data from training data. 
training <- rnaPheno[-testData, ]



# Learn network structure. 
bn <- suppressWarnings(tabu(training, score = "bde", 
                            iss = 10, tabu = 50))

plot(bn)




boot <- boot.strength(training, R = 500, algorithm = "tabu", 
                      algorithm.args = list(score = "bde", iss = 10))

boot[(boot$strength > 0.85) & (boot$direction >= 0.5), ]
avg.boot <- averaged.network(boot, threshold = 0.85)
plot(avg.boot)


nodes <- names(training)
start <- random.graph(nodes = nodes, method = "ic-dag", num = 100, 
                      every = 3)
netlist <- suppressWarnings(lapply(start, function(net){
  tabu(training, score = "bde", tabu = 50, iss = 10)
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)

plot(bn)







# Write csv of network arcs. 
#write.csv(bn$arcs, file = "clustBNarcs.csv")   
plot(bn)

bn.mle <- bn.fit(bn, training, method = "bayes")

bn.fit.barchart(bn.mle$Photo, xlab = "P()")

bn.mle$Photo

bn.fit.barchart(bn.mle$M8, xlab = "P()")

bn.mle$M8

bn.fit.barchart(bn.mle$M3, xlab = "P()")
bn.mle$M3

bn.fit.barchart(bn.mle$Starch, xlab = "P()")
bn.mle$Starch

bn.fit.barchart(bn.mle$Starch, xlab = "P()")
bn.mle$M4
testOrder <- test[with(test, order(M9, M4)), ]
testOrder <- testOrder[, c(12, 22)]
s <- testOrder[testOrder$M4 == "-1" & testOrder$M9 == "-1", ]
l <- testOrder[testOrder$M4 == "-1" & testOrder$M9 == "0", ]
c <- testOrder[testOrder$M4 == "-1" & testOrder$M9 == "1", ]

