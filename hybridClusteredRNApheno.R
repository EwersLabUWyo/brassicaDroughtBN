# hybridClusteredRNApheno.R
# R version 3.3.1 (2016-06-21)
# May 28, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating combined pheno & RNA seq BN for clustered Brassica gene data   
# using bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
library(dplyr)
library(tidyr)
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

# Write csv of cluster dataframe. 
#cl <- cl[, -2]

#write.csv(cl, "clusters.csv")

cmax <- 0

# Look at dim of each cluster
for (i in 1:(max(unique(tr$cluster))))
{
  if (length(cl[,i]) - sum(is.na(cl[,i])) > cmax){
    cmax <- length(cl[,i]) - sum(is.na(cl[,i]))
  }
}

# Note: 488 is the maximum cluster size. (x2 for replicates is 976)
cmax <- cmax * 2 

# Separate clusters to form gene modules. 
c1 <- discRNA[c(cl[1:(length(cl[,1]) - sum(is.na(cl[,1]))), 1]), ]

# Add timepoint and treatment row. 
c1 <- rbind(c1, TP = discRNA[dim(discRNA)[1] - 1,], Trmt = discRNA[dim(discRNA)[1],])

# Transpose. 
c1 <- as.data.frame(t(c1))

# Gather all genes into one gene module, M1. 
m1 <- gather(c1, Gene, M1, -TP, -Trmt)

# Order m1 dataframe by Timepoint and int. 
m1 <- m1[with(m1, order(TP, Trmt)), ]

# Group by time point and treatment. 
m1 <- m1 %>% group_by(TP, Trmt)

# Randomly sample maximum module size * 2 (two replicates; cmax) 
# for each time point and treatment.
m1 <- sample_n(m1, cmax, replace = T)

# Convert to data.frame. 
m1 <- as.data.frame(m1)

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
  
  # Group by time point and treatment. 
  module <- module %>% group_by(TP, Trmt)
  
  # Randomly sample maximum module size * 2 (two replicates; cmax) 
  # for each time point and treatment.
  module <- sample_n(module, cmax, replace = T)
  
  # Rename M to match module number. 
  colnames(module)[which(colnames(module) == "M")] <- paste("M", i, sep = "")
  
  # Convert to data.frame. 
  module <- as.data.frame(module)
  
  # Combine into full dataframe. 
  m1 <- cbind(m1, module[c("Gene", paste("M", i, sep = ""))])
}

# Remove unneccesary dataframes. 
rm(tr)

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
                        method = "hartemink")

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

# Remove unnecesary dataframes.
rm(Pheno)

# To get phenotype data to be of the same length as the RNA data, 
# group by time point and treatment and randomly sample from each section. 


# Group by time point and treatment. 
Pheno <- phenoDisc %>% group_by(TP, INT)

# Randomly sample maximum module size * 2 (two replicates; cmax) 
# for each time point and treatment.
Pheno <- sample_n(Pheno, cmax, replace = T)

# Convert to data.frame. 
Pheno <- as.data.frame(Pheno)

# Combine pheno and RNA data frames.
rnaPheno <- cbind(Pheno, m1)

# Remove unneccesary dataframes. 
rm(discRNA)
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

# Remove Gene column from rnaPheno. 
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL
training$Gene <- NULL

rm(m1)
rm(module)
rm(cl)
rm(c1)
rm(rnaPheno)
rm(ph)
rm(cmax)
rm(testData)

# Convert modules to factors. 
training[, 7:dim(training)[2]] <- lapply(training[, 7:dim(training)[2]], factor)

# Create a whitelist using expert knowledge of physiological interactions. 
wh <- data.frame(from = c("SM", "gs", "Photo"), to = c("gs", "Photo", "fluor"))

# Create a blacklist to soil moisture. 
bl <- tiers2blacklist(list(colnames(training)[5], colnames(training)[-5]))

# Model averaging hybrid
nodes <- names(training)
start <- random.graph(nodes = nodes, method = "ic-dag", num = 500, 
                      every = 50)
netlist <- suppressWarnings(lapply(start, function(net){
  rsmax2(training, whitelist = wh, restrict = "aracne", 
         maximize = "tabu", score = "bde", 
         maximize.args = list(iss = 15))
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)



bn <- rsmax2(training, whitelist = wh, blacklist = bl,
             restrict = "gs", 
             maximize = "tabu", score = "bde", 
             maximize.args = list(iss = 15))
plot(bn)



netlist <- suppressWarnings(lapply(start, function(net){
  aracne(training)
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)



netlist <- suppressWarnings(lapply(start, function(net){
  rsmax2(training, restrict = "si.hiton.pc", test = "x2", maximize = "tabu", score = "bde", 
         maximize.args = list(iss = 5))
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)




boot <- boot.strength(training, R = 250, algorithm = "mmhc")

boot[(boot$strength > 0.85) & (boot$direction >= 0.5), ]
avg.boot <- averaged.network(boot, threshold = 0.85)
plot(avg.boot)

#######################

cpquery(bn.bayes,
        event = (gs == "[0.00895,0.142]"),
        evidence = (M5 == "1") &
          (M6 == "1"))
cp <- cpdist(bn.bayes,
             nodes = "gs",
             evidence = (M5 == "1") &
               (M6 == "1") & (Photo == "(10,14.7]"))
cptab <- as.data.frame(prop.table(table(cp)))




an.mle <- bn.fit(an, training, method = "bayes")

bn.fit.barchart(an.mle$gs, xlab = "P()")



# pheno

ph <- phenoDisc[, 1:6]
nodes <- names(ph)
start <- random.graph(nodes = nodes, method = "melancon", num = 100, 
                      every = 3)
netlist <- suppressWarnings(lapply(start, function(net){
  tabu(training, score = "bde", tabu = 100, iss = 50)
}))



rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)
