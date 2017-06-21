# BHCdeClusteringExplorationFULL.R
# R version 3.3.1 (2016-06-21)
# June 13, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread

#-----------------------------------------------------------------------
library(BHC)
library(bnlearn)
library(data.table)
library(stringr)
library(stringi)
#-----------------------------------------------------------------------

#source("https://bioconductor.org/biocLite.R")
#biocLite("BHC")
#require(BHC)

# Read in RNA data. 
RNA <- read.csv(file = "BrassicaFPKM.csv", row.names = 1)

# Transpose to dataframe for discretize function. 
RNA <- as.data.frame(t(RNA))

# Discretize RNA into quantiles. 
discRNA <- discretize(RNA, method = "interval", breaks = 3)

# Convert RNA to -1, 0, or 1 for low, medium, or high expression.
for (i in 1:dim(discRNA)[2]){
  levels(discRNA[, i]) <- c(-1, 0, 1)
  discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
}

# Transpose data for BHC. 
discRNA <- t(discRNA)

# Remove RNA dataframe. 
rm(RNA)

# Define data dimensions. 
nDataItems <- nrow(discRNA)
nFeatures  <- ncol(discRNA)

# Define genes as item labels. 
itemLabels <- rownames(discRNA)

# Perform clustering with bhc function. 
hct <- bhc(discRNA, itemLabels, numReps=2)

# Write clusters to text file. 
WriteOutClusterLabels(hct, "labels.txt", verbose=TRUE)

# Read in the cluster output. 
clst <- readLines("labels.txt")

# Collapse into one string. 
clst <- paste(clst, sep = "", collapse = " ")

# Convert string to dataframe after splitting strings
# at "---".
clstFrame = as.data.frame(do.call(rbind, str_split(clst, "---")), 
                          stringsAsFactors=FALSE)

# Remove CLUSTER columns which are even numbered. 
clstFrame <- clstFrame[, -c(seq(2, 116, 2))]

# Remove first row, which is empty. 
clstFrame <- clstFrame[, -1]

# Trim whitespace from front and end of character string. 
clstFrame <- lapply(clstFrame, function(x){str_trim(x, side = "both")})

# Split each column string up into a list for each cluster. 
clList <- lapply(clstFrame, function(x) {unlist(str_split(x, " "))})

# Convert the list to a matrix and fill empty values with NA.
m <- stri_list2matrix(clList, fill = NA)

# Write matrix to csv file. 
write.csv(m, "FULLgenesMods.csv")




