# ClustGeneAnnotationBrapaDrought.R
# R version 3.3.1 (2016-06-21)
# March 16, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating annotations for clustered genes from Brassica data.
# Expression data taken from Brassica under well-watered
# and droughted conditions. 

#-----------------------------------------------------------------------
library(data.table)
library(stringr)
#-----------------------------------------------------------------------

# Read in gene annotation file. 
geneAnnot <- fread(file = "/Users/mblai/Documents/Brapa_197_annotation_info.txt", 
                   header = F)

# Convert to data.frame.
geneAnnot <- as.data.frame(geneAnnot)

# Assign Brassica gene as row names. 
rownames(geneAnnot) <- geneAnnot$V1

# Read in cluster genes. 
clusters <- read.csv(file = "clusters.csv",
                    stringsAsFactors = F, row.names = 1)

# Keep only gene names and cluster.
clust <- stack(clusters)

# Remove NA's. 
clust <- clust[!is.na(clust$values), ]

# Subset genes with commas.
commas <- data.frame(Gene = str_subset(clust$values, ","), 
                        Clust = clust[str_detect(clust$values, ",") == T, 2])

# Split genes by commas.
sepCommas <- unlist(str_split(commas$Gene, ","))

# Remove values with duplicates in clust or sepCommas as they could show up 
# in duplicate clusters. This will be fixed later. 

# Identify duplicate values from sepCommas. 
dup <- sepCommas[duplicated(sepCommas) == T]

# Remove from sepCommas list. 
sepCommas <- sepCommas[!sepCommas %in% dup]

# Remove genes from sepCommas list that are duplicated in Gene column 
# of cluster data.frame. 
sepCommas <- sepCommas[!sepCommas %in% clust$value]

# Remove genes with commas from clust list.
clust <- clust[setdiff(clust$values, commas$Gene), ]

# Create a dataframe for the genes with commas and their cluster id. 
# Initialize dataframe for storing variables. 
commaDF <- data.frame(values = rep(NA, length(sepCommas)), ind = rep(NA, length(sepCommas)))

for (i in 1:length(sepCommas)){
  commaDF[i, 1] <- sepCommas[i]
  commaDF[i, 2] <- as.character(clust[str_detect(clust$values, sepCommas[i]) == T, 2])
}

# Add genes split at comma to clust dataframe.
clust <- rbind(clust, commaDF)

# Match to annotation file. 
clustAnnot <- geneAnnot[clust$values, ]

clustAnnot$Cluster <- cbind(clust$ind)
