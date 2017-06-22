# clusteringMaSigProExploration.R
# R version 3.3.1 (2016-06-21)
# March 11, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread

#-----------------------------------------------------------------------
library(maSigPro)
library(bnlearn)
library(stringr)
#-----------------------------------------------------------------------

#### Pre-processing: 

    # Read in the data. 
    dBrassicaFPKM <- read.table(file = "/Users/mblai/Documents/RNAseqData/R500_D_FPKM.txt", header = T)
    wwBrassicaFPKM <- read.table(file = "/Users/mblai/Documents/RNAseqData/R500_WW_FPKM.txt", header = T)
    
    #setwd("E:/Mallory Lai")
    
    #dBrassicaFPKM <- read.table(file = "R500_D_FPKM.txt", header = T)
    #wwBrassicaFPKM <- read.table(file = "R500_WW_FPKM.txt", header = T)
    
    # Remove isoforms. 
    dBrassicaFPKM <- dBrassicaFPKM[duplicated(dBrassicaFPKM$gene_short_name) == 0, ]
    wwBrassicaFPKM <- wwBrassicaFPKM[duplicated(wwBrassicaFPKM$gene_short_name) == 0, ]
    
    # Merge data.frames.
    BrassicaFPKM <- merge(dBrassicaFPKM, wwBrassicaFPKM, by = "gene_short_name")
    
    # Remove unnecessary data.frames. 
    rm(dBrassicaFPKM)
    rm(wwBrassicaFPKM)
    
    # Assign gene names to rownames. 
    rownames(BrassicaFPKM) <- BrassicaFPKM$gene_short_name
    
    # Remove columns without expression values. 
    BrassicaFPKM <- BrassicaFPKM[, -c(1, 2, 27)]
    
    # Reduce data set to only genes with some expression. 
    BrassicaFPKM <- BrassicaFPKM[rowSums(BrassicaFPKM > 0) >= 24,]
    
    # Log2 transform wwBrassicaFPKM data set. 
    log2FPKM <- log2(BrassicaFPKM + .1)
    
    # Scale the dataset.
    log2FPKM <- scale(log2FPKM)
    
    # Format for maSigPro.   
    # Coerce log2FPKM to a matrix. 
    ExprMatrix <- as.matrix(log2FPKM)
    
    # Create design dataframe. 
    Time <- rep(c(rep(c(1:12), each = 2)), 2)
    Replicates <- rep(c(1:24), each = 2)
    D <- rep(c(1,0), each = 24)
    W <- rep(c(0,1), each = 24)
    ExprDesign <- cbind(Time, Replicates, D, W)
    
    # Make rownames the column names of the expression matrix. 
    rownames(ExprDesign) <- colnames(ExprMatrix)
    
    # Remove design dataframe vectors. 
    rm(D)
    rm(Replicates)
    rm(Time)
    rm(W)
    
    # Remove log2FPKM dataframe.
    rm(log2FPKM)

#### Find differentially expressed genes using maSigPro    
    
    # Format the ExprDesign matrix. 
    formatExprDesign <- make.design.matrix(ExprDesign, degree = 11)
    
    # Find differentially expressed genes.   
    # Perform regression fit for time series. 
    fit <- p.vector(ExprMatrix, design = formatExprDesign, 
                    Q = 0.01, counts = F)
      
      # Select regression model by stepwise regression. 
      model <- T.fit(fit, step.method = "backward", alfa = .01)
      
      # Investigation of influential data was performed. 
      
      # Get the significantly expressed genes. 
      sigs <- get.siggenes(model, vars="groups", rsq = 0)
      
      largeDE <- ExprMatrix[rownames(sigs$sig.genes$WvsD$sig.profiles), ]
      
      write.csv(largeDE, "two.ways.forward.05.DE.csv")
      
      DE05FPKM <- BrassicaFPKM[rownames(sigs$sig.genes$WvsD$sig.profiles), ]
      write.csv(DE05FPKM, "DE05FPKM.csv")
      
      ############## Plot
      
      # Function to plot a gene of interest.
      plotGene <- function(i, k){
        PlotGroups(ExprMatrix[rownames(ExprMatrix) == i, ],
                   edesign = ExprDesign, main = i)
      }
      
    # Plot all DEgenes. 
      tiff(filename = "IntGenesDE05%03d.png", width = 8, height = 11,
           units = 'in', res = 300)
      par(mfrow = c(5,3))
      lapply(x, plotGene)
      
      
      dev.off()
      
      DE01FPKM <- BrassicaFPKM[rownames(de01), ]
      write.csv(DE01FPKM, "DE01FPKM.csv")
      DE05FPKM <- BrassicaFPKM[rownames(de05), ]
      write.csv(DE05FPKM, "DE05FPKM.csv")
      
      
      DE01FPKM <- BrassicaFPKM[rownames(de01), ]
      write.csv(DE01FPKM, "DE01FPKM.csv")
      DE05FPKM <- BrassicaFPKM[rownames(largeDE), ]
      write.csv(DE05FPKM, "DE05BackwardFPKM.csv")
      ####################################################

      # Function to plot a gene of interest.
      plotGene <- function(i, k){
        PlotGroups(ExprMatrix[rownames(ExprMatrix) == i, ],
                   edesign = ExprDesign, main = expgenes[i, 2])
      }
      
      
      # Subset genes of interest that are found to be differentially 
      # expressed. 
      intGenesDE05 <- expgenes[which(expgenes$V1 %in% rownames(de05)), ]
    
      # Plot all genes. 
      tiff(filename = "IntGenesDE05%03d.png", width = 8, height = 11,
           units = 'in', res = 300)
      par(mfrow = c(5,3))
      lapply(intGenesDE05$V1, plotGene)
      
      dev.off()
      
      # Subset genes of interest that are not found to be differentially 
      # expressed. 
      intGenesNonDE05 <- expgenes[-which(expgenes$V1 %in% rownames(de05)), ]
      
      # Plot all genes. 
      tiff(filename = "IntGenesNonDE05%03d.png", width = 8, height = 11,
           units = 'in', res = 300)
      par(mfrow = c(5,3))
      lapply(intGenesNonDE05$V1, plotGene)
      
      dev.off()
