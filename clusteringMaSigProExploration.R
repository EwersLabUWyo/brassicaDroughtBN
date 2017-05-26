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
    BrassicaFPKM <- BrassicaFPKM[rowSums(BrassicaFPKM > 10) >= 24,]
    
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
    formatExprDesign <- make.design.matrix(ExprDesign, degree = 9)
    
    # Find differentially expressed genes.   
    # Perform regression fit for time series. 
    fit <- p.vector(ExprMatrix, design = formatExprDesign, 
                    Q = 0.01, counts = F)
      
      # Select regression model by stepwise regression. 
      model <- T.fit(fit, step.method = "two.ways.backward", alfa = .01)
      
      # Investigation of influential data was performed. 
      
      # Get the significantly expressed genes. 
      sigs <- get.siggenes(model, vars="groups", rsq = 0)
      
      largeDE <- ExprMatrix[rownames(sigs$sig.genes$WvsD$sig.profiles), ]
      write.csv(largeDE, "largeDE.csv")
      
      clustgenes <- see.genes(ExprMatrix, ExprDesign, k = 100,
                              cluster.method = "hclust", distance = "cor",
                              agglo.method = "ward.D")  
      
      # Note: distance = "cor" performs considerably worse than "euclidean"
      
      
      # Will try discretizing then clustering. 
      RNA <- as.data.frame(t(ExprMatrix))
      discRNA <- discretize(RNA, method = "quantile", breaks = 3)
      
      
      for (i in 1:7791){
        levels(discRNA[, i]) <- c(-1, 0, 1)
        discRNA[, i] <- as.numeric(as.character(discRNA[, i]))
      }
      
      discRNA <- t(discRNA)
      colnames(discRNA) <- colnames(ExprMatrix)
      
      clustgenes <- see.genes(discRNA, ExprDesign, k = 100,
                              cluster.method = "hclust", distance = "cor",
                              agglo.method = "ward.D")  # This does a terrible job. 
      
      clustgenes <- see.genes(discRNA, ExprDesign, k = 100,
                              cluster.method = "hclust", distance = "euclidean",
                              agglo.method = "ward.D")  
      
      
      
      # Sort D and W.
      d <- sort(paste(discRNA[1, 1:24], sep = ", "))
     
      
      
      
      cgroups <- as.data.frame(clustgenes$cut)
      
      clust1 <- rownames(cgroups)[which(cgroups == 1)]
      
      c1g <- ExprMatrix[clust1, ]

      
      
      
      
      
      
      
      
      
      
            
#### Clustering
    
    # Cluster by genes
    hr <- hclust(dist(de))
    plot(hr)
    rect.hclust(hr, k = 15, border = "red")
    
    gap <- clusGap(dist(de), FUN = kmeans, iter.max = 30, K.max = 20, B = 50, verbose=interactive())
    plot(gap)
    
    # Cluster by sample. 
    hr <- hclust(dist(t(de)))
    plot(hr, xlab = "allTPde")

    gap <- clusGap(t(de), FUN = kmeans, iter.max = 30, K.max = 20, B = 50, verbose=interactive())
    plot(gap)
       
    rect.hclust(hr, k = 6, border = "red")
     
    # Cluster by sample with subset of non-de genes. 
        # Log2 transform wwBrassicaFPKM data set. 
        log2FPKM <- log2(BrassicaFPKM)
        
        # Scale the dataset.
        log2FPKM <- scale(log2FPKM)
        
        # Susbet first 100 genes. 
        genes <- log2FPKM[1:1000, ]

        # Cluster by sample. 
        hr <- hclust(dist(t(genes)))
        plot(hr, xlab = "nonde")
        rect.hclust(hr, k = 6, border = "red")

        # Cluster by genes. 
        # Cluster by genes
        hr <- hclust(dist(genes))
        plot(hr)
        
        gap <- clusGap(genes, FUN = kmeans, iter.max = 30, K.max = 20, B = 50, verbose=interactive())
        plot(gap)
        
        plot(hr)
        rect.hclust(hr, k = 4, border = "red")

degenes <- read.csv(file.choose(), row.names = 1)        
# Cluster by sample. 
hr <- hclust(dist(t(degenes)))
plot(hr, xlab = "de")
rect.hclust(hr, k = 6, border = "red")

d <- scale(degenes)
hr <- hclust(dist(t(d)))
plot(hr, xlab = "scaledDE")
rect.hclust(hr, k = 6, border = "red")


de <- read.csv(file.choose(), row.names = 1)        

