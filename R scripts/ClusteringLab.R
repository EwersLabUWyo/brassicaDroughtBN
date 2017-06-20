# http://jnmaloof.github.io/BIS180L_web/2017/05/25/genetic-networks-1/


# Read in data. 
# Read in differentially expressed genes. 
DE_genes <- read.csv(url("http://jnmaloof.github.io/BIS180L_web/data/DEgenes_GxE.csv"))
# Read in transformed RNA-Seq data. Transformation is to make data normally distributed. 
brass_voom_E <- read.csv(url("http://jnmaloof.github.io/BIS180L_web/data/voom_transform_brassica.csv"))

# Store DE gene names. 
DE_gene_names <- row.names(DE_genes)

# Subset DE genes from brass_voom_E data frame.
GxE_counts <- as.data.frame(brass_voom_E[DE_gene_names,])
# Convert to a matrix. 
GxE_counts <- as.matrix(GxE_counts) # some of the downstream steps require a data matrix

# Cluster by genes
hr <- hclust(dist(GxE_counts))
plot(hr)


hc <- hclust(dist(t(GxE_counts)))
plot(hc)

?rect.hclust
hc <- hclust(dist(t(GxE_counts)))
plot(hc) #redraw the tree everytime before adding the rectangles
rect.hclust(hc, k = 4, border = "red")

plot(hc) #redraw the tree everytime before adding the rectangles
rect.hclust(hc, k = 6, border = "red")

plot(hc) #redraw the tree everytime before adding the rectangles
rect.hclust(hc, h = 25, border = "red")

library(pvclust)
?pvclust #check out the documentation

set.seed(12456) #This ensure that we will have consistent results with one another
# Normally do 1,000+ bootstrap samples
fit <- pvclust(GxE_counts, method.hclust = "ward.D", method.dist = "euclidean", nboot = 50)
plot(fit) # dendogram with p-values

fit <- pvclust(GxE_counts, method.hclust = "ward.D", method.dist = "euclidean", nboot = 500)
plot(fit) # dendogram with p-values

library(gplots) #not to be confused with ggplot2!
heatmap.2(GxE_counts, Rowv = as.dendrogram(hr), scale = "row", density.info="none", trace="none")

library(ggplot2)
prcomp_counts <- prcomp(t(GxE_counts)) #gene wise
scores <- as.data.frame(prcomp_counts$rotation)[,c(1,2)]

set.seed(25) #make this repeatable as kmeans has random starting positions
fit <- kmeans(GxE_counts, 4)
clus <- as.data.frame(fit$cluster)
names(clus) <- paste("cluster")

plotting <- merge(clus, scores, by = "row.names")
plotting$cluster <- as.factor(plotting$cluster)

# plot of observations
ggplot(data = plotting, aes(x = PC1, y = PC2, label = Row.names, color = cluster)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(alpha = 0.8, size = 4, stat = "identity") 



library(cluster)
set.seed(125)
gap <- clusGap(GxE_counts, FUN = kmeans, iter.max = 30, K.max = 20, B = 50, verbose=interactive())
plot(gap, main = "Gap Statistic")
with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))

# http://genomicsclass.github.io/book/pages/pca_svd.html
# https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca



