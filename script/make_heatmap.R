#### Packages ####
library(gplots)
library(RColorBrewer)
library(dplyr)
library(glue)
library(dendextend)

#### Functions ####
remove_na_dist <- function(mat) {
  dist.mat <- as.matrix(dist(mat))
  
  distNAs <- sum(is.na(dist.mat))
  
  if (distNAs > 0) {
    giveNAs <- which(is.na(dist.mat),
                     arr.ind=TRUE)
    
    tab <- sort(table(c(giveNAs)), 
                decreasing=TRUE)
    
    checkNA <- sapply(1:length(tab), function(i) {
      sum(is.na(as.matrix(dist(mat[-as.numeric(names(tab[1:i])), ]))))})
    
    rmv <- names(tab)[1:min(which(checkNA == 0))]
    
    return(mat[-as.numeric(rmv), ])
    
  } else {
    return(mat)
  }
}

#### Parameters ####
alpha <- 0.05
cluster.cols <- FALSE
k.cluster.rows <- FALSE

gse <- 'PROJ1624_1'

heatmap.filename <- glue("~/mrc/project/sita_slam-seq/processed/heatmap_{gse}.png")

delete.columns <- c()

#### Load data ####
fc.table <- read.csv(glue("~/mrc/project/sita_slam-seq/data/{gse}_foldChange.txt"),
                          sep="\t")
padj.table <- read.csv(glue("~/mrc/project/sita_slam-seq/data/{gse}_padj.txt"),
                       sep="\t")

# Delete columns
if (length(delete.columns) > 0) {
  fc.table <- fc.table[, -which(names(fc.table) %in% delete.columns)]
  padj.table <- padj.table[, -which(names(padj.table) %in% delete.columns)]
}

# Convert to matrix
fc.matrix <- data.matrix(fc.table[, 2:ncol(fc.table)])

# Set gene symbol as matrix row names
rownames(fc.matrix) <- fc.table[, 1]

#### Filtering ####
# Convert non-significant gene fold-changes with p < alpha to NA
fc.matrix[which(padj.table[, c(2:ncol(padj.table))] > alpha)] <- 0
fc.matrix[is.na(padj.table[, c(2:ncol(padj.table))])] <- 0
fc.matrix <- fc.matrix[rowSums(abs(fc.matrix) < 5) == ncol(fc.matrix), ]

# Remove rows and columns with all NA or 0
#fc.matrix <- fc.matrix[rowSums(is.na(fc.matrix)) != ncol(fc.matrix),
#                       colSums(is.na(fc.matrix)) != nrow(fc.matrix)]

#fc.matrix <- fc.matrix[rowSums(fc.matrix) != 0, ]

# Remove rows giving NA in distance matrix for clustering
# fc.filt.matrix <- remove_na_dist(fc.matrix)
# print(glue("Genes remaining: {nrow(fc.filt.matrix)} out of {nrow(fc.table)}"))

#### Heatmap ####
# Set up colour palette
palette <- colorRampPalette(c("red", "black", "green"))(n=299)

# Clustering on rows and columns
k.rows <- 3
k.cols <- 2

Rowv <- fc.matrix %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=k.rows) %>% 
  set("labels_col", k=k.rows) %>%
  set("branches_lwd", 4) %>%
  rotate_DendSer(ser_weight=dist(fc.matrix))

Colv <- fc.matrix %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=k.cols) %>% 
  set("branches_lwd", 4) %>%
  rotate_DendSer(ser_weight=dist(t(fc.matrix)))

# Plot and save static heatmap
png(file=heatmap.filename, 
    width=3000, height=9000, res=300)
hm <- heatmap.2(fc.matrix,
                scale="none",
                main=glue("Fold-changes of selected genes (n={nrow(fc.matrix)})
                in PROJ1624 comparisons (N={ncol(fc.matrix)}),
                p < {alpha}"),
                Rowv=if (k.cluster.rows) Rowv else TRUE,
                Colv=cluster.cols,
                dendrogram=if (cluster.cols) "both" else "row",
                srtCol=45, offsetCol=0.5,
                cexRow=1.5, cexCol=2,
                margins=c(12, 14),
                labRow=row.names(fc.matrix),
                labCol=gsub("_vs_", " vs. \n", colnames(fc.matrix)),
                notecol="black",
                density.info="none",
                trace="none",
                col=palette,
                # breaks=seq(-5, 5, length.out=300), 
                symbreaks=TRUE,
                na.color="black",
                key=TRUE, keysize=1, lhei=c(1, 15), key.xlab="log2(foldChange)")
hm
dev.off()
