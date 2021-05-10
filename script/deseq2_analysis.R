#### Packages ####
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(ggrepel)
library(glue)

#### Functions ####
format_condition <- function (colnames) {
  new_cond <- gsub("_[0-9]*$", "", colnames)
  new_cond <- gsub("_rep[0-9]*$", "", new_cond)
  new_cond <- gsub("GSM[0-9]*_", "", new_cond)
  
  return(new_cond)
}

addEnsemblSymbol <- function (table) {
  genes <- row.names(table)
  
  if (grepl("ENSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "hsapiens_gene_ensembl"
    symbol <- "hgnc_symbol"
  } else if (grepl("ENSMUSG", genes[1], fixed=TRUE)) {
    ensemblDataset <- "mmusculus_gene_ensembl"
    symbol <- "mgi_symbol"
  }
  
  mart <- useDataset(ensemblDataset, useMart("ENSEMBL_MART_ENSEMBL", host="http://www.ensembl.org"))
  geneList <- getBM(filters="ensembl_gene_id",
                    attributes=c("ensembl_gene_id", symbol),
                    values=genes,
                    mart=mart)
  
  geneList <- distinct(geneList, ensembl_gene_id, .keep_all=TRUE)
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

deAnalysis <- function(counts.nascent, counts.total, treatment, control, a) {
  countData.nascent <- as.matrix(counts.nascent, rownames=1)
  countData.total <- as.matrix(counts.total, rownames=1)
  
  # Set contrast
  cond <- format_condition(colnames(countData.nascent))
  conditions <- unique(cond)
  
  # Create DESeq dataset object
  design <- data.frame(row.names=colnames(countData.nascent),
                       condition=factor(cond, levels=conditions))
  
  dds.nascent <- DESeqDataSetFromMatrix(countData=countData.nascent,
                                colData=design,
                                design= ~ condition)
  row.names(dds.nascent) <- row.names(countData.nascent)
  
  dds.total <- DESeqDataSetFromMatrix(countData=countData.total,
                                      colData=design,
                                      design= ~ condition)
  row.names(dds.total) <- row.names(countData.total)
  
  # Filter out uninformative rows
  dds.total <- dds.total[rowSums(counts(dds.nascent)) > 0, ]
  dds.nascent <- dds.nascent[rowSums(counts(dds.nascent)) > 0, ]
  
  dds.nascent$condition <- relevel(dds.nascent$condition,
                                   ref=control)
  dds.total$condition <- relevel(dds.total$condition,
                                 ref=control)
  
  # Run DESeq2 main command
  dds.total <- DESeq(dds.total)
  
  sizeFactors(dds.nascent) <- sizeFactors(dds.total) # apply size factors to tc counts
  dds.nascent <- DESeq(dds.nascent)
  
  # Get results
  res.nascent <- results(dds.nascent,
                         contrast=c("condition", treatment, control),
                         alpha=a)

  res.total <- results(dds.total,
                       contrast=c("condition", treatment, control),
                       alpha=a)
  
  # Put into list and return
  resList <- list()
  resList$nascent <- as.data.frame(res.nascent) %>% arrange(padj, desc(log2FoldChange))
  resList$total <- as.data.frame(res.total) %>% arrange(padj, desc(log2FoldChange))
  resList$conditions <- conditions
  
  return(resList)
}

PCAPlot <- function(counts, title) {
  colData <- data.frame(row.names=colnames(counts),
                        condition=format_condition(colnames(counts)))
  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData=colData,
                                design=~condition)
  
  # PCA plot
  rld <- vst(dds, blind=TRUE)
  
  plotPCA(rld, intgroup="condition", title)
}

MAPlot <- function(results, case, control, cutoff) {
  # Extract gene IDs for top 20 deregulated genes for plotting
  results$gene_name <- row.names(results)
  
  results <- results[!is.na(results$padj), ]
  
  dereg <- results[which(results$padj <= cutoff), ] # all significant deregulated genes
  nonsig <- results[!(results$gene_name %in% dereg$gene_name), ] # non-significant deregulated genes
  
  down <- nrow(dereg[dereg$log2FoldChange < 0, ]) # downregulated genes
  up <- nrow(dereg[dereg$log2FoldChange > 0, ]) # upregulated genes
  
  top.dereg <- dereg[order(abs(dereg$log2FoldChange), decreasing=TRUE), 1]
  top.dereg <- top.dereg[!is.na(top.dereg)]
  
  ## Generate summary of MA-plots with additional information & highlights
  # NB: Formatting of axes may shift marginal density plot from scatter plot
  # Only use exported plots for visual inspection
  
  # Generate basic MA-like plot with density coloring
  p <- results %>% ggplot(aes(x=baseMean, y=log2FoldChange)) +
    theme_classic() +
    scale_color_identity() +
    labs(x='baseMean', y='log2FoldChange') +
    ggtitle(paste0(case, ' vs ', control, ' (p < ', cutoff, ')'),
            paste('n = ', nrow(results),
                  '// n(up) = ', up,
                  '// n(down) = ', down)) +
    scale_x_log10() +
    theme(axis.line=element_line(size=0.5),
          axis.text=element_text(size=12),
          axis.ticks=element_line(size=0.5))
  
  # Generate MA-plot with highlights and labeling of significantly deregulated genes
  p.highlight <- p +
    geom_point(data=nonsig,
               aes(x=baseMean, y=log2FoldChange, col='gray60'),
               size=1.3, shape=16) +
    geom_point(data=dereg,
               aes(x=baseMean, y=log2FoldChange, col='red1'),
               size=1.3, shape=16) +
    geom_abline(aes(intercept=-1, slope=0), size=0.8, linetype=3) +
    geom_hline(yintercept=0, size=0.8) +
    geom_abline(aes(intercept=1, slope=0), size=0.8, linetype=3)
    
  p.highlight
}

#########################
#### Start of Script ####
#########################

#### Set working directory ####
setwd("~/mrc/project/sita_slam-seq/processed")

#### Set parameters ####
gse <- "PROJ1791"

control <- "Control"
treatment <- "heatshock"

alpha <- 0.05

delete.columns <- c("Control_rep1", "Control_rep2",
                    "heatshock_rep1", "heatshock_rep2",
                    "LPS_rep1", "LPS_rep2")

#### Load data ####
counts.nascent <- read.table(paste0(gse, "_counts_nascent.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)
counts.total <- read.table(paste0(gse, "_counts_total.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# PCA plot
png(glue("{gse}_PCA_Nascent.png"))
PCAPlot(counts.nascent, "PCA (Nascent)")
dev.off()

png(glue("{gse}_PCA_Total.png"))
PCAPlot(counts.total, "PCA (Total)")
dev.off()

# Delete columns
if (length(delete.columns) > 0) {
  counts.nascent <- counts.nascent[, -which(colnames(counts.nascent) %in% delete.columns)]
  counts.total <- counts.total[, -which(colnames(counts.total) %in% delete.columns)]
}

# DESeq2 analysis
results <- deAnalysis(counts.nascent, counts.total, treatment, control, alpha)

# Add Ensembl symbol
results$nascent <- addEnsemblSymbol(results$nascent)
results$total <- addEnsemblSymbol(results$total)

# Save DESeq2 output
write.table(results$nascent,
            file=glue("{gse}_DEnascent_{treatment}_vs_{control}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(results$total,
            file=glue("{gse}_DEtotal_{treatment}_vs_{control}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# MA plot and save output
png(glue("{gse}_MAplotNascent_{treatment}_vs_{control}_{alpha}.png"))
MAPlot(results$nascent, treatment, control, alpha)
dev.off()

png(glue("{gse}_MAplotTotal_{treatment}_vs_{control}_{alpha}.png"))
MAPlot(results$total, treatment, control, alpha)
dev.off()
