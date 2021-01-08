#### Packages ####
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(ggrepel)
library(glue)

#### Set parameters ####
gse <- "GSE100708"
selectConditions <- c("K562_AID-BRD4_IAA", "K562_AID-BRD4_DMSO") # treatment / control
alpha <- 0.01

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
  
  row.names(geneList) <- geneList[, 1]
  geneList[, 1] <- NULL
  
  table$geneSymbol <- geneList[, 1][match(rownames(table), rownames(geneList))]
  newTable <- table
  
  return(newTable)
}

deAnalysis <- function(counts.nascent, counts.total, selectConditions, a) {
  countData.nascent <- counts.nascent %>% dplyr::select(matches(selectConditions)) %>% as.matrix(rownames=1)
  countData.total <- counts.total %>% dplyr::select(matches(selectConditions)) %>% as.matrix(rownames=1)
  
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
  
  dds.nascent$condition <- relevel(dds.nascent$condition, ref=conditions[2])
  dds.total$condition <- relevel(dds.total$condition, ref=conditions[2])
  
  # Run DESeq2 main command
  dds.total <- DESeq(dds.total)
  sizeFactors(dds.nascent) <- sizeFactors(dds.total) # apply size factors to tc counts
  
  dds.nascent <- DESeq(dds.nascent)
  
  coef <- tail(resultsNames(dds.nascent), n=1)
  
  # Get results
  res.nascent <- results(dds.nascent, name=coef, alpha=a)
  resLFC.nascent <- lfcShrink(dds.nascent, coef=coef, res=res.nascent)
  
  res.total <- results(dds.total, name=coef, alpha=a)
  resLFC.total <- lfcShrink(dds.total, coef=coef, res=res.total)
  
  # Put into list and return
  resList <- list()
  resList$nascent <- as.data.frame(resLFC.nascent) %>% arrange(padj, desc(log2FoldChange))
  resList$total <- as.data.frame(resLFC.total) %>% arrange(padj, desc(log2FoldChange))
  resList$conditions <- conditions
  
  return(resList)
}

MAPlot <- function(results, case, control, cutoff) {
  # Extract gene IDs for top 20 deregulated genes for plotting
  results$gene_name <- row.names(results)
  
  dereg <- results[which(results$padj <= cutoff), ] # all significant deregulated genes
  nonsig <- results[!(results$gene_name %in% dereg$gene_name), ] # non-significant deregulated genes
  
  down <- nrow(dereg[dereg$log2FoldChange < 0, ]) # downregulated genes
  up <- nrow(dereg[dereg$log2FoldChange > 0, ]) # upregulated genes
  
  top.dereg <- dereg[order(abs(dereg$log2FoldChange), decreasing=TRUE)[1:20], 1]
  top.dereg <- top.dereg[!is.na(top.dereg)]
  
  ## Generate summary of MA-plots with additional information & highlights
  # NB: Formatting of axes may shift marginal density plot from scatter plot
  # Only use exported plots for visual inspection
  
  # Generate basic MA-like plot with density coloring
  p <- results %>% ggplot(aes(x=log10(baseMean), y=log2FoldChange)) +
    theme_classic() +
    scale_color_identity() +
    labs(x='baseMean', y='log2FoldChange') +
    ggtitle(paste0(case, ' vs ', control),
            paste('n = ', nrow(results),
                  '// n(up) = ', up,
                  '// n(down) = ', down)) +
    scale_x_log10() +
    theme(axis.line=element_line(size=0.5),
          axis.text=element_text(size=12),
          axis.ticks=element_line(size=0.5))
  
  # Generate MA-plot with highlights and labeling of significantly deregulated genes
  if (nrow(dereg) > 0) {
    p.highlight <- p +
      geom_point(data=nonsig,
                 aes(x=baseMean, y=log2FoldChange, col='gray60'),
                 size=1.3, shape=16) +
      geom_point(data=dereg,
                 aes(x=baseMean, y=log2FoldChange, col='red1'),
                 size=1.3, shape=16) +
      geom_abline(aes(intercept=-1, slope=0), size=0.8, linetype=3) +
      geom_hline(yintercept=0, size=0.8) +
      geom_abline(aes(intercept=1, slope=0), size=0.8, linetype=3) +
      scale_x_log10()
    
    p.highlight
  }
}

#### Set working directory ####
setwd("/Users/Pomato/mrc/project/sita_slam-seq/processed")

#########################
#### Start of Script ####
#########################

counts.nascent <- read.table(paste0(gse, "_counts_nascent.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)
counts.total <- read.table(paste0(gse, "_counts_total.txt"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# DESeq2 analysis
resultsLFC <- deAnalysis(counts.nascent, counts.total, selectConditions, alpha)

# Add Ensembl symbol
resultsLFC$nascent <- addEnsemblSymbol(resultsLFC$nascent)
resultsLFC$total <- addEnsemblSymbol(resultsLFC$total)

write.table(resultsLFC$nascent,
            file=glue("{gse}_deNascent_{resultsLFC$conditions[1]}vs{resultsLFC$conditions[2]}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t")
write.table(resultsLFC$total,
            file=glue("{gse}_deTotal_{resultsLFC$conditions[1]}vs{resultsLFC$conditions[2]}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t")

# MA plot
png(glue("{gse}_MAplotNascent_{resultsLFC$conditions[1]}.{resultsLFC$conditions[2]}.png"))
MAPlot(resultsLFC$nascent, resultsLFC$conditions[1], resultsLFC$conditions[2], alpha)
dev.off()

png(glue("{gse}_MAplotTotal_{resultsLFC$conditions[1]}.{resultsLFC$conditions[2]}.png"))
MAPlot(resultsLFC$total, resultsLFC$conditions[1], resultsLFC$conditions[2], alpha)
dev.off()