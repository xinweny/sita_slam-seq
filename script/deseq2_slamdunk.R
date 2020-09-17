#### Parsing command line ####
library(getopt)

spec <- matrix(c(
  'help', 'h', 0, 'logical', 'print the usage of the command',
  'group', 't', 2, 'character', 'group',
  'design', 'd', 2, 'character', 'path to design file',
  'countFolder', 'c', 2, 'character', 'path to count file folder',
  'pvalue', 'p', 2, 'numeric', 'p-value cutoff',
  'output', 'O', 2, 'character', 'path to output folder'
), ncol=5, byrow=TRUE)

opt <- getopt(spec)

if (!is.null(opt$help) || length(opt) == 3) {
  # Get script name
  cmd <- commandArgs(FALSE)
  self <- strsplit(cmd[grep('--file', cmd)], '=')[[1]][2]
  cat(basename(self), ': compute different transcriptional output from SLAM-seq data.\n\n')
  cat(getopt(spec, command=self, usage=TRUE))
  q(status=1)
}

if (is.null(opt$group)) { stop('arg group must be specified') }
if (is.null(opt$design)) { stop('arg design must be specified') }
if (is.null(opt$countFolder)) { stop('arg countFolder must be specified') }
if (is.null(opt$pvalue)) { opt$pvalue <- 0.1 }
if (is.null(opt$output)) { opt$output <- opt$group }

#### Packages ####
library(tidyverse)
library(DESeq2)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)

#### Functions ####
# Function to parse count table
parseSample <- function(file) {
  # Get sample name
  connection <- file(file, 'r')
  name <- readLines(con=connection, n=1)
  name <- sub('.*name:', '', name)
  close(connection)
  
  # Get raw and normalised counts
  countTab <- read_tsv(file, comment='#')
  
  countTab <- countTab %>%
    mutate(RPM = readCount / sum(readCount) * 1e6) %>%
    dplyr::select(gene_name, readCount, tcReadCount, RPM)
  
  names(countTab) <- c('gene_name', paste0(name, '_total'), paste0(name, '_tc'), paste0(name), '_RPM')
  
  return(countTab)
}

# Function to run DESeq2
deAnalysis <- function(counts, design) {
  countData.total <- counts %>%
    as.matrix()
  
  row.names(countData.total) <- counts$gene_name
  
  countData <- counts %>%
    dplyr::select(contains('_tc')) %>%
    as.matrix()
  
  rownames(countData) <- counts$gene_name
  
  sampleOrder <- sub('_tc', '', colnames(countData))
  
  design <- design[match(sampleOrder, design$name), ]
  
  dds <- DESeqDataSetFromMatrix(countData=countData,
                                colData=design,
                                design= ~ condition)
  
  row.names(dds) <- row.names(countData)
  
  dds.total <- DESeqDataSetFromMatrix(countData=countData.total,
                                      colData=design,
                                      design= ~ condition)
  
  dds.total <- dds.total[rowSums(counts(dds)) > 0, ] # remove uninformative rows
  dds <- dds[rowSums(counts(dds)) > 0, ] # remove uninformative rows
  
  # Run DESeq2 main command
  dds.total <- DESeq(dds.total)
  sizeFactors(dds) <- sizeFactors(dds.total) # apply size factors to tc counts
  
  dds <- DESeq(dds)
  
  return(dds)
}

# Function to get DESeq2 contrast
getContrast <- function(counts, dds, case, control) {
  results <- results(dds, contrast=c('condition', case, control))
  results <- data.frame(gene_name=as.character(rownames(dds)),
                        log2FC_deseq2=results$log2FoldChange,
                        padj=results=padj)
  
  results <- results[which(results$padj != 'NA'), ]
  
  # Calculate mean baseline expression (mean RPM of untreated samples)
  ctrl <- control
  ctrl <- design %>%
    dplyr::filter(condition == ctrl) %>%
    .$name %>%
    paste0(., '_RPM')
  
  avg.RPM <- data.frame(gene_name=as.character(counts$gene_name),
                        avg.RPM.ctrl=counts %>%
                          dplyr::select(all_of(ctrl)) %>%
                          rowMeans())
  
  # Intersect with baseline expression and gene symbols, then export as table
  export.deseq2 <- left_join(results, avg.RPM) %>%
    dplyr::select(gene_name, everything()) %>% # reorder gene_name as first column
    arrange(padj, log2FC_deseq2)
  
  return(export.deseq2)
}

# Function to generate MA-plot of differential gene expression calling
MAPlot <- function(export.deseq2, case, control, cutoff) {
  # Extract gene IDs for top 20 deregulated genes for plotting
  dereg <- export.deseq2[which(export.deseq2$padj <= cutoff), ] # all significant deregulated genes
  down <- nrow(dereg[dereg$log2FC_deseq2 <= 0, ]) # downregulated genes
  up <- nrow(dereg[dereg$log2FC_deseq2 >= 0, ]) # upregulated genes
  top.dereg <- dereg[order(abs(dereg$log2FC_deseq2), decreasing=TRUE)[1:20], 1]
  top.dereg <- top.dereg[!is.na(top.dereg)]
  
  # Create generic data frame for plotting by ggplot2
  x <- log10(export.deseq2$avg.RPM.ctrl) # x-axis: baseline mRNA expression
  y <- export.deseq2$log2FC_deseq2 # y-axis: log2 fold-change (treated / ctrl)
  
  gene_name <- export.deseq2$gene_name
  d <- densCols(x, y, nbin=100,
                colramp=colorRampPalette((brewer.pal(9, 'Greys')[-c(1:5)]))) # d - gray-scale density coloring for scatter plots
  df <- data.frame(x, y, d, gene_name)
  
  ## Generate unlabeled plot with accurate marginal density for use in figures
  # Generate basic MA-like plot with density coloring
  p <- ggplot(df, aes(x=x, y=y, label=gene_name)) +
    theme_classic() +
    scale_color_identity() +
    labs(x='average expression total mRNA (RPM)', y='log2FC')
  
  p.basic <- p +
    geom_point(aes(x, y, col=d), size=1.3, shape=16) +
    geom_abline(aes(intercept=-1, slope=0), size=0.8, linetype=3) +
    geom_hline(yintercept=0, size=0.8) +
    geom_abline(aes(intercept=1, slope=0), size=0.8, linetype=3)
  
  # Generate marginal density plot of fold-changes
  p.right <- ggplot(df, aes(y)) +
    geom_density(alpha=0.5, fill='gray40', size=0.8) +
    labs(x='', y='') +
    coord_flip() +
    geom_vline(xintercept=0) +
    theme_classic()
  
  grid.arrange(p.basic, p.right,
               ncol=2, widths=c(4, 1)) # assemble plots
  
  ## Generate summary of MA-plots with additional information & highlights
  # NB: Formatting of axes may shift marginal density plot from scatter plot
  # Only use exported plots for visual inspection
  
  # Generate basic MA-like plot with density coloring
  p <- ggplot(df, aes(x=x, y=y, label=gene_name)) +
    theme_classic() +
    scale_color_identity() +
    labs(x='average expression total mRNA (RPM)', y='log2FC') +
    ggtitle(paste0(case, ' vs ', control),
            paste('n = ', nrow(df),
                  '// n(up) = ', up,
                  '// n(down) = ', down)) +
  theme(axis.line=element_line(size=0.5),
        axis.text=element_text(size=12),
        axis.ticks=element_line(size=0.5))
  
  p.basic <- p +
    geom_point(aes(x, y, col=d), size=1.3, shape=16) +
    geom_abline(aes(intercept=-1, slope=0), size=0.8, linetype=3) +
    geom_hline(yintercept=0, size=0.8) +
    geom_abline(aes(intercept=1, slope=0), size=0.8, linetype=3)
  
  p.right <- ggplot(df, aes(y)) +
    geom_density(alpha=0.5, fill='gray40', size=0.8) +
    labs(x='', y='') +
    coord_flip() +
    geom_vline(xintercept=0) +
    theme_classic() +
    ggtitle('', '') +
    theme(legend.position='none',
          text=element_text(color='white'),
          axis.text=element_text(color='white'),
          axis.ticks.x=element_line(color='white'),
          axis.line.x=element_line(color='white'))
  
  # Generate MA-plot with highlights and labeling of significantly deregulated genes
  if (nrow(dereg) > 0) {
    p.highlight <- p +
      geom_point(data=df[!(df$gene_name %in% dereg$gene_name), ],
                 aes(x, y, col='gray60'), size=1.3, shape=16) +
      geom_point(data=df[df$gene_name %in% dereg$gene_name, ],
                 aes(x, y, col='red1'), size=1.3, shape=16) +
      geom_abline(aes(intercept=-1, slope=0), size=0.8, linetype=3) +
      geom_hline(yintercept=0, size=0.8) +
      geom_abline(aes(intercept=1, slope=0), size=0.8, linetype=3)
    
    p.highlight.2 <- p.highlight +
      geom_label_repel(data=df[df$gene_name %in% top.dereg, ])
    
    # Export plot
    grid.arrange(p.basic, p.right,
                 ncol=2, widths=c(4, 1)) # assemble plots
    grid.arrange(p.highlight, p.right,
                 ncol=2, widths=c(4, 1)) # add highlights
    grid.arrange(p.highlight.2, p.right,
                 ncol=2, widths=c(4, 1)) # add labels
  }
}

#### Read design and counts ####
design <- read_tsv(opt$design)

design <- design %>%
  filter(group == opt$group)

if (nrow(design) == length(unique(design$condition))) {
  quit(save='no', status=0, runLast=TRUE) # stop if no replicates
}

countFiles <- list.files(opt$countFolder, pattern='.csv')
counts <- parseSample(file.path(opt$countFolder, countFiles[1]))

if(length(countFiles) > 1) {
  for (i in 2:length(countFiles)) {
    counts <- counts %>% inner_join(parseSample(file.path(opt$countFolder, countFiles[i])),
                                    by='gene_name')
  }
}

#### Run DESeq2 ####
dds <- deAnalysis(counts, design)

# Export PCA plot (default DESeq2 PCA on 500 most variable genes)
if (!dir.exists(opt$output)) { dir.create(opt$output) }

pdf(file.path(opt$output, 'PCA.pdf'))
plotPCA(varianceStabilizingTransformation(dds), intgroup='condition')
dev.off()

ctrl <- design %>%
  filter(control == 1) %>%
  .$condition %>%
  unique()

cases <- design %>%
  filter(control == 0) %>%
  .$condition %>%
  unique()

for (case in cases) {
  if (!dir.exists(file.path(opt$output, case))) { dir.create(file.path(opt$output, case)) }
  
  #### Extract contrasts ####
  export.deseq2 <- getContrast(counts, dds, case, ctrl)
  
  write.tsv(export.deseq2,
            file.path(opt$output, case, 'DESeq2.txt'))
  
  #### Plot results ####
  pdf(file.path(opt$output, case, 'MAPlot.pdf'))
  MAPlot(export.deseq2, case, ctrl, opt$pvalue)
  dev.off()
}
