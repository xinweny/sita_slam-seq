#### Parsing command line ####
library(getopt)

spec <- matrix(c(
  'help', 'h', 0, 'logical', 'print the usage of the command',
  'gse', 'g', 1, 'character', 'gse accession number',
  'countFolder', 'c', 1, 'character', 'path to count file folder'
), ncol=5, byrow=TRUE)

opt <- getopt(spec)

if (!is.null(opt$help)) {
  # Get script name
  cmd <- commandArgs(FALSE)
  self <- strsplit(cmd[grep('--file', cmd)], '=')[[1]][2]
  cat(basename(self), ': create count matrix for total and nascent transcripts for all conditions in a GSE.\n\n')
  cat(getopt(spec, command=self, usage=TRUE))
  q(status=1)
}

#### Packages ####
library(tidyverse)

#### Functions ####
parseSample <- function(file, type) {
  # Get sample name
  name <- basename(file)
  name <- sub('.csv', '', name)
  
  # Get raw and normalised counts
  countTab <- read.csv(file, header=TRUE, sep="\t")
  
  if (type == 'total') {
    countTab <- countTab %>% dplyr::select(gene_name, readCount)
  } else if (type == 'nascent') {
    countTab <- countTab %>% dplyr::select(gene_name, tcReadCount)
  }
  
  names(countTab) <- c('gene_name', name)
  
  return(countTab)
}

#### Script ####
countFiles <- list.files(opt$countFolder, pattern='.csv')
counts.total <- parseSample(file.path(opt$countFolder, countFiles[1]), type='total')
counts.nascent <- parseSample(file.path(opt$countFolder, countFiles[1]), type='nascent')

if(length(countFiles) > 1) {
  for (i in 2:length(countFiles)) {
    counts.total <- counts.total %>% inner_join(parseSample(file.path(opt$countFolder, countFiles[i]), type='total'),
                                                by='gene_name')
    counts.nascent <- counts.nascent %>% inner_join(parseSample(file.path(opt$countFolder, countFiles[i]), type='nascent'),
                                                    by='gene_name')
  }
}

if (!dir.exists(file.path('data', opt$gse, 'out', 'processed'))) { dir.create(file.path('data', opt$gse, 'out', 'processed')) }

write.table(counts.total,
            file.path('data', opt$gse, 'out', 'processed', paste0(opt$gse, '_counts_total.txt')),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(counts.nascent,
            file.path('data', opt$gse, 'out', 'processed', paste0(opt$gse, '_counts_nascent.txt')),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
