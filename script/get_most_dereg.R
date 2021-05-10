#### Packages ####
library(dplyr)
library(glue)

#### Config ####
# Set working directory
setwd("~/mrc/project/sita_slam-seq/processed")

# Parameters
proj <- "PROJ1791"
set <- "SET2"
treatment <- "heatshock"
control <- "Control"

# Load DESeq2 output tables
deseq.nascent <- read.table(glue("{proj}_{set}_DEnascent_{treatment}_vs_{control}.txt"), 
                            header=TRUE, sep="\t", row.names=1, check.names=FALSE)
deseq.total <- read.table(glue("{proj}_{set}_DEtotal_{treatment}_vs_{control}.txt"),
                          header=TRUE, sep="\t", row.names=1, check.names=FALSE)

downreg.nascent <- deseq.nascent %>% filter(log2FoldChange < 0)
downreg.total <- deseq.total %>% filter(log2FoldChange < 0)

write.table(downreg.nascent,
            file=glue("{proj}_{set}_nascentDOWNREG_{treatment}_vs_{control}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(downreg.total,
            file=glue("{proj}_{set}_totalDOWNREG_{treatment}_vs_{control}.txt"),
            row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)



