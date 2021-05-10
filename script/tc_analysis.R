setwd("/Users/pomato/mrc/project/sita_slam-seq/")

#### Packages ####
require(ggplot2)

#### Load data ####
tc <- read.table("./processed/PROJ1791_overallTCrates.txt",
                 header=TRUE, sep="\t", check.names=FALSE)

#### Visualise data ####
plot <- ggplot(data=tc, aes(x=group, y=TC, color=factor(rep))) +
  geom_point(position=position_dodge(width=0.7)) +
  scale_colour_manual(values=rep("black", 5)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar",width=0.05, colour="red") +
  stat_summary(fun=mean, fun.max=mean, fun.min=mean, colour="red", geom="crossbar", width=0.2) +
  theme(legend.position="none") +
  xlab("") +
  ylab("Overall % T>C conversion rate")

png("./processed/PROJ1791_overallTCrates.png")
plot
dev.off()