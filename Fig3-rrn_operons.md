# Analysis of 16S copy numbers and ratio
```
#read in data with ratios of each taxa and the avg copy # for the phylum (from rrnDB)
#have removed candidate phyla as they're not represented in rrnDB

library(readr)
#this file is zipped on GitHub
ratio_copy <- read_delim("ratio_copy.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#spearman correlation
cor.test(ratio_copy$Phylum_mean_copyno, ratio_copy$value, method='spearman')
#r=-.00388, p=0.6

#pearson correlation
cor.test(ratio_copy$Phylum_mean_copyno, ratio_copy$value)
#r=-0.011, p=0.26


#log transform the 16S ratio values for plotting
ratio_copy$value<-log(ratio_copy$value)

#Purge rows with -Inf
library(dplyr)
ratio_copy<-filter(ratio_copy, value != '-Inf')

ratio_copy$Phylum<-gsub('p__', '', ratio_copy$Phylum)

pal26<-c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")
#plot it with ggplot
library(ggplot2)
ggplot(ratio_copy, aes(Phylum_mean_copyno, value, colour=Phylum))+
  geom_point()+
  xlab("Average Copy Number at Phylum level")+
  ylab("Log10 16S ratio")+
  coord_cartesian(ylim=c(-10,10))+
  theme_bw()+
  scale_colour_manual(values=pal26)
```
