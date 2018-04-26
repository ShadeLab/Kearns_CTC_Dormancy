# Rarefaction curves

```
#read in unrarefied data
library(readr)
bean_table <- read_delim("otu_table_tax_filt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#extracttaxonomy and OTU labels
bean_otus<-bean_table[,1]
str(bean_otus)
View(bean_otus)
bean_taxonomy<-bean_table[,50]
str(bean_taxonomy)

#remove taxonomy/OTU names
bean_table<-bean_table[,c(-1,-50)]
str(bean_table)


#determine the lowest sampling depth
min(colSums(bean_table))
max(colSums(bean_table))
#lowest depth is 37123

#rarefy
library(vegan)
bean_table_t<-t(bean_table)
#rarefy table
bean_rare_high<-rrarefy(bean_table_t, sample=37123)
#make rarefaction curves with Phyloseq
library(phyloseq)
library(tidyverse)

bean_curves2<-rarecurve(bean_table_t, step=100)

rare <- lapply(bean_curves2, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- rownames(bean_table_t)

rare2 <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

#add factors to data frame
rare2$name<-rare2$sample
rare2$name<-gsub("cDNA", "", rare2$name)
rare2$name<-gsub("DNA", "", rare2$name)
rare2$name<-gsub("[[:digit:]]+", "", rare2$name)
rare2$name<-as.factor(rare2$name)

#add type (DNA/RNA)
rare2$type<-rare2$sample
rare2$type<-gsub("cDNA", "RRA", rare2$type)
rare2$type<-gsub("DNA", "ZA", rare2$type)
rare2$type<-gsub("[[:digit:]]+", "", rare2$type)
rare2$type<-gsub("C", "", rare2$type)
rare2$type<-gsub("D", "", rare2$type)
rare2$type<-gsub("N", "", rare2$type)
rare2$type<-gsub("RRA", "RNA", rare2$type)
rare2$type<-gsub("ZA", "DNA", rare2$type)

#add treatment
rare2$treat<-rare2$sample
rare2$treat<-gsub("cDNA", "", rare2$treat)
rare2$treat<-gsub("DNA", "", rare2$treat)
rare2$treat<-gsub("[[:digit:]]+", "", rare2$treat)

rare2$treat<-gsub("C", "Control", rare2$treat)
rare2$treat<-gsub("D", "Drought", rare2$treat)
rare2$treat<-gsub("N", "Nutrients", rare2$treat)

rare2$type<-gsub("C", "", rare2$type)
rare2$type<-gsub("D", "", rare2$type)
rare2$type<-gsub("N", "", rare2$type)

library(ggplot2)
#all
ggplot(rare2, aes(raw.read, OTU, colour=sample))+
  geom_line()+
  theme_bw()+
  ylab("OTUs Observed")+
  xlab("Sequencing Depth")+
  geom_vline(xintercept=37123)
#faceted by DNA/RNA
ggplot(rare2, aes(raw.read, OTU, colour=sample))+
  geom_line()+
  theme_bw()+
  ylab("OTUs Observed")+
  xlab("Sequencing Depth")+
  facet_wrap(~type)+
  geom_vline(xintercept=37123)

#split rarefaction by treatment
ggplot(rare2, aes(raw.read, OTU, colour=sample))+
  geom_line()+
  theme_bw()+
  ylab("OTUs Observed")+
  xlab("Sequencing Depth")+
  facet_wrap(~treat)+
  geom_vline(xintercept=37123)+
  theme(legend.position="none")
```
