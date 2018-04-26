## Worflow for creating histograms of phylum-level distribution of 16S ratio
### R (v.3.4.1)

```
#read in unrarefied data
library(readr)
bean_table <- read_delim("otu_table_tax_filt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#extract taxonomy and OTU labels
bean_otus<-bean_table[,1]
bean_taxonomy<-bean_table[,50]


#remove taxonomy/OTU names
bean_table<-bean_table[,c(-1,-50)]
str(bean_table)


#determine the lowest sampling depth
min(colSums(bean_table))
max(colSums(bean_table))
#lowest depth is 37123, max is 64163

#rarefy
library(vegan)
bean_table_t<-t(bean_table)
#rarefy table
bean_rare_high<-rrarefy(bean_table_t, sample=37123)
bean_rare_high<-data.frame(t(bean_rare_high))
str(bean_rare_high)

#add OTU names and taxonomy
bean_rare_high2<-cbind(bean_otus, bean_taxonomy, bean_rare_high)


#remove rows with all 0's
bean_table2<-bean_rare_high2[which(rowSums(bean_rare_high) > 0),] 
str(bean_table2)

#write taxonomy/OTU name to vectors
bean_taxonomy2<-bean_table2[,2]
bean_otus2<-bean_table2[,1]
#remove taxonomy/OTU names
bean_table2<-bean_table2[,c(-1,-2)]


#extract 16S rRNA and rRNA gene data, write to new frame
bean_rna<-as.data.frame(bean_table2[,grep('cDNA', names(bean_table2))])
bean_dna<-as.data.frame(bean_table2[,-grep('cDNA', names(bean_table2))])
str(bean_rna)
str(bean_dna)

#reorder tables to be in same way
bean_rna<-bean_rna[,order(names(bean_rna))]
bean_dna<-bean_dna[,order(names(bean_dna))]

#dd 1's to zeroes in DNA table
bean_dna<-replace(bean_dna, bean_dna == 0, 1)

#calculate ratio
bean.ratio<-bean_rna/bean_dna
bean.ratio<-cbind(bean_taxonomy2, bean_otus2, bean.ratio)

#purge 'cDNA' from the table
colnames(bean.ratio)<-sub("cDNA", "", colnames(bean.ratio))

#split ratio table by treatment
nitrogen.ratio<-as.data.frame(bean.ratio[,grep('N', names(bean.ratio))])
nitrogen.ratio<-cbind(bean_taxonomy2, bean_otus2, nitrogen.ratio)

control.ratio<-as.data.frame(bean.ratio[,grep('C', names(bean.ratio))])
control.ratio<-cbind(bean_taxonomy2, bean_otus2, control.ratio)

drought.ratio<-as.data.frame(bean.ratio[,grep('D', names(bean.ratio))])
drought.ratio<-cbind(bean_taxonomy2, bean_otus2, drought.ratio)

#melt data frames for histogram
library(reshape2)
nitrogen.ratio.m<-melt(nitrogen.ratio)
drought.ratio.m<-melt(drought.ratio)
control.ratio.m<-melt(control.ratio)

#mean ratios
mean(nitrogen.ratio.m$value) 
#0.73
mean(drought.ratio.m$value) 
#0.79
mean(control.ratio.m$value)
#0.83
mean(c(0.73, 0.79, 0.83))
#0.78

#split columns by taxonomic ID
library(tidyr)
unique(nitrogen.ratio.m$Phylum)
nitrogen.ratio.m<-separate(nitrogen.ratio.m, bean_taxonomy2, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
max(nitrogen.ratio.m$value)
#406 is max ratio

control.ratio.m<-separate(control.ratio.m, bean_taxonomy2, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
max(control.ratio.m$value)
#max ratio is 649

drought.ratio.m<-separate(drought.ratio.m, bean_taxonomy2, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
max(drought.ratio.m$value)
#max ratio is 381

#full dataset
full.ratio.m<-rbind(nitrogen.ratio.m,control.ratio.m,drought.ratio.m)
full.ratio.m$variable<-gsub('[[:digit:]]+', '', full.ratio.m$variable)
#write.csv(full.ratio.m, 'full.ratio.m.csv')

#collapse the rare phyla down to 'others' 
full.ratio.collap<-full.ratio.m
unique(full.ratio.collap$Phylum)

full.ratio.collap$Phylum<-gsub('Crenarchaeota', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Parvarchaeota', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('FBP', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Fibrobacteres', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('BRC1', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('BHI80-139', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Chlorobi', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Elusimicrobia', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('OP3', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Kazan-3B-28', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Euryarchaeota', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('GN02', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Thermi', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Chlamydiae', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Nitrospirae', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('WS2', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('OD1', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('NKB19', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('Tenericutes', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('WPS-2', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('TM6', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('OP11', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('TM7', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('WS3', '', full.ratio.collap$Phylum)
full.ratio.collap$Phylum<-gsub('[[]]', '', full.ratio.collap$Phylum)


unique(full.ratio.collap$Phylum)

hist_colour<-rainbow(37, s=.6, v=.9)[sample(1:37,37)]

hist_colour2<-c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

ggplot(full.ratio.collap, aes(value, fill=Phylum))+
  geom_histogram(binwidth=0.5)+
  xlab("Log10 16S Ratio")+
  ylab("Count")+
  theme_bw()+
  scale_x_log10()+
  scale_fill_manual(values=hist_colour2)+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(~variable, scales='free')


#purge NAs
full.ratio.collap<-filter(full.ratio.collap, Phylum != 'NA')
hist_colour2<-c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#plot faceted by phylum
  ggplot(full.ratio.collap, aes(value, fill=Phylum))+
  geom_histogram(binwidth=0.5)+
  xlab("Log10 16S Ratio")+
  ylab("Count")+
  theme_bw()+
  scale_x_log10()+
  scale_fill_manual(values=hist_colour)+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(~Phylum, scales='free_y')
```
