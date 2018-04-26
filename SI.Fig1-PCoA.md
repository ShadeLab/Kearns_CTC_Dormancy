# PCoA for 16S rRNA and rRNA gene
```
#read in principal coordinates
library(readr)
stress_data <- read_delim("stress_data.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
library(ggplot2)
#plot PCoA, based on weighted UniFrac
ggplot(stress_data, aes(WU_Axis1, WU_Axis2, shape=Type, colour=Treatment))+
  geom_point(size=5)+
  theme_bw()+
  xlab("Axis1- 49.24%")+
  ylab("Axis2- 24.32%")+
  theme(axis.text = element_text(colour="black",size=16), 
        axis.title = element_text(colour="black",size=20))
```

# Adonis for DNA/RNA difference
library(readr)
bean_table <- read_delim("otu_table_tax_filt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
bean_table <-bean_table[,-1]
bean_table <-bean_table[,-49]

#rarefy table
library(vegan)
bean_table_t<-t(bean_table)
bean_table_rare<-rrarefy(bean_table_t, sample=37000)

#get sample names
samps<-row.names(bean_table_rare)
samps<-gsub('cDNA', 'R', samps)
samps<-gsub('DNA', 'R2', samps)
samps<-gsub('D', '', samps)
samps<-gsub('C', '', samps)
samps<-gsub('N', '', samps)
samps<-gsub('R2', 'D', samps)
samps<-gsub('[[:digit:]]+', '', samps)

#adonis
adonis(bean_table_rare ~ samps, permutations=10000)
#p=9.9e-5, F=12.337, R2=0.47
```
# Adonis for treatment
```
#split rarefied table into DNA and RNA tables
bean_table_rare2<-as.data.frame(t(bean_table_rare))
bean_dna<-data.frame(bean_table_rare2[,-grep('cDNA', names(bean_table_rare2))])
bean_dna_t<-t(bean_dna)

bean_rna<-data.frame(bean_table_rare2[,grep('cDNA', names(bean_table_rare2))])
bean_rna_t<-t(bean_rna)

#get names for variables (C,N,D)
rna_names<-row.names(bean_rna_t)
rna_names<-gsub('[[:digit:]]+', '', rna_names)
rna_names<-gsub('cDNA', '', rna_names)

dna_names<-row.names(bean_dna_t)
dna_names<-gsub('[[:digit:]]+', '', dna_names)
dna_names<-gsub('DNA', '', dna_names)

#adonis for RNA
library(vegan)
adonis(bean_rna_t ~ rna_names, permutations=10000)
#F=9.59,p=9.9e-5,R2=0.4774

#adonis for DNA
adonis(bean_dna_t ~ dna_names, permutations=10000)
#F=1.54,p=0.01,R2=0.128
```
