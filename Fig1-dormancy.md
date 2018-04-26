# Fig. 1A- Box-plot of CTC and 16S dormancy
```
#read in CTC dormancy data
stress_dorm_ctc16 <- read.delim("~/stress_dorm_ctc16.txt")
library(ggplot2)
library(reshape2)
#remove non-dormancy things
stress_dormancy<-stress_dorm_ctc16[,c(-3,-4,-5,-8)]

#melt and plot
stress_dorm_m<-melt(stress_dormancy)
stress_dorm_m$active<-100-stress_dorm_m$value


act1<-ggplot(stress_dorm_m, aes(Type, active, fill=variable))+
  geom_boxplot()+
  theme_bw()+
  ylab("Percent Active")+
  xlab("")+
  coord_cartesian(ylim=c(0,100))+
  scale_fill_manual(values=c('grey', 'white'))

#for CTC
bartlett.test(CTC_Dorm ~ Type, data=stress_dormancy)
#p=0.0003
pairwise.t.test(stress_dormancy$CTC_Dorm, stress_dormancy$Type, p.adjust='hochberg')
#            Control Drought
#Drought   < 2e-16 -      
#Nutrients 6.6e-15     0.34   

t.test(stress_dormancy$CTC_Dorm)


#for 16S rRNA
bartlett.test(X16s_Dorm ~ Type, data=stress_dormancy)
#p=0.91
S16_aov<-aov(X16s_Dorm ~ Type, data=stress_dormancy)
summary(S16_aov)
#p=2e-16, F=120.5
```

# Fig. 1B threshold activity
```
#read in unrarefied data
setwd("/Users/patty/Dropbox/R/bean_nitrogen_drought/")
library(readr)
bean_table <- read_delim("otu_table_tax_filt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#remove taxonomy and OTU labels
bean_otus<-bean_table[,grep('OTU', names(bean_table))]
bean_taxonomy<-bean_table[,grep('taxonomy', names(bean_table))]
bean_table<-bean_table[,c(-1)]
bean_table<-bean_table[,c(-49)]

#rarefy data
library(vegan)
bean_table_t<-t(bean_table)
bean_table2<-rrarefy(bean_table_t, sample=37123)
bean_table2<-data.frame(t(bean_table2))

#extract 16S rRNA and rRNA gene data, write to new frame
bean_rna<-data.frame(bean_table[,grep('cDNA', names(bean_table))])
#str(bean_rna)
bean_dna<-data.frame(bean_table[,-grep('cDNA', names(bean_table))])
#str(bean_dna)

#reorder tables to be in same way
bean_rna<-bean_rna[,order(names(bean_rna))]
bean_dna<-bean_dna[,order(names(bean_dna))]

#function
dormancy_calc<-function(S16rRNA, S16rRNAG, num){
  #S16rRNA is a table of 16S rRNA data
  #S16rRNAG is a table of 16S rRNA gene data
  #num is the ratio threshold for a tax to be active, there's no default
  #add 1 to each DNA OTU, so no 0's in denominator
  S16rRNAG<-replace(S16rRNAG, S16rRNAG == 0, 1)
  #calculate ratio for each OTU
  v.rat<-S16rRNA/S16rRNAG
  #calculate the no. OTUs for each column
  total.otus<-apply((S16rRNAG), 2, function(x) length(which(x>1)))
  #count active OTUs
  active.otus<-apply(v.rat, 2, function(x) length(which(x>num)))
  #calculate percent active
  per.active<-100*(active.otus/total.otus)
  #calculate percent dormant
  per.dorm<-100-per.active
  #add threshold to the output table
  ratio_threshold<-rep(c(num), times=length(per.dorm))
  #calculate sampling depth
  samp.depth<-rep(min(colSums(bean_table)), times=length(per.dorm))
  #calculate average ratio
  avg.ratio<-apply(v.rat, 2, function(x) mean(x))
  #max ratio
  max.ratio<-apply(v.rat, 2, function(x) max(x))
  #min ratio
  min.ratio<-apply(v.rat, 2, function(x) min(x))
  #standard deviation
  sd.ratio<-apply(v.rat, 2, function(x) sd(x))
  #standard error
  se.ratio<-apply(v.rat, 2, function(x) sd(x)/sqrt(length(x)))
  #median ratio
  median.ratio<-apply(v.rat, 2, function(x) median(x))
  #print the results!!
  as.data.frame(cbind(total.otus,active.otus, per.active,per.dorm, ratio_threshold, samp.depth, avg.ratio, median.ratio, max.ratio, min.ratio, sd.ratio, se.ratio))
}

bean_dormancy1<-dormancy_calc(bean_rna, bean_dna, 1)
bean_dormancy2<-dormancy_calc(bean_rna, bean_dna, 2)
bean_dormancy3<-dormancy_calc(bean_rna, bean_dna, 3)
bean_dormancy4<-dormancy_calc(bean_rna, bean_dna, 4)
bean_dormancy5<-dormancy_calc(bean_rna, bean_dna, 5)
bean_dormancy6<-dormancy_calc(bean_rna, bean_dna, 6)
bean_dormancy7<-dormancy_calc(bean_rna, bean_dna, 7)
bean_dormancy8<-dormancy_calc(bean_rna, bean_dna, 8)
bean_dormancy9<-dormancy_calc(bean_rna, bean_dna, 9)
bean_dormancy10<-dormancy_calc(bean_rna, bean_dna, 10)
bean_dormancy15<-dormancy_calc(bean_rna, bean_dna, 15)
bean_dormancy20<-dormancy_calc(bean_rna, bean_dna, 20)
bean_dormancy25<-dormancy_calc(bean_rna, bean_dna, 25)
bean_dormancy30<-dormancy_calc(bean_rna, bean_dna, 30)
bean_dormancy35<-dormancy_calc(bean_rna, bean_dna, 35)
bean_dormancy40<-dormancy_calc(bean_rna, bean_dna, 40)
bean_dormancy45<-dormancy_calc(bean_rna, bean_dna, 45)
bean_dormancy50<-dormancy_calc(bean_rna, bean_dna, 50)

#catenate the files
bean_dorm_titration<-as.data.frame(rbind(bean_dormancy1, bean_dormancy2, bean_dormancy3, bean_dormancy4, bean_dormancy5, 
                                         bean_dormancy6, bean_dormancy7, bean_dormancy8, bean_dormancy9, bean_dormancy10,
                                         bean_dormancy15, bean_dormancy20, bean_dormancy25, bean_dormancy30, bean_dormancy35, 
                                         bean_dormancy40, bean_dormancy45, bean_dormancy50))

#make variables for treatments
bean_dorm_titration$treatment<-row.names(bean_dorm_titration)
bean_dorm_titration$treatment<-gsub('cDNA', '', bean_dorm_titration$treatment)
bean_dorm_titration$treatment<-gsub('DNA', '', bean_dorm_titration$treatment)
bean_dorm_titration$treatment<-gsub('[[:digit:]]+', '', bean_dorm_titration$treatment)

#plot the data in ggplot
library(ggplot2)
ratio_titration<- ggplot(bean_dorm_titration, aes(ratio_threshold, per.active, colour=treatment))+
  geom_smooth()+
  ylab("Percent Active")+
  coord_cartesian(ylim=c(0,100))+
  theme_bw()+
  xlab("16S Ratio Threshold")
```
# Plot both together
```
#plot bot together
library(gridExtra)
grid.arrange(act1, ratio_titration, ncol=2)

```
