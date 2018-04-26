# Effect of sequencing depth of % active
```
############For all samples
setwd("/Users/patty/Dropbox/R/bean_nitrogen_drought/")
#read in unrarefied data
library(readr)
bean_table <- read_delim("otu_table_tax_filt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#extracttaxonomy and OTU labels
bean_otus<-bean_table[,1]
str(bean_otus)
bean_taxonomy<-bean_table[,50]
str(bean_taxonomy)

#remove taxonomy/OTU names
bean_table<-bean_table[,c(-1,-50)]
str(bean_table)

#determine the lowest/max sampling depth
min(colSums(bean_table))
max(colSums(bean_table))
#lowest depth is 37123, max is 64163


#function to split table by DNA/RNA, rarefy the datam and calculate percent dormant/active and stats
dormancy_calc_rare<-function(full.table, num, dep){
  #full.table is a OTU table, no tax/names
  #dep si sampling depth to be used
  #num is the ratio threshold for a tax to be active, there's no default
  #transpose table for rarefaction
  table.t<-t(full.table)
  #rarefy table to depth 'dep', a number that cannot be higher than your lowest sequencing depth
  library(vegan)
  table.r<-rrarefy(table.t, sample=dep)
  #re-transpose the table
  table.r.t<-data.frame(t(table.r))
  #purge all OTUs with all zeroes
  table.1<-table.r.t[which(rowSums(table.r.t) > 0),] 
  #extract 16S rRNA and rRNA gene data, write to new frame
  table.rna<-as.data.frame(table.1[,grep('cDNA', names(table.1))])
  table.dna<-as.data.frame(table.1[,-grep('cDNA', names(table.1))])
  #reorder tables to be in sampe order
  table.dna<-table.dna[,order(names(table.dna))]
  table.rna<-table.rna[,order(names(table.rna))]
  #add 1 to each DNA OTU, so no 0's in denominator
  table.dna2<-replace(table.dna, table.dna == 0, 1)
  #calculate ratio for each OTU
  v.rat<-table.rna/table.dna2
  #calculate the no. OTUs for each column
  total.otus<-apply((table.dna), 2, function(x) length(which(x>1)))
  #count active OTUs
  active.otus<-apply(v.rat, 2, function(x) length(which(x>num)))
  #calculate percent active
  per.active<-100*(active.otus/total.otus)
  #calculate percent dormant
  per.dorm<-100-per.active
  #add threshold to the output table
  ratio_threshold<-rep(c(num), times=length(per.dorm))
  #calculate sampling depth
  samp.depth<-rep(min(colSums(table.dna)), times=length(per.dorm))
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
set.seed(5001)

#calculate with a ratio threshold of 1
test<-dormancy_calc_rare(bean_table, 1, 100)
test0<-dormancy_calc_rare(bean_table, 1, 500)
test1<-dormancy_calc_rare(bean_table, 1, 1000)
test2<-dormancy_calc_rare(bean_table, 1, 2000)
test3<-dormancy_calc_rare(bean_table, 1, 3000)
test4<-dormancy_calc_rare(bean_table, 1, 4000)
test5<-dormancy_calc_rare(bean_table, 1, 5000)
test6<-dormancy_calc_rare(bean_table, 1, 6000)
test7<-dormancy_calc_rare(bean_table, 1, 7000)
test8<-dormancy_calc_rare(bean_table, 1, 8000)
test9<-dormancy_calc_rare(bean_table, 1, 9000)
test10<-dormancy_calc_rare(bean_table, 1, 10000)
test11<-dormancy_calc_rare(bean_table, 1, 11000)
test12<-dormancy_calc_rare(bean_table, 1, 12000)
test13<-dormancy_calc_rare(bean_table, 1, 13000)
test14<-dormancy_calc_rare(bean_table, 1, 14000)
test15<-dormancy_calc_rare(bean_table, 1, 15000)
test16<-dormancy_calc_rare(bean_table, 1, 16000)
test17<-dormancy_calc_rare(bean_table, 1, 17000)
test18<-dormancy_calc_rare(bean_table, 1, 18000)
test19<-dormancy_calc_rare(bean_table, 1, 19000)
test20<-dormancy_calc_rare(bean_table, 1, 20000)
test21<-dormancy_calc_rare(bean_table, 1, 21000)
test22<-dormancy_calc_rare(bean_table, 1, 22000)
test23<-dormancy_calc_rare(bean_table, 1, 23000)
test24<-dormancy_calc_rare(bean_table, 1, 24000)
test25<-dormancy_calc_rare(bean_table, 1, 25000)
test26<-dormancy_calc_rare(bean_table, 1, 26000)
test27<-dormancy_calc_rare(bean_table, 1, 27000)
test28<-dormancy_calc_rare(bean_table, 1, 28000)
test29<-dormancy_calc_rare(bean_table, 1, 29000)
test30<-dormancy_calc_rare(bean_table, 1, 30000)
test31<-dormancy_calc_rare(bean_table, 1, 31000)
test32<-dormancy_calc_rare(bean_table, 1, 32000)
test33<-dormancy_calc_rare(bean_table, 1, 33000)
test34<-dormancy_calc_rare(bean_table, 1, 34000)
test35<-dormancy_calc_rare(bean_table, 1, 35000)
test36<-dormancy_calc_rare(bean_table, 1, 36000)
test37<-dormancy_calc_rare(bean_table, 1, 37000)
test372<-dormancy_calc_rare(bean_table, 1, 37123)

#threshold of 2
test02<-dormancy_calc_rare(bean_table, 2, 100)
test020<-dormancy_calc_rare(bean_table, 2, 500)
test021<-dormancy_calc_rare(bean_table, 2, 1000)
test022<-dormancy_calc_rare(bean_table, 2, 2000)
test023<-dormancy_calc_rare(bean_table, 2, 3000)
test024<-dormancy_calc_rare(bean_table, 2, 4000)
test025<-dormancy_calc_rare(bean_table, 2, 5000)
test026<-dormancy_calc_rare(bean_table, 2, 6000)
test027<-dormancy_calc_rare(bean_table, 2, 7000)
test028<-dormancy_calc_rare(bean_table, 2, 8000)
test029<-dormancy_calc_rare(bean_table, 2, 9000)
test0210<-dormancy_calc_rare(bean_table, 2, 10000)
test0211<-dormancy_calc_rare(bean_table, 2, 11000)
test0212<-dormancy_calc_rare(bean_table, 2, 12000)
test0213<-dormancy_calc_rare(bean_table, 2, 13000)
test0214<-dormancy_calc_rare(bean_table, 2, 14000)
test0215<-dormancy_calc_rare(bean_table, 2, 15000)
test0216<-dormancy_calc_rare(bean_table, 2, 16000)
test0217<-dormancy_calc_rare(bean_table, 2, 17000)
test0218<-dormancy_calc_rare(bean_table, 2, 18000)
test0219<-dormancy_calc_rare(bean_table, 2, 19000)
test0220<-dormancy_calc_rare(bean_table, 2, 20000)
test0221<-dormancy_calc_rare(bean_table, 2, 21000)
test0222<-dormancy_calc_rare(bean_table, 2, 22000)
test0223<-dormancy_calc_rare(bean_table, 2, 23000)
test0224<-dormancy_calc_rare(bean_table, 2, 24000)
test0225<-dormancy_calc_rare(bean_table, 2, 25000)
test0226<-dormancy_calc_rare(bean_table, 2, 26000)
test0227<-dormancy_calc_rare(bean_table, 2, 27000)
test0228<-dormancy_calc_rare(bean_table, 2, 28000)
test0229<-dormancy_calc_rare(bean_table, 2, 29000)
test0230<-dormancy_calc_rare(bean_table, 2, 30000)
test0231<-dormancy_calc_rare(bean_table, 2, 31000)
test0232<-dormancy_calc_rare(bean_table, 2, 32000)
test0233<-dormancy_calc_rare(bean_table, 2, 33000)
test0234<-dormancy_calc_rare(bean_table, 2, 34000)
test0235<-dormancy_calc_rare(bean_table, 2, 35000)
test0236<-dormancy_calc_rare(bean_table, 2, 36000)
test0237<-dormancy_calc_rare(bean_table, 2, 37000)
test02372<-dormancy_calc_rare(bean_table, 2, 37123)

#threshold of 5
test05<-dormancy_calc_rare(bean_table, 5, 100)
test050<-dormancy_calc_rare(bean_table, 5, 500)
test051<-dormancy_calc_rare(bean_table, 5, 1000)
test052<-dormancy_calc_rare(bean_table, 5, 2000)
test053<-dormancy_calc_rare(bean_table, 5, 3000)
test054<-dormancy_calc_rare(bean_table, 5, 4000)
test055<-dormancy_calc_rare(bean_table, 5, 5000)
test056<-dormancy_calc_rare(bean_table, 5, 6000)
test057<-dormancy_calc_rare(bean_table, 5, 7000)
test058<-dormancy_calc_rare(bean_table, 5, 8000)
test059<-dormancy_calc_rare(bean_table, 5, 9000)
test0510<-dormancy_calc_rare(bean_table, 5, 10000)
test0511<-dormancy_calc_rare(bean_table, 5, 11000)
test0512<-dormancy_calc_rare(bean_table, 5, 12000)
test0513<-dormancy_calc_rare(bean_table, 5, 13000)
test0514<-dormancy_calc_rare(bean_table, 5, 14000)
test0515<-dormancy_calc_rare(bean_table, 5, 15000)
test0516<-dormancy_calc_rare(bean_table, 5, 16000)
test0517<-dormancy_calc_rare(bean_table, 5, 17000)
test0518<-dormancy_calc_rare(bean_table, 5, 18000)
test0519<-dormancy_calc_rare(bean_table, 5, 19000)
test0520<-dormancy_calc_rare(bean_table, 5, 20000)
test0521<-dormancy_calc_rare(bean_table, 5, 21000)
test0522<-dormancy_calc_rare(bean_table, 5, 22000)
test0523<-dormancy_calc_rare(bean_table, 5, 23000)
test0524<-dormancy_calc_rare(bean_table, 5, 24000)
test0525<-dormancy_calc_rare(bean_table, 5, 25000)
test0526<-dormancy_calc_rare(bean_table, 5, 26000)
test0527<-dormancy_calc_rare(bean_table, 5, 27000)
test0528<-dormancy_calc_rare(bean_table, 5, 28000)
test0529<-dormancy_calc_rare(bean_table, 5, 29000)
test0530<-dormancy_calc_rare(bean_table, 5, 30000)
test0531<-dormancy_calc_rare(bean_table, 5, 31000)
test0532<-dormancy_calc_rare(bean_table, 5, 32000)
test0533<-dormancy_calc_rare(bean_table, 5, 33000)
test0534<-dormancy_calc_rare(bean_table, 5, 34000)
test0535<-dormancy_calc_rare(bean_table, 5, 35000)
test0536<-dormancy_calc_rare(bean_table, 5, 36000)
test0537<-dormancy_calc_rare(bean_table, 5, 37000)
test05372<-dormancy_calc_rare(bean_table, 5, 37123)


#combine data
testy<-rbind(test,test0, test1, test2, test3, test4, test6, test7, test8, test9, test10, test11, test12, test13, test14, test15, test16, test17, test18, test19, test20, test21, test22, test23, test24, test25, test26, test27, test28, test29, test30, test21, test32, test33, test34, test35, test36, test37, test372)
testy2<-rbind(test02,test020, test021, test022, test023, test024, test026, test027, test028, test029, test0210, test0211, test0212, test0213, test0214, test0215, test0216, test0217, test0218, test0219, test0220, test0221, test0222, test0223, test0224, test0225, test0226, test0227, test0228, test0229, test0230, test0221, test0232, test0233, test0234, test0235, test0236, test0237, test02372)
testy5<-rbind(test05,test050, test051, test052, test053, test054, test056, test057, test058, test059, test0510, test0511, test0512, test0513, test0514, test0515, test0516, test0517, test0518, test0519, test0520, test0521, test0522, test0523, test0524, test0525, test0526, test0527, test0528, test0529, test0530, test0521, test0532, test0533, test0534, test0535, test0536, test0537, test05372)
testy_full<-rbind(testy, testy2, testy5)

#add treatment to the table of data 
testy_full$treatment<-row.names(testy_full)
testy_full$treatment<-gsub('cDNA', '', testy_full$treatment)
testy_full$treatment<-gsub('DNA', '', testy_full$treatment)
testy_full$treatment<-gsub('[[:digit:]]+', '', testy_full$treatment)
testy_full$treatment<-gsub('N', 'Nutrients', testy_full$treatment)
testy_full$treatment<-gsub('C', 'Control', testy_full$treatment)
testy_full$treatment<-gsub('D', 'Drought', testy_full$treatment)

#change ratio # to '16S ratio=#'
testy_full$ratio_threshold<-gsub('1', '16S ratio = 1', testy_full$ratio_threshold)
testy_full$ratio_threshold<-gsub('2', '16S ratio = 2', testy_full$ratio_threshold)
testy_full$ratio_threshold<-gsub('5', '16S ratio = 5', testy_full$ratio_threshold)

#plot it
library(ggplot2)
ggplot(testy_full, aes(samp.depth, per.active, colour=treatment))+
  geom_smooth()+
  ylab("Percent Active")+
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))+
  xlab("Sequencing Depth")+
  theme_bw()+
  facet_wrap(~ratio_threshold)


```
