# Syto24 cell counts per gram soil
```
stress_dorm_ctc16 <- read.delim("~/stress_dorm_ctc16.txt")
cyto_counts<-stress_dorm_ctc16[,c(2,5)]

library(ggplot2)
ggplot(cyto_counts, aes(Type, Syto_Counts_per_gram))+
  geom_boxplot()+
  theme_bw()+
  ylab("Cells g-soil-1")+
  xlab("")
  
#test for variance difference
bartlett.test(Syto_Counts_per_gram ~ Type, data=cyto_counts)
#Bartlett's K-squared =45.266, df = 2, p-value =0.1482

#ANOVA
syto_aov<-aov(Syto_Counts_per_gram ~ Type, data=cyto_counts)
summary(syto_aov)

#Df    Sum Sq
#Type          2 1.811e+12
#Residuals   141 1.539e+14
#Mean Sq F value
#Type        9.053e+11   0.829
#Residuals   1.092e+12        
#Pr(>F)
#Type         0.438
#Residuals   
```
