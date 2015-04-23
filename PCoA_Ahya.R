#-------------------
# Principal Coordinate Analysis
#-------------------

library(vegan)
library(rgl)
library(ape)

setwd("~/Desktop/diseaseScript_Final/pcoa/")
dat=read.csv("VSDandPVALS_disease.csv")
names(dat)
head(dat)

data=dat[,2:25]
row.names(data)=data$X
head(data)

#-------------Set experimental conditions
names(data)
disease=c(rep("AL",8), rep("H",8), rep("D",8))
indiv= c("8", "7", "6", "5", "4", "3", "2", "1", "16", "15", "14", "13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1")
conditions=data.frame(cbind(disease,indiv))
conditions

#-------------Calulate principal coordinates
dd.veg=vegdist(t(data), "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

#-------------First and second axes
quartz()
plot(scores[,1], scores[,2],col=as.numeric(conditions$disease), pch=16)
ordihull(scores,disease,label=T)

#-------------Second and third axes
quartz()
plot(scores[,2], scores[,3],col=as.numeric(conditions$disease),pch=16)
ordihull(scores[,2:3],disease,label=T)

#-------------PERMANOVA
adonis(t(data)~disease+indiv,data=conditions,method="manhattan")  
adonis(t(data)~disease,data=conditions,method="manhattan")  

pco1=scores[,1]
TukeyHSD(aov(pco1~disease))

##########-----Unadjustd pvalue < 0.05 Only------##############

#-------------Subset for significant
head(dat)
DHdat=row.names(dat[dat$pvalDH<0.05 & !is.na(dat$pvalDH),])
AHdat=row.names(dat[dat$pvalAH<0.05 & !is.na(dat$pvalAH),])
DAdat=row.names(dat[dat$pvalDA<0.05 & !is.na(dat$pvalDA),])

sdat=union(DHdat,AHdat)
sdat=union(sdat,DAdat)
length(sdat)

sdata=dat[(row.names(dat) %in% sdat),]
head(sdata)
names(sdata)

data=sdata[,2:25]
row.names(data)=sdata$X
head(data)

#-------------Calulate principal coordinates
dd.veg=vegdist(t(data), "manhattan")
div.dd.veg=dd.veg/1000
head(div.dd.veg)

dd.pcoa=pcoa(div.dd.veg) 
head(dd.pcoa)
scores=dd.pcoa$vectors

#-------------First and second axes
quartz()
plot(scores[,1], scores[,2],col=as.numeric(conditions$disease), pch=16)
ordihull(scores,disease,label=T)

#-------------Second and third axes 
quartz()
plot(scores[,2], scores[,3],col=as.numeric(conditions$disease),pch=16)
ordihull(scores[,2:3],disease,label=T)


#-------------PERMANOVA
adonis(t(data)~disease+indiv,data=conditions,method="manhattan")  
adonis(t(data)~disease,data=conditions,method="manhattan")  

pco1=scores[,1]
TukeyHSD(aov(pco1~disease))
