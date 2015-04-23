# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

setwd("~/Desktop/diseaseScript_Final/deseq2/")
library("DESeq2")

#------Read in counts and simplify sample names
countData <- read.table("25july_diseasecounts.tab")
head(countData)
names(countData)=gsub("V3BC_1[3-8][0-9]_[ABCDEFGH][1-6]_", "", names(countData))
names(countData)=sub("_", "", names(countData))

#------Sum technical replicates
names(countData)
countData$AL8 <- rowSums(countData[c(1,9)])
countData$AL7 <- rowSums(countData[c(2,10)])
countData$AL6 <- rowSums(countData[c(3,11)])
countData$AL5 <- rowSums(countData[c(4,12)])
countData$AL4 <- rowSums(countData[c(5,13)])
countData$AL3 <- rowSums(countData[c(6,14)])
countData$AL2 <- rowSums(countData[c(7,15)])
countData$AL1 <- rowSums(countData[c(8,16)])
countData$H16 <- rowSums(countData[c(17,25)])
countData$H15 <- rowSums(countData[c(18,26)])
countData$H14 <- rowSums(countData[c(19,27)])
countData$H13 <- rowSums(countData[c(20,28)])
countData$H12 <- rowSums(countData[c(21,29)])
countData$H11 <- rowSums(countData[c(22,30)])
countData$H10 <- rowSums(countData[c(23,31)])
countData$H9 <- rowSums(countData[c(24,32)])
countData$D8 <- rowSums(countData[c(33,41)])
countData$D7 <- rowSums(countData[c(34,42)])
countData$D6 <- rowSums(countData[c(35,43)])
countData$D5 <- rowSums(countData[c(36,44)])
countData$D4 <- rowSums(countData[c(37,45)])
countData$D3 <- rowSums(countData[c(38,46)])
countData$D2 <- rowSums(countData[c(39,47)])
countData$D1 <- rowSums(countData[c(40,48)])

names(countData)

countData=countData[c(49:72)]
head(countData)

write.csv(countData, "countData_summedTechReps.csv", quote=F)

totalCounts=colSums(countData)
totalCounts

#AL8    AL7    AL6    AL5    AL4    AL3    AL2    AL1    H16    H15    H14    H13    H12    H11    H10     #H9 
#140453 268841 211257 241892  95056 127526 126173  86817 187561 245394 106807 273151 115532 131165 #137994  76800 
#D8     D7     D6     D5     D4     D3     D2     D1 
#154829 146660 125608 191101 126408 149862 190247 102477 

min(totalCounts) #76800
max(totalCounts) #273151
mean(totalCounts) #156650

#------Construct conditions
disease=c("AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", "H", "H", "H", "H", "H", "H", "H", "H", "D", "D", "D", "D", "D", "D", "D", "D") 
indiv= c("8", "7", "6", "5", "4", "3", "2", "1", "16", "15", "14", "13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1")

g=data.frame(disease, indiv)
g
colData<- g

#------Make big dataframe
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ disease) 

#------Sample outlier detection with arrayQualityMetrics
library(arrayQualityMetrics)
library(Biobase)

vsd=varianceStabilizingTransformation(dds, blind=T)
e=ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("disease"), force=T)

#-----No sample outliers detected

######################   ######################    ######################
#                             DESeq
######################   ######################    ######################

dds<-DESeq(dds)

#------See results
res<-results(dds)
head(res)
summary(res)

mcols(res,use.names=TRUE)

#------Histogram of pvalues
quartz()
hist(res$pvalue)

#------Table NAs
table(is.na(res$padj))

#------Table FDR<10%
table(res$padj<0.1)

#------Table Unadjusted pval < 0.05
table(res$pvalue<0.05)

#------Independent filtering
attr(res, "filterThreshold")

#------Pvalues by normalized counts
quartz()
plot(res$baseMean+1, -log2(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[2](pvalue)),
     ylim=c(0,15),
     cex=.4, col=rgb(0,0,0,.3))
    abline(h=-log2(0.05), col="red")

######################   ######################    ######################
#         Use contrasts to compare levels of disease
######################   ######################    ######################

#------AL to healthy
resAH=results(dds, contrast=c("disease", "AL", "H"))
table(resAH$padj<0.1)
table(resAH$pvalue<0.05)

#------Disease to healthy
resDH=results(dds, contrast=c("disease","D","H"))
table(resDH$padj<0.1)
table(resDH$pvalue<0.05)

#------Disease to AL
resDA=results(dds, contrast=c("disease","D","AL"))
table(resDA$padj<0.1)
table(resDA$pvalue<0.05)

save(dds, res, resAH, resDH, resDA, file="DESeq2Results.Rdata")


######################   ######################    ######################
#                         Visualizations
######################   ######################    ######################

library("RColorBrewer");
library("gplots");
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100);

#------MA plot
quartz()
plotMA(dds,ylim=c(-2,2),main="DESeq2");

#------Dispersion estimates
quartz()
plotDispEsts(dds);

######################   ######################    ######################
#                         Exporting
######################   ######################    ######################

#######-----------get VSD
vsd=getVarianceStabilizedData(dds) 
head(vsd)
colnames(vsd)=paste(g$disease, g$indiv, sep="")

###--------------get pvals from contrasts
head(resDH)
valsDH=cbind(resDH$pvalue, resDH$padj)
head(valsDH)
colnames(valsDH)=c("pvalDH", "padjDH")

head(resAH)
valsAH=cbind(resAH$pvalue, resAH$padj)
head(valsAH)
colnames(valsAH)=c("pvalAH", "padjAH")

head(resDA)
valsDA=cbind(resDA$pvalue, resDA$padj)
head(valsDA)
colnames(valsDA)=c("pvalDA", "padjDA")

results=cbind(valsDH,valsAH,valsDA)
head(results)

######-------------make vsd and pvals table
vsdpvals=cbind(vsd,results)
head(vsdpvals)

write.csv(vsdpvals, "VSDandPVALS_disease.csv", quote=F)

####---------------get data for WGCNA
head(vsdpvals)
vsdp=as.data.frame(vsdpvals)

AH01=vsdp[vsdp$pvalAH<0.1 & !is.na(vsdp$pvalAH),]
length(AH01[,1]) 

DA01=vsdp[vsdp$pvalDA<0.1 & !is.na(vsdp$pvalDA),]
length(DA01[,1]) 

DH01=vsdp[vsdp$pvalDH<0.1 & !is.na(vsdp$pvalDH),]
length(DH01[,1]) 

degs01=union(row.names(AH01),row.names(DA01))
degs01=union(degs01,row.names(DH01))
length(degs01)

wdegs01=vsd[(row.names(vsd) %in% degs01),]
head(wdegs01)

write.csv(wdegs01, "genes4WGCNA.csv", quote=F)


######----------log GO
head(resDA)
logs=data.frame(cbind("gene"=row.names(resDA),"logP"=round(-log(resDA$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resDA$log2FoldChange<0]=-1  ##change to correct model
table(sign)
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_DA_logP.csv",sep=",")

head(resAH)
logs=data.frame(cbind("gene"=row.names(resAH),"logP"=round(-log(resAH$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resAH$log2FoldChange<0]=-1  ##change to correct model
table(sign)
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_AH_logP.csv",sep=",")

head(resDH)
logs=data.frame(cbind("gene"=row.names(resDH),"logP"=round(-log(resDH$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resDH$log2FoldChange<0]=-1  ##change to correct model
table(sign)
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_DH_logP.csv",sep=",")

#-------------write results tables; includes log2fc and pvals
write.table(results(dds), file="DESeq.results.txt", quote=FALSE, sep="\t");  
write.table(resDA, file="DESeq.results.DA.txt", quote=F, sep="\t")
write.table(resAH, file="DESeq.results.AH.txt", quote=F, sep="\t")
write.table(resDH, file="DESeq.results.DH.txt", quote=F, sep="\t")

#-------------write annotated results tables
gg=read.delim("~/Documents/genomes/ahya_annotations_may25_2014/ahya2digNvec_plus_iso2gene.tab", sep="\t")
head(gg)

resDA=as.data.frame(resDA)
resDA$X=row.names(resDA)
names(resDA)
resDA=resDA[c(7,1:6)]
head(resDA)
resDAannot=merge(resDA,gg,by=1)
head(resDAannot)
resDAannot <- resDAannot[order(resDAannot$padj),]
write.table(resDAannot, "annotated_results_DA.txt", sep="\t", quote=F, row.names=F)

resDH=as.data.frame(resDH)
resDH$X=row.names(resDH)
names(resDH)
resDH=resDH[c(7,1:6)]
resDHannot=merge(resDH,gg,by=1)
head(resDHannot)
resDHannot <- resDHannot[order(resDHannot$padj),]
write.table(resDHannot, "annotated_results_DH.txt", sep="\t", quote=F, row.names=F)

resAH=as.data.frame(resAH)
resAH$X=row.names(resAH)
names(resAH)
resAH=resAH[c(7,1:6)]
resAHannot=merge(resAH,gg,by=1)
head(resAHannot)
resAHannot <- resAHannot[order(resAHannot$padj),]
write.table(resAHannot, "annotated_results_AH.txt", sep="\t", quote=F, row.names=F)
