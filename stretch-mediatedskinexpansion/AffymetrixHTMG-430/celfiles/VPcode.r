#read data and pre-processing
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
library(GEOquery)
gset <- getGEO("GSE126231", destdir = ".")
names(gset)
gset_one <- gset[[1]]
class(gset_one)
##[1] "ExpressionSet"
##attr(,"package")
##[1] "Biobase"
exprSet <- exprs(gset_one)
dim(exprSet)
##[1] 45141     8
pdata <- pData(gset_one)
dim(pdata)
##[1]  8 33
download.file(pdata$supplementary_file[1],strsplit(pdata$supplementary_file[1],split='/')[[1]][9])
download.file(pdata$supplementary_file[2],strsplit(pdata$supplementary_file[2],split='/')[[1]][9])
download.file(pdata$supplementary_file[3],strsplit(pdata$supplementary_file[3],split='/')[[1]][9])
download.file(pdata$supplementary_file[4],strsplit(pdata$supplementary_file[4],split='/')[[1]][9])
download.file(pdata$supplementary_file[5],strsplit(pdata$supplementary_file[5],split='/')[[1]][9])
download.file(pdata$supplementary_file[6],strsplit(pdata$supplementary_file[6],split='/')[[1]][9])
download.file(pdata$supplementary_file[7],strsplit(pdata$supplementary_file[7],split='/')[[1]][9])
download.file(pdata$supplementary_file[8],strsplit(pdata$supplementary_file[8],split='/')[[1]][9])
library(GEOquery)
#rsrc <- NetAffxResource("zliu39-c@my.cityu.edu.hk", "13661385845a")
#grep("mg-430", names(rsrc), ignore.case = TRUE, value = TRUE)
#annos <- rsrc[["HT_MG-430_PM"]]
#sapply(affxAnnotation(annos), force)[c(4,11)]
#readAnnotation(rsrc, annotation = affxAnnotation(annos)[[4]], content = FALSE)
#readAnnotation(rsrc, annotation = affxAnnotation(annos)[[11]], content = FALSE)
#unzip("HT_MG-430_PM.probe_tab.zip")
#unzip("/HT_MG-430_PM.cdf.zip")
#gunzip("GSM3594010_661_2017.ga.cel.gz")
#seed <- new("AffyExpressionPDInfoPkgSeed", cdfFile = "HT_MG-430_PM.cdf", celFile = "GSM3594010_661_2017.ga.cel", tabSeqFile = "HT_MG-430_PM.probe_tab")
#makePdInfoPackage(seed)
#install.packages("pd.ht.mg.430.pm/", repos = NULL, type = "source")
library("pd.ht.mg.430a")
file_CELs <- list.celfiles("../celfiles", listGzipped = TRUE, full.name = TRUE)
rawAffy <- read.celfiles(filenames = file_CELs, phenoData = phenoData(gset_one), sampleNames = rownames(pData(gset_one)))
head(exprs(rawAffy), n = 2)
#芯片常用的标准化方法也有好几种（如：MAS5），这里就不一一尝试了，就拿最常用的RMA（Robust Multichip Average algorithm），引用芯片教程，RMA标准化过程主要分为3步：
#1. Background correction (removes array auto-fluorescence ) 背景矫正
#2. Quantile normalization (makes all intensity distributions identical) 标准化
#3. Probeset summarization (calculates one representative value per probeset) 汇总
class(rawAffy)
##[1] "ExpressionFeatureSet"
##attr(,"package")
##[1] "oligoClasses"
eset <- rma(rawAffy)
exprSet <- exprs(eset)
dim(exprSet)
#[1] 45141     8
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430")
#install.packages("htmg430pmmmentrezg.db_13.0.1.tar.gz", repos = NULL, type="source")
library("htmg430pmmmentrezg.db")
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
library(Biobase)
gpl <- getGEO('GPL11180', destdir=".")
colnames(Table(gpl)) 
head(Table(gpl)[,c(1,11)]) ## you need to check this , which column do you need
write.csv(Table(gpl)[,c(1,11)],"GPL11180.csv")
probeset <- rownames(exprSet)
probe2symbol=Table(gpl)[,c(1,11)]
mydata<-cbind(exprSet,row.names(exprSet))
colnames(mydata)<-c("CTRL_1","CTRL_2","CTRL_3","EXPD4_1","EXPD4_2","EXPD4_3","TPA_1","TPA_2","ID")
mydata<-merge(mydata,probe2symbol,order="ID")
head(mydata)
dim(mydata)
write.csv(mydata,"GSE126231.csv",row.names=F)
#use mean
data<-read.csv("GSE126231.csv",header=T,row.names=1)
exprSet_symbol1 <- aggregate(x = data[,1:(ncol(data)-1)],
                             by = list(data$Gene.Symbol),
                             FUN = mean)
head(exprSet_symbol1)
exprSet_symbol2 <-  as.data.frame(t(sapply(split(data,data$Gene.Symbol),
                                           function(x) colMeans(x[,1:(ncol(x)-1)]))))
head(exprSet_symbol2)
write.csv(exprSet_symbol2,"GSE126231_uniq.csv")

exprSet_symbol2<-read.table("GSE126231_uniq.csv",sep=",",header=T,row.names = 1)
head(exprSet_symbol2)
emat.rma.log2 <- exprSet
emat.rma.nologs <- 2^emat.rma.log2
class(emat.rma.log2)
head(emat.rma.log2, 1)
emat.rma.nologs <- 2^emat.rma.log2
head(emat.rma.nologs, 1)
#计算平均值，并做对数转换
results.rma <- data.frame((emat.rma.log2[,c(1,4)] + emat.rma.log2[,c(2,5)] + emat.rma.log2[,c(3,6)])/3)
results.rma$GSM3594016 <- (emat.rma.log2[,7] + emat.rma.log2[,8])/2
results.rma$fc.EXPD4.CTRL <- results.rma[,2]-results.rma[,1]
results.rma$fc.TPA.CTRL <- results.rma[,3]-results.rma[,1]
head(results.rma)
subset.logic <- results.rma$fc.EXPD4.CTRL>0
subset.data <- results.rma[subset.logic,]
length(subset.logic); nrow(results.rma)
head(subset.logic)
apply(abs(results.rma[,4:5]), 2, max)
p.value <- apply(emat.rma.log2, 1, function(x){t.test(x[1:3], x[4:6])$p.value})
results.rma$EXPD4.p.value <- p.value
names(results.rma)
results.rma <- results.rma[, c(1,2,4,6)]
probeset <- rownames(results.rma)
library("htmg430pmmmentrezg.db")
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
library(Biobase)
gpl <- getGEO('GPL11180', destdir=".")
probe2symbol=Table(gpl)[,c(1,11)]
mydata<-cbind(results.rma,row.names(results.rma))
colnames(mydata)<-c("CTRL","EXPD4","FC","Pvalue","ID")
mydata<-merge(mydata,probe2symbol,order="ID")
head(mydata)
dim(mydata)
write.csv(mydata,"allgenes.EXPD4vsCTRL.name.csv")
mydata[which(mydata$Pvalue %in% NA),'sig'] <- 'no diff'
mydata[which(mydata$FC >= 1 & mydata$Pvalue < 0.05),'sig'] <- 'up (p.adj < 0.05, log2FC >= 1)'
mydata[which(mydata$FC <= -1 & mydata$Pvalue < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -1)'
mydata[which(abs(mydata$FC) < 1 | mydata$Pvalue >= 0.05),'sig'] <- 'no diff'
library(ggplot2)
pdf("volcano_p.EXPD4.pdf")
volcano_p <- ggplot(mydata, aes(FC, -log(Pvalue, 10))) +
  geom_point(aes(color = sig),alpha = 0.6, size = 1) +
  scale_color_manual(values = c('blue2', 'gray30', 'red2')) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.position = c(0.26, 0.92)) +
  theme(legend.title = element_blank(), 
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value', color = NA) +
  xlim(-8, 8)+ylim(0, 5)
volcano_p
dev.off()

