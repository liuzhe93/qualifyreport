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
#[1] 45141    10
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
sum(abs(results.rma[,"fc.EXPD4.CTRL"])>=log2(2))
#4069
results.st <- results.rma[abs(results.rma$fc.EXPD4.CTRL)>=log2(2),]
sel.genes <- row.names(results.st)
#t测验，并选出p<0.05的差异表达基因：
p.value <- apply(emat.rma.log2[sel.genes,], 1, function(x){t.test(x[1:3], x[4:6])$p.value})
results.st$p.value <- p.value
names(results.st)
results.st <- results.st[, c(1,2,4,6)]
results.st <- results.st[p.value<0.05,]
head(results.st, 2)
nrow(results.st)
#1037
write.csv(results.st,"DEGs.EXPD4vsCTRL.csv")
library("pheatmap")
heatmap_data<-emat.rma.log2[row.names(results.st),1:6]
colnames(heatmap_data)<-c("CTRL_rep1","CTRL_rep2","CTRL_rep3","EXPD4_rep1","EXPD4_rep2","EXPD4_rep3")
pheatmap(heatmap_data, scale = "row", show_rownames = F)


probeset <- rownames(results.st)
library("htmg430pmmmentrezg.db")
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
library(Biobase)
gpl <- getGEO('GPL11180', destdir=".")
probe2symbol=Table(gpl)[,c(1,11)]
mydata<-cbind(results.st,row.names(results.st))
colnames(mydata)<-c("CTRL","EXPD4","FC","Pvalue","ID")
mydata<-merge(mydata,probe2symbol,order="ID")
head(mydata)
dim(mydata)
#[1] 1037    6
write.csv(mydata,"DEGs.EXPD4vsCTRL.name.csv")


library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
keytypes(org.Mm.eg.db) 
library("biomaRt")
library("curl")
mart<-useDataset("mmusculus_gene_ensembl",useMart("ensembl"))
x_bg <- read.table("background.gene.txt",header=FALSE)
x_bg$V1 <- as.character(x_bg$V1) #需要character格式，然后进行ID转化
bg = bitr(x_bg$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
head(bg,2)
x_cluster1 <- read.csv("DEGs.EXPD4vsCTRL.name.csv",header=T,row.names=1) #单列基因名文件
x_cluster1[which(x_cluster1$Pvalue %in% NA),'sig'] <- 'no diff'
x_cluster1[which(x_cluster1$FC >= 1 & x_cluster1$Pvalue < 0.05),'sig'] <- 'up (p.adj < 0.05, log2FC >= 1)'
x_cluster1[which(x_cluster1$FC <= -1 & x_cluster1$Pvalue < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -1)'
x_cluster1[which(abs(x_cluster1$FC) < 1 | x_cluster1$Pvalue >= 0.05),'sig'] <- 'no diff'

x<-x_cluster1$Gene.Symbol
test1 = bitr(x, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",level = 3,readable = TRUE)
ego_all <- enrichGO(gene = test1$ENTREZID, universe = bg$ENTREZID,OrgDb = org.Mm.eg.db,ont = "All", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_all1 <- setReadable(ego_all, OrgDb = org.Mm.eg.db)
write.csv(summary(ego_all),"ALL-enrich_EXPD4.csv",row.names =FALSE)
# BP
ego_BP <- enrichGO(gene = test1$ENTREZID, universe = bg$ENTREZID,OrgDb = org.Mm.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mm.eg.db)
pdf("EXPD4.pdf")
clusterProfiler::dotplot(ego_BP,showCategory=20,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
barplot(ego_BP, showCategory=20,title="EnrichmentGO_BP")#条状图，按p从小到大排，绘制前20个Term
plotGOgraph(ego_BP)
dev.off()
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'mmu', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)
pdf("KEGG.pdf")
dotplot(kk,title="Enrichment KEGG_dot")
dev.off()
x1<-cbind(x_cluster1$Gene.Symbol,x_cluster1$FC)
colnames(x1)<-c("SYMBOL","log2FoldChange")
gene_list<-merge(x1,test1,order="SYMBOL")
temp<-subset(gene_list,select=c("ENTREZID","log2FoldChange"))
write.table(temp,"geneList.EXPD4.txt",sep="\t",quote=F,row.names=F)
temp<-read.table("geneList.EXPD4.txt",sep="\t",header=T)
library(dplyr)
temp <- temp %>%
  group_by(ENTREZID) %>%
  summarise_all(mean)
value<-t(t(temp$log2FoldChange))
rownames(value)<-as.character(temp$ENTREZID)
mmu04141 <- pathview(gene.data = value,
                     pathway.id = "mmu04141", #上述结果中的mmu04141通路
                     species = "mmu")

x_up <- x_cluster1[x_cluster1$FC>0,]
dim(x_up)
#[1] 933   7
x_down <- x_cluster1[x_cluster1$FC<0,]
dim(x_down)
#[1] 104   7
x_up_genes<-x_up$Gene.Symbol
test1 = bitr(x_up_genes, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",level = 3,readable = TRUE)
ego_all <- enrichGO(gene = test1$ENTREZID, universe = bg$ENTREZID,OrgDb = org.Mm.eg.db,ont = "All", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_all1 <- setReadable(ego_all, OrgDb = org.Mm.eg.db)
write.csv(summary(ego_all),"ALL-enrich_EXPD4.up.csv",row.names =FALSE)
# BP
ego_BP <- enrichGO(gene = test1$ENTREZID, universe = bg$ENTREZID,OrgDb = org.Mm.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mm.eg.db)
write.csv(ego_BP1,"BP.enrich_EXPD4.up.csv",row.names =FALSE)
top20BP<-head(ego_BP1,20)
top20genes<-top20BP$geneID
g1<-strsplit(top20genes[1],split='/')
g2<-strsplit(top20genes[2],split='/')
g3<-strsplit(top20genes[3],split='/')
g4<-strsplit(top20genes[4],split='/')
g5<-strsplit(top20genes[5],split='/')
g6<-strsplit(top20genes[6],split='/')
g7<-strsplit(top20genes[7],split='/')
g8<-strsplit(top20genes[8],split='/')
g9<-strsplit(top20genes[9],split='/')
g10<-strsplit(top20genes[10],split='/')
g11<-strsplit(top20genes[11],split='/')
g12<-strsplit(top20genes[12],split='/')
g13<-strsplit(top20genes[13],split='/')
g14<-strsplit(top20genes[14],split='/')
g15<-strsplit(top20genes[15],split='/')
g16<-strsplit(top20genes[16],split='/')
g17<-strsplit(top20genes[17],split='/')
g18<-strsplit(top20genes[18],split='/')
g19<-strsplit(top20genes[19],split='/')
g20<-strsplit(top20genes[20],split='/')
g1<-as.data.frame(g1)
g2<-as.data.frame(g2)
g3<-as.data.frame(g3)
g4<-as.data.frame(g4)
g5<-as.data.frame(g5)
g6<-as.data.frame(g6)
g7<-as.data.frame(g7)
g8<-as.data.frame(g8)
g9<-as.data.frame(g9)
g10<-as.data.frame(g10)
g11<-as.data.frame(g11)
g12<-as.data.frame(g12)
g13<-as.data.frame(g13)
g14<-as.data.frame(g14)
g15<-as.data.frame(g15)
g16<-as.data.frame(g16)
g17<-as.data.frame(g17)
g18<-as.data.frame(g18)
g19<-as.data.frame(g19)
g20<-as.data.frame(g20)
colnames(g1)<-"geneName"
colnames(g2)<-"geneName"
colnames(g3)<-"geneName"
colnames(g4)<-"geneName"
colnames(g5)<-"geneName"
colnames(g6)<-"geneName"
colnames(g7)<-"geneName"
colnames(g8)<-"geneName"
colnames(g9)<-"geneName"
colnames(g10)<-"geneName"
colnames(g11)<-"geneName"
colnames(g12)<-"geneName"
colnames(g13)<-"geneName"
colnames(g14)<-"geneName"
colnames(g15)<-"geneName"
colnames(g16)<-"geneName"
colnames(g17)<-"geneName"
colnames(g18)<-"geneName"
colnames(g19)<-"geneName"
colnames(g20)<-"geneName"
merge_genes<-rbind(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20)
table(merge_genes)
counts<-table(merge_genes)


x_cluster1[x_cluster1$Gene.Symbol=="Akt1",]


pdf("EXPD4.up.BP.pdf")
clusterProfiler::dotplot(ego_BP,showCategory=20,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
barplot(ego_BP, showCategory=20,title="EnrichmentGO_BP")#条状图，按p从小到大排，绘制前20个Term
plotGOgraph(ego_BP)
dev.off()
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'mmu', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
write.csv(summary(kk),"KEGG-enrich.up.csv",row.names =FALSE)
pdf("KEGG.up.pdf")
dotplot(kk,title="Enrichment KEGG_dot")
dev.off()
x1<-cbind(x_up$Gene.Symbol,x_up$FC)
colnames(x1)<-c("SYMBOL","log2FoldChange")
gene_list<-merge(x1,test1,order="SYMBOL")
temp<-subset(gene_list,select=c("ENTREZID","log2FoldChange"))
write.table(temp,"geneList.EXPD4.up.txt",sep="\t",quote=F,row.names=F)


x_down_genes<-x_down$Gene.Symbol
test1 = bitr(x_down_genes, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",level = 3,readable = TRUE)
ego_all <- enrichGO(gene = test1$ENTREZID, universe = bg$ENTREZID,OrgDb = org.Mm.eg.db,ont = "All", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_all1 <- setReadable(ego_all, OrgDb = org.Mm.eg.db)
write.csv(summary(ego_all),"ALL-enrich_EXPD4.down.csv",row.names =FALSE)
# BP
ego_BP <- enrichGO(gene = test1$ENTREZID, universe = bg$ENTREZID,OrgDb = org.Mm.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mm.eg.db)
pdf("EXPD4.down.BP.pdf")
clusterProfiler::dotplot(ego_BP,showCategory=20,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
barplot(ego_BP, showCategory=20,title="EnrichmentGO_BP")#条状图，按p从小到大排，绘制前20个Term
plotGOgraph(ego_BP)
dev.off()
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'mmu', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
write.csv(summary(kk),"KEGG-enrich.down.csv",row.names =FALSE)
pdf("KEGG.down.pdf")
dotplot(kk,title="Enrichment KEGG_dot")
dev.off()
x1<-cbind(x_down$Gene.Symbol,x_down$FC)
colnames(x1)<-c("SYMBOL","log2FoldChange")
gene_list<-merge(x1,test1,order="SYMBOL")
temp<-subset(gene_list,select=c("ENTREZID","log2FoldChange"))
write.table(temp,"geneList.EXPD4.down.txt",sep="\t",quote=F,row.names=F)


























# adhesion
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
adhesion<-read.table("adhesion.txt")
colnames(adhesion)<-"X"
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
exprSet_symbol2<-read.csv("GSE126231_uniq.csv")
expdata<-exprSet_symbol2[,2:9]
row.names(expdata)<-exprSet_symbol2[,1]
adhesion<-merge(adhesion,exprSet_symbol2,order=X)
library("ggplot2")
library("ggsignif")
## adhesion
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
ADH<-t(adhesion)
ADH<-ADH[2:dim(ADH)[1],]
ADH<-cbind(ADH,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(ADH)<-c(adhesion$X,"group")
ADH<-as.data.frame(ADH)
ADH<-ADH[1:6,]
compaired <- list(c("Control", "Expandation"))
pdf("adhesion.Itga3.pdf",width=8, height=8)
ADH$Itga3_val<-as.double(ADH$Itga3)
ggplot(ADH,aes(x=group,y=Itga3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Cadm1.pdf",width=8, height=8)
ADH$Cadm1_val<-as.double(ADH$Cadm1)
ggplot(ADH,aes(x=group,y=Cadm1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Itga6.pdf",width=8, height=8)
ADH$Itga6_val<-as.double(ADH$Itga6)
ggplot(ADH,aes(x=group,y=Itga6_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Itgb1.pdf",width=8, height=8)
ADH$Itgb1_val<-as.double(ADH$Itgb1)
ggplot(ADH,aes(x=group,y=Itgb1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Lad1.pdf",width=8, height=8)
ADH$Lad1_val<-as.double(ADH$Lad1)
ggplot(ADH,aes(x=group,y=Lad1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Lamc2.pdf",width=8, height=8)
ADH$Lamc2_val<-as.double(ADH$Lamc2)
ggplot(ADH,aes(x=group,y=Lamc2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Parva.pdf",width=8, height=8)
ADH$Parva_val<-as.double(ADH$Parva)
ggplot(ADH,aes(x=group,y=Parva_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Pxn.pdf",width=8, height=8)
ADH$Pxn_val<-as.double(ADH$Pxn)
ggplot(ADH,aes(x=group,y=Pxn_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("adhesion.Tln1.pdf",width=8, height=8)
ADH$Tln1_val<-as.double(ADH$Tln1)
ggplot(ADH,aes(x=group,y=Tln1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()

## cytoskeleton

cytoskeleton<-read.table("Cytoskeleton.txt")
colnames(cytoskeleton)<-"X"
cytoskeleton<-merge(cytoskeleton,exprSet_symbol2,order=X)
CYT<-t(cytoskeleton)
CYT<-CYT[2:dim(CYT)[1],]
CYT<-cbind(CYT,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(CYT)<-c(cytoskeleton$X,"group")
CYT<-as.data.frame(CYT)
CYT<-CYT[1:6,]
compaired <- list(c("Control", "Expandation"))
pdf("cytoskeleton.Actn1.pdf",width=8, height=8)
CYT$Actn1_val<-as.double(CYT$Actn1)
ggplot(CYT,aes(x=group,y=Actn1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Actr3.pdf",width=8, height=8)
CYT$Actr3_val<-as.double(CYT$Actr3)
ggplot(CYT,aes(x=group,y=Actr3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Arpc4.pdf",width=8, height=8)
CYT$Arpc4_val<-as.double(CYT$Arpc4)
ggplot(CYT,aes(x=group,y=Arpc4_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Cnbp.pdf",width=8, height=8)
CYT$Cnbp_val<-as.double(CYT$Cnbp)
ggplot(CYT,aes(x=group,y=Cnbp_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Flna.pdf",width=8, height=8)
CYT$Flna_val<-as.double(CYT$Flna)
ggplot(CYT,aes(x=group,y=Flna_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Myh9.pdf",width=8, height=8)
CYT$Myh9_val<-as.double(CYT$Myh9)
ggplot(CYT,aes(x=group,y=Myh9_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Myo1b.pdf",width=8, height=8)
CYT$Myo1b_val<-as.double(CYT$Myo1b)
ggplot(CYT,aes(x=group,y=Myo1b_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Myo5a.pdf",width=8, height=8)
CYT$Myo5a_val<-as.double(CYT$Myo5a)
ggplot(CYT,aes(x=group,y=Myo5a_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Nav1.pdf",width=8, height=8)
CYT$Nav1_val<-as.double(CYT$Nav1)
ggplot(CYT,aes(x=group,y=Nav1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Pdlim5.pdf",width=8, height=8)
CYT$Pdlim5_val<-as.double(CYT$Pdlim5)
ggplot(CYT,aes(x=group,y=Pdlim5_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Tmsb10.pdf",width=8, height=8)
CYT$Tmsb10_val<-as.double(CYT$Tmsb10)
ggplot(CYT,aes(x=group,y=Tmsb10_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Tubb3.pdf",width=8, height=8)
CYT$Tubb3_val<-as.double(CYT$Tubb3)
ggplot(CYT,aes(x=group,y=Tubb3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cytoskeleton.Twf1.pdf",width=8, height=8)
CYT$Twf1_val<-as.double(CYT$Twf1)
ggplot(CYT,aes(x=group,y=Twf1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()


## prolioferation.txt
prolioferation<-read.table("prolioferation.txt")
colnames(prolioferation)<-"X"
prolioferation<-merge(prolioferation,exprSet_symbol2,order=X)
PRO<-t(prolioferation)
colnames(PRO)<-PRO[1,]
PRO<-PRO[2:dim(PRO)[1],]
PRO<-cbind(PRO,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
PRO<-as.data.frame(PRO)
colnames(PRO)<-c(colnames(PRO)[1:length(colnames(PRO))-1],"group")
PRO<-PRO[1:6,]
compaired <- list(c("Control", "Expandation"))
pdf("PRO.E2f1.pdf",width=8, height=8)
PRO$E2f1_val<-as.double(PRO$E2f1)
ggplot(PRO,aes(x=group,y=E2f1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Egfr.pdf",width=8, height=8)
PRO$Egfr_val<-as.double(PRO$Egfr)
ggplot(PRO,aes(x=group,y=Egfr_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Egr1.pdf",width=8, height=8)
PRO$Egr1_val<-as.double(PRO$Egr1)
ggplot(PRO,aes(x=group,y=Egr1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Klk10.pdf",width=8, height=8)
PRO$Klk10_val<-as.double(PRO$Klk10)
ggplot(PRO,aes(x=group,y=Klk10_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Krt6a.pdf",width=8, height=8)
PRO$Krt6a_val<-as.double(PRO$Krt6a)
ggplot(PRO,aes(x=group,y=Krt6a_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Map2k3.pdf",width=8, height=8)
PRO$Map2k3_val<-as.double(PRO$Map2k3)
ggplot(PRO,aes(x=group,y=Map2k3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Mapk1.pdf",width=8, height=8)
PRO$Mapk1_val<-as.double(PRO$Mapk1)
ggplot(PRO,aes(x=group,y=Mapk1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Rras2.pdf",width=8, height=8)
PRO$Rras2_val<-as.double(PRO$Rras2)
ggplot(PRO,aes(x=group,y=Rras2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.S100a8.pdf",width=8, height=8)
PRO$S100a8_val<-as.double(PRO$S100a8)
ggplot(PRO,aes(x=group,y=S100a8_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Sprr1a.pdf",width=8, height=8)
PRO$Sprr1a_val<-as.double(PRO$Sprr1a)
ggplot(PRO,aes(x=group,y=Sprr1a_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Mapk1.pdf",width=8, height=8)
PRO$Mapk1_val<-as.double(PRO$Mapk1)
ggplot(PRO,aes(x=group,y=Mapk1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Smad2.pdf",width=8, height=8)
PRO$Smad2_val<-as.double(PRO$Smad2)
ggplot(PRO,aes(x=group,y=Smad2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Smad3.pdf",width=8, height=8)
PRO$Smad3_val<-as.double(PRO$Smad3)
ggplot(PRO,aes(x=group,y=Smad3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Ccnb3.pdf",width=8, height=8)
PRO$Ccnb3_val<-as.double(PRO$Ccnb3)
ggplot(PRO,aes(x=group,y=Ccnb3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Ccnb1.pdf",width=8, height=8)
PRO$Ccnb1_val<-as.double(PRO$Ccnb1)
ggplot(PRO,aes(x=group,y=Ccnb1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Ccnb2.pdf",width=8, height=8)
PRO$Ccnb2_val<-as.double(PRO$Ccnb2)
ggplot(PRO,aes(x=group,y=Ccnb2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Cdc7.pdf",width=8, height=8)
PRO$Cdc7_val<-as.double(PRO$Cdc7)
ggplot(PRO,aes(x=group,y=Cdc7_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Myc.pdf",width=8, height=8)
PRO$Myc_val<-as.double(PRO$Myc)
ggplot(PRO,aes(x=group,y=Myc_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Cav1.pdf",width=8, height=8)
PRO$Cav1_val<-as.double(PRO$Cav1)
ggplot(PRO,aes(x=group,y=Cav1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Igfbp5.pdf",width=8, height=8)
PRO$Igfbp5_val<-as.double(PRO$Igfbp5)
ggplot(PRO,aes(x=group,y=Igfbp5_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Itgb3.pdf",width=8, height=8)
PRO$Itgb3_val<-as.double(PRO$Itgb3)
ggplot(PRO,aes(x=group,y=Itgb3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Fn1.pdf",width=8, height=8)
PRO$Fn1_val<-as.double(PRO$Fn1)
ggplot(PRO,aes(x=group,y=Fn1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Itga2.pdf",width=8, height=8)
PRO$Itga2_val<-as.double(PRO$Itga2)
ggplot(PRO,aes(x=group,y=Itga2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Postn.pdf",width=8, height=8)
PRO$Postn_val<-as.double(PRO$Postn)
ggplot(PRO,aes(x=group,y=Postn_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Fam83d.pdf",width=8, height=8)
PRO$Fam83d_val<-as.double(PRO$Fam83d)
ggplot(PRO,aes(x=group,y=Fam83d_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Kif14.pdf",width=8, height=8)
PRO$Kif14_val<-as.double(PRO$Kif14)
ggplot(PRO,aes(x=group,y=Kif14_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("PRO.Cdk1.pdf",width=8, height=8)
PRO$Cdk1_val<-as.double(PRO$Cdk1)
ggplot(PRO,aes(x=group,y=Cdk1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()



## immune
immune<-read.table("immune.txt")
colnames(immune)<-"X"
immune<-merge(immune,exprSet_symbol2,order=X)
IMM<-t(immune)
colnames(IMM)<-IMM[1,]
IMM<-IMM[2:dim(IMM)[1],]
IMM<-cbind(IMM,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
IMM<-as.data.frame(IMM)
colnames(IMM)<-c(colnames(IMM)[1:length(colnames(IMM))-1],"group")
IMM<-IMM[1:6,]
compaired <- list(c("Control", "Expandation"))
pdf("IMM.Bcap29.pdf",width=8, height=8)
IMM$Bcap29_val<-as.double(IMM$Bcap29)
ggplot(IMM,aes(x=group,y=Bcap29_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Bcl2l14.pdf",width=8, height=8)
IMM$Bcl2l14_val<-as.double(IMM$Bcl2l14)
ggplot(IMM,aes(x=group,y=Bcl2l14_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Ccng1.pdf",width=8, height=8)
IMM$Ccng1_val<-as.double(IMM$Ccng1)
ggplot(IMM,aes(x=group,y=Ccng1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Ifnar1.pdf",width=8, height=8)
IMM$Ifnar1_val<-as.double(IMM$Ifnar1)
ggplot(IMM,aes(x=group,y=Ifnar1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Isg20.pdf",width=8, height=8)
IMM$Isg20_val<-as.double(IMM$Isg20)
ggplot(IMM,aes(x=group,y=Isg20_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tnfaip8l1.pdf",width=8, height=8)
IMM$Tnfaip8l1_val<-as.double(IMM$Tnfaip8l1)
ggplot(IMM,aes(x=group,y=Tnfaip8l1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tnfaip8l2.pdf",width=8, height=8)
IMM$Tnfaip8l2_val<-as.double(IMM$Tnfaip8l2)
ggplot(IMM,aes(x=group,y=Tnfaip8l2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Ifnar2.pdf",width=8, height=8)
IMM$Ifnar2_val<-as.double(IMM$Ifnar2)
ggplot(IMM,aes(x=group,y=Ifnar2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Myd88.pdf",width=8, height=8)
IMM$Myd88_val<-as.double(IMM$Myd88)
ggplot(IMM,aes(x=group,y=Myd88_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tnf.pdf",width=8, height=8)
IMM$Tnf_val<-as.double(IMM$Tnf)
ggplot(IMM,aes(x=group,y=Tnf_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Il1b.pdf",width=8, height=8)
IMM$Il1b_val<-as.double(IMM$Il1b)
ggplot(IMM,aes(x=group,y=Il1b_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Il6.pdf",width=8, height=8)
IMM$Il6_val<-as.double(IMM$Il6)
ggplot(IMM,aes(x=group,y=Il6_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Ifna4.pdf",width=8, height=8)
IMM$Ifna4_val<-as.double(IMM$Ifna4)
ggplot(IMM,aes(x=group,y=Ifna4_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Cxcl9.pdf",width=8, height=8)
IMM$Cxcl9_val<-as.double(IMM$Cxcl9)
ggplot(IMM,aes(x=group,y=Cxcl9_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Irf3.pdf",width=8, height=8)
IMM$Irf3_val<-as.double(IMM$Irf3)
ggplot(IMM,aes(x=group,y=Irf3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tlr1.pdf",width=8, height=8)
IMM$Tlr1_val<-as.double(IMM$Tlr1)
ggplot(IMM,aes(x=group,y=Tlr1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tlr5.pdf",width=8, height=8)
IMM$Tlr5_val<-as.double(IMM$Tlr5)
ggplot(IMM,aes(x=group,y=Tlr5_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tlr4.pdf",width=8, height=8)
IMM$Tlr4_val<-as.double(IMM$Tlr4)
ggplot(IMM,aes(x=group,y=Tlr4_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tlr2.pdf",width=8, height=8)
IMM$Tlr2_val<-as.double(IMM$Tlr2)
ggplot(IMM,aes(x=group,y=Tlr2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Tlr6.pdf",width=8, height=8)
IMM$Tlr6_val<-as.double(IMM$Tlr6)
ggplot(IMM,aes(x=group,y=Tlr6_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Cd14.pdf",width=8, height=8)
IMM$Cd14_val<-as.double(IMM$Cd14)
ggplot(IMM,aes(x=group,y=Cd14_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Cd25.pdf",width=8, height=8)
IMM$Cd25_val<-as.double(IMM$Cd25)
ggplot(IMM,aes(x=group,y=Cd25_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Foxp3.pdf",width=8, height=8)
IMM$Foxp3_val<-as.double(IMM$Foxp3)
ggplot(IMM,aes(x=group,y=Foxp3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Cd45.pdf",width=8, height=8)
IMM$Cd45_val<-as.double(IMM$Cd45)
ggplot(IMM,aes(x=group,y=Cd45_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Il33.pdf",width=8, height=8)
IMM$Il33_val<-as.double(IMM$Il33)
ggplot(IMM,aes(x=group,y=Il33_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Il25.pdf",width=8, height=8)
IMM$Il25_val<-as.double(IMM$Il25)
ggplot(IMM,aes(x=group,y=Il25_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Il4.pdf",width=8, height=8)
IMM$Il4_val<-as.double(IMM$Il4)
ggplot(IMM,aes(x=group,y=Il4_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Il13.pdf",width=8, height=8)
IMM$Il13_val<-as.double(IMM$Il13)
ggplot(IMM,aes(x=group,y=Il13_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Pik3ca.pdf",width=8, height=8)
IMM$Pik3ca_val<-as.double(IMM$Pik3ca)
ggplot(IMM,aes(x=group,y=Pik3ca_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Pik3cd.pdf",width=8, height=8)
IMM$Pik3cd_val<-as.double(IMM$Pik3cd)
ggplot(IMM,aes(x=group,y=Pik3cd_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Pik3r2.pdf",width=8, height=8)
IMM$Pik3r2_val<-as.double(IMM$Pik3r2)
ggplot(IMM,aes(x=group,y=Pik3r2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Pik3r1.pdf",width=8, height=8)
IMM$Pik3r1_val<-as.double(IMM$Pik3r1)
ggplot(IMM,aes(x=group,y=Pik3r1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Irak1.pdf",width=8, height=8)
IMM$Irak1_val<-as.double(IMM$Irak1)
ggplot(IMM,aes(x=group,y=Irak1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Fadd.pdf",width=8, height=8)
IMM$Fadd_val<-as.double(IMM$Fadd)
ggplot(IMM,aes(x=group,y=Fadd_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Nfkb1.pdf",width=8, height=8)
IMM$Nfkb1_val<-as.double(IMM$Nfkb1)
ggplot(IMM,aes(x=group,y=Nfkb1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("IMM.Ikbke.pdf",width=8, height=8)
IMM$Ikbke_val<-as.double(IMM$Ikbke)
ggplot(IMM,aes(x=group,y=Ikbke_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()


setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
apone<-read.table("AP1.txt")
colnames(apone)<-"X"
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
exprSet_symbol2<-read.csv("GSE126231_uniq.csv")
expdata<-exprSet_symbol2[,2:9]
row.names(expdata)<-exprSet_symbol2[,1]
apone<-merge(apone,exprSet_symbol2,order=X)
library("ggplot2")
library("ggsignif")
## AP1 regulate
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
APO<-t(apone)
APO<-APO[2:dim(APO)[1],]
APO<-cbind(APO,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(APO)<-c(apone$X,"group")
APO<-as.data.frame(APO)
APO<-APO[1:6,]
compaired <- list(c("Control", "Expandation"))
pdf("APO.Fos.pdf",width=8, height=8)
APO$Fos_val<-as.double(APO$Fos)
ggplot(APO,aes(x=group,y=Fos_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("APO.Jund.pdf",width=8, height=8)
APO$Jund_val<-as.double(APO$Jund)
ggplot(APO,aes(x=group,y=Jund_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("APO.Fosb.pdf",width=8, height=8)
APO$Fosb_val<-as.double(APO$Fosb)
ggplot(APO,aes(x=group,y=Fosb_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("APO.Fosl1.pdf",width=8, height=8)
APO$Fosl1_val<-as.double(APO$Fosl1)
ggplot(APO,aes(x=group,y=Fosl1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()

setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
tead<-read.table("TEAD.txt")
colnames(tead)<-"X"
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
exprSet_symbol2<-read.csv("GSE126231_uniq.csv")
expdata<-exprSet_symbol2[,2:9]
row.names(expdata)<-exprSet_symbol2[,1]
tead<-merge(tead,exprSet_symbol2,order=X)
library("ggplot2")
library("ggsignif")
## TEAD regulate
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
TEA<-t(tead)
TEA<-TEA[2:dim(TEA)[1],]
TEA<-cbind(TEA,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(TEA)<-c(tead$X,"group")
TEA<-as.data.frame(TEA)
TEA<-TEA[1:6,]
compaired <- list(c("Control", "Expandation"))
pdf("TEA.Cyr61.pdf",width=8, height=8)
TEA$Cyr61_val<-as.double(TEA$Cyr61)
ggplot(TEA,aes(x=group,y=Cyr61_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("TEA.Tead1.pdf",width=8, height=8)
TEA$Tead1_val<-as.double(TEA$Tead1)
ggplot(TEA,aes(x=group,y=Tead1_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("TEA.Tead2.pdf",width=8, height=8)
TEA$Tead2_val<-as.double(TEA$Tead2)
ggplot(TEA,aes(x=group,y=Tead2_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("TEA.Tead3.pdf",width=8, height=8)
TEA$Tead3_val<-as.double(TEA$Tead3)
ggplot(TEA,aes(x=group,y=Tead3_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("TEA.Tead4.pdf",width=8, height=8)
TEA$Tead4_val<-as.double(TEA$Tead4)
ggplot(TEA,aes(x=group,y=Tead4_val,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()



# export up-regulated genes both in expandation and TPA
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
cellsurvival<-read.table("cellsurvival.txt")
colnames(cellsurvival)<-"X"
cellularstress<-read.table("cellularstress.txt")
colnames(cellularstress)<-"X"
Cytoskeleton<-read.table("Cytoskeleton.txt")
colnames(Cytoskeleton)<-"X"
ECM<-read.table("ECM.txt")
colnames(ECM)<-"X"
adhesion<-read.table("adhesion.txt")
colnames(adhesion)<-"X"
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
exprSet_symbol2<-read.csv("GSE126231_uniq.csv")
expdata<-exprSet_symbol2[,2:9]
row.names(expdata)<-exprSet_symbol2[,1]
compaired <- list("Expandation","TPA")
pdf("adhesion.Itga3.pdf",width=8, height=8)
ggplot(adhesion,aes(group,Itga3,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()











compaired <- list(c("Control", "Expandation"),c("Control","TPA"))
pdf("cellsurvival.Aurkb.pdf",width=8, height=8)
ggplot(CS,aes(group,Itga3,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()

setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
cellsurvival<-merge(cellsurvival,exprSet_symbol2,order=X)
cellularstress<-merge(cellularstress,exprSet_symbol2,order=X)
Cytoskeleton<-merge(Cytoskeleton,exprSet_symbol2,order=X)
ECM<-merge(ECM,exprSet_symbol2,order=X)
library("ggplot2")
library("ggsignif")
## cellsurvival
CS<-t(cellsurvival)
CS<-CS[2:dim(CS)[1],]
CS<-cbind(CS,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(CS)<-c(cellsurvival$X,"group")
write.table(CS,"CS_t.txt",sep="\t",quote=F,row.names=F)
CS<-read.table("CS_t.txt",sep="\t",header=T)
compaired <- list(c("Control", "Expandation"),c("Control","TPA"))
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
pdf("cellsurvival.Aurkb.pdf",width=8, height=8)
ggplot(CS,aes(group,Aurkb,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellsurvival.Ccna2.pdf",width=8, height=8)
ggplot(CS,aes(group,Ccna2,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellsurvival.Ccnb1.pdf",width=8, height=8)
ggplot(CS,aes(group,Ccnb1,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellsurvival.Ccnb2.pdf",width=8, height=8)
ggplot(CS,aes(group,Ccnb2,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellsurvival.Cdk1.pdf",width=8, height=8)
ggplot(CS,aes(group,Cdk1,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellsurvival.Cdk6.pdf",width=8, height=8)
ggplot(CS,aes(group,Cdk6,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellsurvival.Mapk6.pdf",width=8, height=8)
ggplot(CS,aes(group,Mapk6,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
## cellular stress
CS<-t(cellularstress)
CS<-CS[2:dim(CS)[1],]
CS<-cbind(CS,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(CS)<-c(cellularstress$X,"group")
write.table(CS,"CS_t.txt",sep="\t",quote=F,row.names=F)
CS<-read.table("CS_t.txt",sep="\t",header=T)
compaired <- list(c("Control", "Expandation"),c("Control","TPA"))
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
pdf("cellularstress.Krt16.pdf",width=8, height=8)
ggplot(CS,aes(group,Krt16,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellularstress.Krt6a.pdf",width=8, height=8)
ggplot(CS,aes(group,Krt6a,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellularstress.Krt6b.pdf",width=8, height=8)
ggplot(CS,aes(group,Krt6b,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellularstress.S100a8.pdf",width=8, height=8)
ggplot(CS,aes(group,S100a8,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellularstress.S100a9.pdf",width=8, height=8)
ggplot(CS,aes(group,S100a9,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellularstress.Sprr1a.pdf",width=8, height=8)
ggplot(CS,aes(group,Sprr1a,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("cellularstress.Sprr1b.pdf",width=8, height=8)
ggplot(CS,aes(group,Sprr1b,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
## Cytoskeleton
CS<-t(Cytoskeleton)
CS<-CS[2:dim(CS)[1],]
CS<-cbind(CS,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(CS)<-c(Cytoskeleton$X,"group")
write.table(CS,"CS_t.txt",sep="\t",quote=F,row.names=F)
CS<-read.table("CS_t.txt",sep="\t",header=T)
compaired <- list(c("Control", "Expandation"),c("Control","TPA"))
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
pdf("Cytoskeleton.Actr3.pdf",width=8, height=8)
ggplot(CS,aes(group,Actr3,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("Cytoskeleton.Flna.pdf",width=8, height=8)
ggplot(CS,aes(group,Flna,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
## ECM
CS<-t(ECM)
CS<-CS[2:dim(CS)[1],]
CS<-cbind(CS,c(rep("Control",3),rep("Expandation",3),rep("TPA",2)))
colnames(CS)<-c(ECM$X,"group")
write.table(CS,"CS_t.txt",sep="\t",quote=F,row.names=F)
CS<-read.table("CS_t.txt",sep="\t",header=T)
compaired <- list(c("Control", "Expandation"),c("Control","TPA"))
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
pdf("ECM.Fscn1.pdf",width=8, height=8)
ggplot(CS,aes(group,Fscn1,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("ECM.Macf1.pdf",width=8, height=8)
ggplot(CS,aes(group,Macf1,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("ECM.Mmp9.pdf",width=8, height=8)
ggplot(CS,aes(group,Mmp9,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()
pdf("ECM.Pdpn.pdf",width=8, height=8)
ggplot(CS,aes(group,Pdpn,fill=group))+geom_boxplot(width=0.5)+theme(plot.title=element_text(size = 25),axis.text.x=element_text(size=15,angle=0),axis.text.y=element_text(size=15),axis.title.x=element_text(size = 23),axis.title.y=element_text(size = 23))+labs(x='Gene', y= 'RNA expression level')+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = t.test)
dev.off()

library(Hmisc)
capitalize(y)
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/celfiles")
exprSet_symbol2<-read.csv("GSE126231_uniq.csv")
expdata<-exprSet_symbol2[,2:9]
row.names(expdata)<-exprSet_symbol2[,1]
setwd("E:/stretch-mediatedskinexpansion/AffymetrixHTMG-430/geneanno")
adhesion<-read.table("adhesion.genes.uniq.txt")
GTPase<-read.table("GTPase.genes.uniq.txt")
cytoskeleton<-read.table("cytoskeleton.genes.uniq.txt")
adhesion<-as.matrix(adhesion)
adhesion_gene<-capitalize(tolower(adhesion))
adhesion_exp<-expdata[adhesion_gene,]
adhesion_exp_rmNA<-na.omit(adhesion_exp)

GTPase<-as.matrix(GTPase)
GTPase_gene<-capitalize(tolower(GTPase))
GTPase_exp<-expdata[GTPase_gene,]
GTPase_exp_rmNA<-na.omit(GTPase_exp)


cytoskeleton<-as.matrix(cytoskeleton)
cytoskeleton_gene<-capitalize(tolower(cytoskeleton))
cytoskeleton_exp<-expdata[cytoskeleton_gene,]
cytoskeleton_exp_rmNA<-na.omit(cytoskeleton_exp)






