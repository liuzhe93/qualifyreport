##1.ReadData##
#https://www.jianshu.com/p/87fa537e479a
#add bioconductor source
#source("https://bioconductor.org/biocLite.R")
#install affy package from bioconductor
#biocLite("affy")
#load affy package
library("affy")
#set the directory that you are working with,this can be replaced by your own CEL files path
setwd("E:/stretch-mediatedskinexpansion/Affymetrix")
#set the .CEL files path,this can be replaced by your own CEL files path
celpath="E:/stretch-mediatedskinexpansion/Affymetrix/"
#read .CEL files from directory of celfile.path
data = ReadAffy(celfile.path =celpath)
#replace the sample names of the data by trim off the ".CEL" suffix
sampleNames(data) = sub("\\.cel$", "", sampleNames(data))
#read the phenoData of the CEL file.The "type.csv" file is consitituded with .csv format after you comprehensed the samples grouping on the GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11787
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126231
samplestype <- read.csv(paste(celpath,"type.csv",sep = ""),header = TRUE,row.names = 1,na.strings = "NA",sep = ',',stringsAsFactors = F,as.is = !default.stringsAsFactors())
#replace the rownames of the samplestype by its sampleID
rownames(samplestype) = samplestype$sampleID
#match the rownames of your phenoData to the samplenames of experiments data,which will return the matched index
mt = match(rownames(samplestype), sampleNames(data))
#give a completed column description of the phenoData
vmd = data.frame(labelDescription = c("Sample ID", "Samples type: Mechanisms of stretch-mediated skin expansion at single cell resolution"))
#make the matched sampletype rows and varMetadata to your data's phenoData
phenoData(data) = new("AnnotatedDataFrame", data = samplestype[mt, ], varMetadata = vmd)
#erase samples which has no type
data = data[,!is.na(data$type)]

##2.Quality Control##
#install the essential packages from Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("affyQCReport")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("simpleaffy")
#load the pakages of affyQCReport and simpleaffy
library("affyQCReport")
library("simpleaffy")
#Execute the Quality Control
saqc = qc(data)
#plot the quality control result
plot(saqc)

