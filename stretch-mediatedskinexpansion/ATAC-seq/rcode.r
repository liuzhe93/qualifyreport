library("knitr")
library("rmdformats")
library("dplyr")
library("DT")
library("tidyr")
library("ggplot2")
library("magrittr")
library("Rsamtools")
library("GenomicAlignments")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("soGGi")
library("rtracklayer")
library("ChIPQC")
library("ChIPseeker")
library("rGREAT")
library("limma")
library("DESeq2")
library("tracktables")
library("clusterProfiler")
library("org.Mm.eg.db")
library("MotifDb")
library("Biostrings")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Mmusculus.UCSC.mm9")
# # Finally we need development version of soGGi (named here 1.10.4) # not
# version on Bioconductor (1.10.0)
# devtools::install_github('ThomasCarroll/soGGi')
library(Rsubread)
sortedBAM_ctl<-"Ctl_S18.last.bam"
sortedBAM_exp<-"Exp-D2.last.bam"
pmapped_ctl <- propmapped(sortedBAM_ctl)
pmapped_exp <- propmapped(sortedBAM_exp)
save(pmapped_ctl,file="Rdata/pmapped_ctl.RData")
save(pmapped_exp,file="Rdata/pmapped_exp.RData")
library(Rsamtools)
library(ggplot2)
library(magrittr)
pdf("figures/Distributionofmappedreads_ctl.pdf")
idxstatsBam(sortedBAM_ctl) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()
dev.off()
pdf("figures/Distributionofmappedreads_exp.pdf")
idxstatsBam(sortedBAM_exp) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()
dev.off()
library(GenomicAlignments)
atacReads_ctl <- readGAlignmentPairs(sortedBAM_ctl, param = ScanBamParam(mapqFilter = 1, 
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
        "mapq", "isize")))
length(atacReads_ctl)
atacReads_ctl
atacReads_exp <- readGAlignmentPairs(sortedBAM_exp, param = ScanBamParam(mapqFilter = 1, 
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
        "mapq", "isize")))
length(atacReads_exp)
atacReads_exp
atacReads_ctl_read1 <- GenomicAlignments::first(atacReads_ctl)
insertSizes <- abs(elementMetadata(atacReads_ctl_read1)$isize)
head(insertSizes)
library(magrittr)
library(dplyr)
library(ggplot2)
pdf("figures/Plottingtinsertsizes_ctl.pdf")
fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
    Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
    Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
    geom_line()
fragLenPlot + theme_bw()
dev.off()
pdf("figures/PlottingtinsertsizesYlog2_ctl.pdf")
fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
dev.off()
pdf("figures/Greenleaf_ctl.pdf")
fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 
    437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
    theme_bw()
dev.off()
pdf("figures/GreenleafYlog2_ctl.pdf")
fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 
    247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
    geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
dev.off()
atacReads_exp_read1 <- GenomicAlignments::first(atacReads_exp)
insertSizes <- abs(elementMetadata(atacReads_exp_read1)$isize)
head(insertSizes)
pdf("figures/Plottingtinsertsizes_exp.pdf")
fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
    Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
    Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
    geom_line()
fragLenPlot + theme_bw()
dev.off()
pdf("figures/PlottingtinsertsizesYlog2_exp.pdf")
fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
dev.off()
pdf("figures/Greenleaf_exp.pdf")
fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 
    437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
    theme_bw()
dev.off()
pdf("figures/GreenleafYlog2_exp.pdf")
fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 
    247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
    geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
dev.off()



library("TxDb.Mmusculus.UCSC.mm9.knownGene")
TSSs <- resize(genes(TxDb.Mmusculus.UCSC.mm9.knownGene), fix = "start", 1)
TSSs
library(soGGi)
#Plotting ATAC-seq signal of TSSs (Creating open, mono- and di-nucleosome signal profiles)
# Nucleosome free
nucFree_ctl <- regionPlot(bamFile = sortedBAM_ctl, testRanges = TSSs, style = "point", 
    format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100, 
    forceFragment = 50)
# Mononucleosome
monoNuc_ctl <- regionPlot(bamFile = sortedBAM_ctl, testRanges = TSSs, style = "point", 
    format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240, 
    forceFragment = 80)
# Dinucleosome
diNuc_ctl <- regionPlot(bamFile = sortedBAM_ctl, testRanges = TSSs, style = "point", 
    format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437, 
    forceFragment = 160)
save(monoNuc_ctl,nucFree_ctl,diNuc_ctl,file='Rdata/soGGiResults_ctl.RData')

#Plotting ATAC-seq signal of TSSs (Creating open, mono- and di-nucleosome signal profiles)
# Nucleosome free
nucFree_exp <- regionPlot(bamFile = sortedBAM_exp, testRanges = TSSs, style = "point", 
    format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100, 
    forceFragment = 50)
# Mononucleosome
monoNuc_exp <- regionPlot(bamFile = sortedBAM_exp, testRanges = TSSs, style = "point", 
    format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240, 
    forceFragment = 80)
# Dinucleosome
diNuc_exp <- regionPlot(bamFile = sortedBAM_exp, testRanges = TSSs, style = "point", 
    format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437, 
    forceFragment = 160)
save(monoNuc_exp,nucFree_exp,diNuc_exp,file='Rdata/soGGiResults_exp.RData')

library(soGGi)
load(file = "Rdata/soGGiResults_ctl.RData")
#Plotting ATAC-seq signal of TSSs (Plotting open, mono- and di-nucleosome signal profiles)
pdf("figures/PlottingATAC-seqsignalofTSSs_nucFree_ctl.pdf")
plotRegion(nucFree_ctl, outliers = 0.01)
dev.off()
pdf("figures/PlottingATAC-seqsignalofTSSs_monoNuc_ctl.pdf")
plotRegion(monoNuc_ctl, outliers = 0.01)
dev.off()
pdf("figures/PlottingATAC-seqsignalofTSSs_diNuc_ctl.pdf")
plotRegion(diNuc_ctl, outliers = 0.01)
dev.off()

load(file = "Rdata/soGGiResults_exp.RData")
#Plotting ATAC-seq signal of TSSs (Plotting open, mono- and di-nucleosome signal profiles)
pdf("figures/PlottingATAC-seqsignalofTSSs_nucFree_exp.pdf")
plotRegion(nucFree_exp, outliers = 0.01)
dev.off()
pdf("figures/PlottingATAC-seqsignalofTSSs_monoNuc_exp.pdf")
plotRegion(monoNuc_exp, outliers = 0.01)
dev.off()
pdf("figures/PlottingATAC-seqsignalofTSSs_diNuc_exp.pdf")
plotRegion(diNuc_exp, outliers = 0.01)
dev.off()
insertSizes <- abs(elementMetadata(atacReads_ctl_read1)$isize)
atacReads_Open_ctl <- atacReads_ctl[insertSizes < 100, ]
atacReads_MonoNuc_ctl <- atacReads_ctl[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc_ctl <- atacReads_ctl[insertSizes > 315 & insertSizes < 437, ]
openRegionBam_ctl <- gsub("\\.bam", "_openRegions\\.bam", sortedBAM_ctl)
monoNucBam_ctl <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM_ctl)
diNucBam_ctl <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM_ctl)
export(atacReads_Open_ctl, openRegionBam_ctl, format = "bam")
export(atacReads_MonoNuc_ctl, monoNucBam_ctl, format = "bam")
# export(atacReads_Open_ctl,diNucBam_ctl,format = 'bam')
openRegionBigWig_ctl <- gsub("\\.bam", "_openRegions\\.bw", sortedBAM_ctl)
openRegionRPMBigWig_ctl <- gsub("\\.bam", "_openRegionsRPM\\.bw", sortedBAM_ctl)
atacFragments_Open_ctl <- granges(atacReads_Open_ctl)
export.bw(coverage(atacFragments_Open_ctl), openRegionBigWig_ctl)
insertSizes <- abs(elementMetadata(atacReads_exp_read1)$isize)
atacReads_Open_exp <- atacReads_exp[insertSizes < 100, ]
atacReads_MonoNuc_exp <- atacReads_exp[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc_exp <- atacReads_exp[insertSizes > 315 & insertSizes < 437, ]
openRegionBam_exp <- gsub("\\.bam", "_openRegions\\.bam", sortedBAM_exp)
monoNucBam_exp <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM_exp)
diNucBam_exp <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM_exp)
export(atacReads_Open_exp, openRegionBam_exp, format = "bam")
export(atacReads_MonoNuc_exp, monoNucBam_exp, format = "bam")
# export(atacReads_Open_exp,diNucBam_exp,format = 'bam')
openRegionBigWig_exp <- gsub("\\.bam", "_openRegions\\.bw", sortedBAM_exp)
openRegionRPMBigWig_exp <- gsub("\\.bam", "_openRegionsRPM\\.bw", sortedBAM_exp)
atacFragments_Open_exp <- granges(atacReads_Open_exp)
export.bw(coverage(atacFragments_Open_exp), openRegionBigWig_exp)

library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)

blkList <- import.bed("ATAC_blacklists/mm9-blacklist.bed")
openRegionPeaks <- "peakcalling/ATAC_openRegions_peaks.narrowPeak"
qcRes <- ChIPQCsample("Exp-D2.last_openRegions.bam", 
    peaks = openRegionPeaks, annotation = "mm9",
    verboseT = FALSE)
QCmetrics(qcRes) %>% t %>% data.frame %>% dplyr:::select(Reads, starts_with(c("Filt")), 
    starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)
flagtagcounts(qcRes) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
    100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)
MacsCalls_chr20 <- granges(qcRes[seqnames(qcRes) %in% "chr20"])

data.frame(Blacklisted = sum(MacsCalls_chr20 %over% blkList), Not_Blacklisted = sum(!MacsCalls_chr20 %over% 
    blkList))
MacsCalls_chr20_filtered <- MacsCalls_chr20[!MacsCalls_chr20 %over% blkList]

library(ChIPseeker)
MacsCalls_chr20_filteredAnno <- annotatePeak(MacsCalls_chr20_filtered, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
MacsCalls_chr20_filteredAnno

plotAnnoPie(MacsCalls_chr20_filteredAnno)

plotAnnoBar(MacsCalls_chr20_filteredAnno)

upsetplot(MacsCalls_chr20_filteredAnno)

MacsGranges_Anno <- as.GRanges(MacsCalls_chr20_filteredAnno)
TSS_MacsGranges_Anno <- MacsGranges_Anno[abs(MacsGranges_Anno$distanceToTSS) < 
    500]
TSS_MacsGranges_Anno
library(rGREAT)
seqlevelsStyle(MacsCalls_chr20_filtered) <- "UCSC"

great_Job <- submitGreatJob(MacsCalls_chr20_filtered, species = "hg19")
availableCategories(great_Job)

library(ATACseqQC)
bamfile <- "Ctl_S18.last.bam"
bamfile.labels <- gsub(".bam", "", basename(bamfile))
## generate fragement size distribution
fragSize <- fragSizeDist(bamfile, bamfile.labels)
## shift the coordinates of 5'ends of alignments in the bam file
library(BSgenome.Hsapiens.UCSC.hg19)
Ctl_S18 <- readBamFile(bamfile, asMates=TRUE)
Ctl_S181 <- shiftGAlignmentsList(Ctl_S18)
shiftedBamfile <- file.path(outPath, "shifted.bam")
export(Ctl_S181, shiftedBamfile)

