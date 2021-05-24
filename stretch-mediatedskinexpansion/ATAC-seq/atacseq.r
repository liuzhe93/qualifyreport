setwd("E:/stretch-mediatedskinexpansion/ATAC-seq")
library(ATACseqQC)
bamfile <- "Ctl_S18.last.bam"
bamfile.labels <- "Ctl_S18"
## generate fragement size distribution
fragSize <- fragSizeDist(bamfile, bamfile.labels)

library(BSgenome.Hsapiens.UCSC.hg19)
gal <- readBamFile(bamfile, asMates=TRUE)
gal1 <- shiftGAlignmentsList(gal)
shiftedBamfile <- file.path(outPath, "shifted.bam")
export(gal1, shiftedBamfile)

