# rm(list = ls())
library(GenomicRanges)
library(ggplot2)
source("/home/rohit/Projects/2/Mmu/utilities.R")
source("/home/rohit/Projects/2/Mmu/find_closest_gene.R")
flpath = '/mnt/cube/scratch/brachyury/Mmu/ChIP/mm10/'

# Read Peakzilla files of two antibodies
flatb1 = read.table(paste0(flpath,"SRR1168501-antibody1.mm10.c1s0f200.tsv.bed"),  sep="\t", header=TRUE)
flatb2 = read.table(paste0(flpath,"SRR1168502-antibody2.mm10.c1s0f200.tsv.bed"),  sep="\t", header=TRUE)

common.peak.dataframe <- function(peak1, peak2, score=0, foldenrichment=0){
  # Filter datasets accordigly
  atbd1 <- peak1[(peak1$FoldEnrichment >= 2) & (peak1$Score >= 0),]
  atbd2 <- peak2[(peak2$FoldEnrichment >= 2) & (peak2$Score >= 0),]
  # Convert data frames to Genomic Ranges Objects
  peak1GR <- makeGRangesFromDataFrame(atbd1, keep.extra.columns=TRUE, 
                                      ignore.strand=TRUE, seqinfo=NULL, 
                                      seqnames.field=c("Chromosome"),
                                      start.field=c("Start"), 
                                      end.field=c("End"),
                                      starts.in.df.are.0based=FALSE)
  peak2GR <- makeGRangesFromDataFrame(atbd2, keep.extra.columns=TRUE, 
                                      ignore.strand=TRUE, seqinfo=NULL, 
                                      seqnames.field=c("Chromosome"),
                                      start.field=c("Start"), 
                                      end.field=c("End"),
                                      starts.in.df.are.0based=FALSE)

  # Find overlapping genomic ranges
  common.peaks <- subsetByOverlaps(peak1GR, peak2GR, minoverlap=1)
  return(common.peaks)
}

common.sig.peaks <- common.peak.dataframe(flatb1,flatb2, score=0, foldenrichment = 2)
# Wite to a file
write.table(common.sig.peaks,
            file=paste0(flpath,"Mmu_Bra_peaks_common_peaks.txt"),
            sep="\t",append=FALSE,quote=FALSE,row.names=FALSE)

############################################## END ################################################