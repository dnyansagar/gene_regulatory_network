#rm(list = ls())
library(GenomicRanges)
#find closest gene source files 
source("/home/rohit/Projects/2/Mmu/utilities.R")
source("/home/rohit/Projects/2/Mmu/find_closest_gene.R")
# Read Peakzilla files of two antibodies
flatb1 = read.table("/mnt/cube/scratch/brachyury/Mmu/ChIP/mm10/SRR1168501-antibody1.mm10.c1s0f200.tsv.bed",  sep="\t", header=TRUE)
flatb2 = read.table("/mnt/cube/scratch/brachyury/Mmu/ChIP/mm10/SRR1168502-antibody2.mm10.c1s0f200.tsv.bed",  sep="\t", header=TRUE)
mmu_genes <- read.table(file="/home/rohit/Projects/2/Mmu/mm10_ensembl_compr.protein_coding.supergene.srt.bed", 
                        comment.char="#", sep="\t", header=FALSE)

nemve <- read.table("/home/rohit/Projects/2/Nve/nematostella_target_genes_with.peaks_score_cutoff_7_170523.txt", sep = '\t', header = TRUE)
nve_genes <- read.table(file="/home/rohit/Projects/2/Nve/nveGenes.good.130208.longCDS.bed", comment.char="#", sep="\t", header=FALSE)
  
xentr <- read.table("/home/rohit/Projects/2/Xtr/xenopus_Bra_peaks_score_cutoff_2_160108.txt", sep = '\t', header = TRUE)
xtr_genes <- read.table(file="/home/rohit/Projects/2/Xtr/XenTro4.1_ensembl_compr.protein_coding.supergene.bed", comment.char="#", sep="\t", header=FALSE)


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
peak.dataframe        <- function(selected_peaks){
  # make GRanges object from peaks
  selected_peaks <- makeGRangesFromDataFrame(selected_peaks, keep.extra.columns=TRUE, 
                                             ignore.strand=TRUE, seqinfo=NULL, 
                                             seqnames.field=c("seqname"),
                                             start.field=c("start"), 
                                             end.field=c("end"),
                                             starts.in.df.are.0based=FALSE)
  return(selected_peaks)
}
get.nearest.tss       <- function(peak_fl, genes_fl){
  # Make data frames from gene bed file to Genomic Ranges Objects
  genes.features <- makeGRangesFromDataFrame(genes_fl,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             seqinfo=NULL,
                                             seqnames.field=c("V1"),
                                             start.field=c("V2"),
                                             end.field=c("V3"),
                                             strand.field="V6",
                                             starts.in.df.are.0based=FALSE)
  # Get transcription start site
  tra.start.site <- flank(genes.features,1)
  # Peak closest to transcription start site
  nearest.target.with.peaks <- getNearestFeatureIndicesAndDistances(peak_fl, tra.start.site)
  # nearest.target.peaks contain index, distance, and peak info
  targets <- tra.start.site[nearest.target.with.peaks$index]
  # Targets contain all columns of genes bed file
  # Make data frame of 'Target', 'seqnames', 'start', 'end', 'strand'
  target_df <- data.frame(GeneID = elementMetadata(targets)[,"V4"],
                          seqnames=seqnames(targets),
                          start=start(targets),end=end(targets),
                          strand=strand(targets))
  peak_df <- data.frame(seqnames=seqnames(nearest.target.with.peaks$peak),
                     start=start(nearest.target.with.peaks$peak),
                     end=end(nearest.target.with.peaks$peak), 
                     Summit=elementMetadata(nearest.target.with.peaks$peak)[,"Summit"],
                     Score=elementMetadata(nearest.target.with.peaks$peak)[,"Score"], 
                     Enrichment =elementMetadata(nearest.target.with.peaks$peak)[,"FoldEnrichment"],
                     distance=nearest.target.with.peaks$distance)
  target_peaks <- cbind(target_df,peak_df)
  return(target_peaks)
}
process               <- function(bed.dataframe, peak.dataframe){
  # Get first target
  first <- get.nearest.tss(peak.dataframe, bed.dataframe)
  # Remove first targets from the bed file
  bed.dataframe.fr <- bed.dataframe[!(bed.dataframe$V4 %in% first$GeneID),]
  # Get second target
  second <- get.nearest.tss(peak.dataframe, bed.dataframe.fr)
  # Merge both targets info 
  both.targets <- merge(first, second, by= "Summit", all.x = TRUE)
  # Get selected columns
  targets.selected.columns <- both.targets[
    c("seqnames.x", "start.x", "end.x", "Summit", "Score.x", "Enrichment.x","GeneID.x",
      "seqnames.x", "start.x", "end.x", "strand.x", "distance.x", "GeneID.y",	
      "seqnames.y", "start.y",	"end.y", "strand.y", "distance.y")]
  return(targets.selected.columns)
}


common.sig.peaks <- common.peak.dataframe(flatb1,flatb2, score=0, foldenrichment = 2)
mouse_both <- process(genes, common.sig.peaks)

xtr_peaks <- peak.dataframe(xentr)
write.table(selected.columns,
            file="/home/rohit/Projects/2/Mmu/Mmu_Bra_peaks_common_peaks_targets_trial.txt",
            sep="\t",append=FALSE,quote=FALSE,row.names=FALSE)

############################################## END ################################################