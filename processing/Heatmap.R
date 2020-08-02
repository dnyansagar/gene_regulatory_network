
# Installation of packages ------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Rsubread")
# BiocManager::install("csaw")
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# devtools::install_github("jokergoo/ComplexHeatmap")
# install.packages("~/Dropbox/Dnyansagar_Data/Manuscript_Brachyury/brachyury_data/MiniChip_0.0.0.9000.tar.gz",repos=NULL)

help(package = "MiniChip", help_type = "html")

# Load the packages -------------------------------------------------------


pacman::p_load(MiniChip, GenomicRanges, ComplexHeatmap, Cairo, tidyverse)


# Make genome ranges of the peaks -----------------------------------------


peaks_fl1 <- read.table(file="~/evo/Bra/processing/00_data/nematostella_Bra_peaks_score_cutoff_7_160108.bed", 
                        sep="\t", header=TRUE)
peaks <- makeGRangesFromDataFrame(peaks_fl1, 
              keep.extra.columns=TRUE, 
              ignore.strand=TRUE, 
              seqinfo=NULL, 
              seqnames.field=c("seqnames"),
              start.field=c("start"), 
              end.field=c("end"), 
              starts.in.df.are.0based=FALSE)


# Schweiger et. al. enhancer data -----------------------------------------


ovr_enh <- "~/evo/Bra/processing/enhancer_overlap/overlap_michi_enhancers"
overlap_enhancer <- read.table(ovr_enh, sep = "\t", header = TRUE)

ovr_pro <- "~/evo/Bra/processing/enhancer_overlap/overlap_michi_promoters"
overlap_promoter <- read.table(ovr_pro, sep = "\t", header = TRUE)

michi <- read.table('~/evo/Bra/processing/enhancer_overlap/Michi_Supplemental_Table_3.txt',
           header = T, 
           sep ='\t')

michi_gastrula_planula <- michi %>% filter(stage =='gastrula' | stage =='gastrula and planula')
michi_gastrula <- michi %>% filter(stage =='gastrula')  



# List BAM files ----------------------------------------------------------


all.bamFiles <- list.files("/Users/dnyansagar/evo/Bra/BAM/histone/", 
                           full.names=TRUE,
                           pattern="*bam$")

bamNames <- sub('\\..*$', '', basename(all.bamFiles))
span <- 3025
step <- 50
summits <- peaks
start(summits) <- start(summits) 
end(summits) <- end(summits)
names(summits) <- elementMetadata(summits)[,"Summit"]


# Generate Count around the sumits ----------------------------------------


counts <- SummitHeatmap(summits, 
                        all.bamFiles, 
                        bamNames, 
                        span, step, useCPM=TRUE, 
                        readShiftSize=100, minMQS = 0)

nicebamNames <- bamNames

# Generate HeatMaps -------------------------------------------------------


ht_list <- DrawSummitHeatmaps(counts[c(1,2,3,4)], 
                              bamNames[c(1,2,3,4)], 
                              nicebamNames[c(1,2,3,4)], 
                              orderSample=1, 
                              use.log=TRUE,
                              splitHM = overlap$Summit,
                              bottomCpm=c(0,0,0,0), 
                              medianCpm = c(1,1,1,1), 
                              topCpm = c(4,4,4,4), 
                              plotcols = c("seagreen","navy","firebrick3","indianred3"),
                              orderWindows=2)
draw(ht_list, padding = unit(c(3, 8, 8, 2), "mm"))


# End ---------------------------------------------------------------------


