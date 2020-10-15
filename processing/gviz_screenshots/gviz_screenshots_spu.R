# rm(list = ls())
# 
# This scipt uses large memory therefore it is advised that it should prferably be run on cluster 
# 
# Load packages -----------------------------------------------------------

pacman::p_load(Gviz, grid, mgcv,Hmisc, data.table,GenomicFeatures)

# Load Files --------------------------------------------------------------
# Set working directory where the bedGraph files are (optional) 
setwd("/home/dnyansagar/Bra/processing/txdb")
gff <- '/export/home/dnyansagar/Bra/processing/latestSpu/GCF_000002235.5_Spur_5.0_genomic.gff'
peak_fl <- '/export/home/dnyansagar/Bra/processing/00_data/spur5_bra_chip_24hpf.targets2.txt'

# Make TxDb object from gff3 file  

txdbs <- list(spur=makeTxDbFromGFF(gff))
peak.files <- list(spur=read.table(peak_fl,
                  header=T,sep='\t',fill=T))



# read all the bedgraph files in and convert them to data frames ----------
#
# Generate bedGraph files as follows
# for f in data/*.srt.bam; do echo "bedtools genomecov -ibam $f -bg -g genome.fasta |gzip - > $(basename $f .srt.bam).bedGraph.gz" ;done 
#
read.bg <- function(short.id) {
  df <- read.table(gzfile(paste0(short.id,'.bedGraph.gz')), 
                   col.names = c('chromosome', 'start', 'end', 'value'),
                   colClasses=c('character','numeric','numeric','numeric'))
  # convert to a data table afterward. so far fread doesn't like gzipped files
  setDT(df) 
  df
}
read.all <- function(species) {
  lapply(c('CHIP_DNA_BRA_24ha_spur5','CHIP_DNA_BRA_24hM_spur5', 
           'INPUT_DNA_24ha_spur5', 'INPUT_DNA_24hM_spur5'),function(x) read.bg(x))
}
bedGraphData <- sapply(c('spur'),read.all,simplify = F,USE.NAMES = T)

# Functions ---------------------------------------------------------------

get.data.track <- function(bg,species,replicate,chrom, shade, tkname) {
  print(tkname)
  DataTrack(
    range = bg[[replicate]][bg[[replicate]]$chromosome == chrom,],
    name=tkname, type = "histogram", fill.histogram = shade,
    col.histogram = shade, window=-2, windowSize=2, 
    cex.axis=0.8,
    rotation.title=1,
    cex.title = 1.2,
    ylim =c(0,10),   genome = 'species')
}


plot.species <- function(species,txid,flank=1000,add=T, gene) {
  print(gene)
  txdb <- txdbs[[species]]
  peak.file <- peak.files[[species]]
  bg <- bedGraphData[[species]]
  txs <- transcripts(txdb,filter=list(tx_name=txid))
  peak <- peak.file[peak.file$GeneID == txid | peak.file$GeneID.1 == txid,]
  start = min(min(peak$start),start(txs))-flank
  end = max(max(peak$end),end(txs))+flank
  chrom = peak$seqnames[1]
  options(ucscChromosomeNames = F)
  
  dtrack1 <- get.data.track(bg,species,1,chrom, 'black', 'AB1')
  dtrack2 <- get.data.track(bg,species,2,chrom, 'black', 'AB2')
  dtrack3 <- get.data.track(bg,species,3,chrom, 'grey50', 'Input1')
  dtrack4 <- get.data.track(bg,species,4,chrom, 'grey50', 'Input2')
  
  if (missing(gene)){ name=''} else {name=paste0(' (', gene, ')')}
  
  ptrack <- AnnotationTrack(start=peak$start,end=peak$end, chromosome = chrom,
                            lineheight=0, min.height=1, fill="red", 
                            rotation.title=1,
                            name = "Peaks",
                            cex.title=1)
  gtrack <- GenomeAxisTrack(col="black")
  
  grtrack <- GeneRegionTrack( txdb, chromosome = chrom, #start = start, end = end,
                              fill = "navy", col="navy",
                              showId = TRUE,name = "Gene",
                              cex.title=1.25,
                              cex.group=1.25,
                              rotation.title=1,
                              fontface.group=1.25,
                              background.panel="white",
                              fontcolor.group="black",
                              col.line="navy"
  )
  
  plotTracks(list(gtrack, dtrack2, dtrack4,dtrack1, dtrack3, ptrack, grtrack),  
             from = start, to = end, background.title='white',
             col.axis='black', col.title='black', add=add,
             title.width = 0.9,
             transformation = function(x) {x/median(x)},
             main=paste0(capitalize(species), name)
  )
}

# Plot genes of interest -----------------------------------------

plot.species('spur', 'XM_030994545.1',add=F,flank=30000, 'brachyury')
plot.species('spur', 'XM_003724498.3',add=F,flank=10000, 'Fzd8')
plot.species('spur', 'XM_003724513.3',add=F,flank=10000, 'Fzd1')
plot.species('spur', 'XM_030981289.1',add=F,flank=10000, 'Dach1')
plot.species('spur', 'NM_001172052.1',add=F,flank=10000, 'Fgfrl1')

plot.species('spur', 'NM_001123500.1',add=F,flank=15000, 'Wnt1')
plot.species('spur', 'XM_030974737.1',add=F,flank=30000, 'Wnt3')
plot.species('spur', 'XM_030986594.1',add=F,flank=10000, 'Wnt4')
plot.species('spur', 'XM_784984.5',add=F,flank=30000, 'Wnt6')
plot.species('spur', 'XM_781958.5',add=F,flank=10000, 'Wnt7b')
plot.species('spur', 'NM_214667.1',add=F,flank=10000, 'Wnt8')
plot.species('spur', 'XM_030975014.1',add=F,flank=10000, 'Wnt9a')
plot.species('spur', 'XM_791523.4',add=F,flank=10000, 'Wnt16')
plot.species('spur', 'XM_030998663.1',add=F,flank=10000, 'Meis3')
plot.species('spur', 'XM_011668971.2',add=F,flank=10000, 'Gata4')
plot.species('spur', 'NM_214613.1',add=F,flank=10000, 'Msx2')
plot.species('spur', 'NM_001124764.1',add=F,flank=10000, 'Fgf16')
plot.species('spur', 'XM_030985467.1',add=F,flank=10000, 'Fgf17')
plot.species('spur', 'XM_011684108.2',add=F,flank=10000, 'Prickle1')
plot.species('spur', 'XM_011684110.2',add=F,flank=10000, 'Prickle2')
plot.species('spur', 'XM_030981057.1',add=F,flank=10000, 'Elf2')
plot.species('spur', 'XM_030981538.1',add=F,flank=10000, 'Dkc1')
plot.species('spur', 'XM_030981782.1',add=F,flank=10000, 'Zic3')
plot.species('spur', 'XM_030981866.1',add=F,flank=10000, 'Fndc3b')


# END ---------------------------------------------------------------------

