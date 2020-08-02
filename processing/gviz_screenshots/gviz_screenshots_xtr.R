# rm(list = ls())

# Load packages -----------------------------------------------------------

pacman::p_load(Gviz, grid, mgcv,Hmisc, data.table,GenomicFeatures)


# Load Files --------------------------------------------------------------

setwd('/home/dnyansagar/Bra/processing/txdb')
txdbs <- list(
  Xtro=makeTxDbFromGFF('GCF_000004195.3_Xenopus_tropicalis_v9.1_genomic.gff'))
peak.files <- list(
  Xtro=read.table('Xenopus_Bra_peaks_score_cutoff_181121.targets1.fil.bed',
                 header=T,sep='\t',fill=T))


# read all the bedgraph files in and convert them to data frames ----------

read.bg <- function(short.id) {
  df <- read.table(gzfile(paste0(short.id,'.aln.bedgraph.gz')), 
                   col.names = c('chromosome', 'start', 'end', 'value'),
                   colClasses=c('character','numeric','numeric','numeric'))
  # convert to a data table afterward. so far fread doesn't like gzipped files
  setDT(df) 
  df
}
read.all <- function(species) {
  lapply(c('I1','C1','I2', 'C2'),function(x) read.bg(paste0(species,x)))
}
bedGraphData <- sapply(c('Xtro'),read.all,simplify = F,USE.NAMES = T)


# Functions ---------------------------------------------------------------

get.data.track <- function(bg,species,replicate,chrom, shade, tkname) {
  print(tkname)
  DataTrack(
    range = bg[[replicate]][bg[[replicate]]$chromosome == chrom,],
    name=tkname, type = "histogram", fill.histogram = shade,
    col.histogram = shade, window=-1, windowSize=1, 
    cex.axis=.8,
    cex.title=3, 
    cex=0.7,
    cex.legend=0.9,
    fontsize.legend=0.9,
    ylim =c(0,100),   genome = 'species')
}
# 
plot.species <- function(species, txid, protid, flank=1000, add=T, gene) {
  txdb <- txdbs[[species]]
  peak.file <- peak.files[[species]]
  bg <- bedGraphData[[species]]
  txs <- transcripts(txdb,filter=list(tx_name=txid))
  peak <- peak.file[peak.file$GeneID == protid | peak.file$GeneID.1 == protid,]
  start = min(min(peak$start), start(txs))-flank
  end =   max(max(peak$end), end(txs))+flank
  chrom = peak$seqnames[1]
  options(ucscChromosomeNames = F)
  
  dtrack1 <- get.data.track(bg, species, 1, chrom, 'blue3', 'Input1') 
  dtrack2 <- get.data.track(bg, species, 2, chrom, 'black', 'AB1') 
  dtrack3 <- get.data.track(bg, species, 3, chrom, 'blue3', 'Input2') 
  dtrack4 <- get.data.track(bg, species, 4, chrom, 'black', 'AB') 
  
  if (missing(gene)){ name=''} else {name=paste0(' (', gene, ')')}
  ptrack <- AnnotationTrack(start=peak$start,end=peak$end,
                            background.panel="white",
                            #fontcolor.title="red",
                            chromosome = chrom,
                            lineheight=0,
                            min.height=1, 
                            fill="red", cex.title=3 , name = "Peaks")
  
  gtrack <- GenomeAxisTrack(col="black") 
  
  grtrack <- GeneRegionTrack( txdb, chromosome = chrom, start = start, end = end,
                              fill = "navy", col="navy", 
                              showId = TRUE, name = "Gene", 
                              cex.title=3,
                              cex.group=2,
                              fontface.group=2,
                              background.panel="white",
                              fontcolor.group="black",
                              col.line="navy",
  )
  
  plotTracks(list(gtrack,  dtrack4,  ptrack, grtrack), #  dtrack2, 
             from = start, to = end, background.title='white',
             col.axis='black', col.title='black', add=add, #slategrey
             main=paste0(capitalize(species), name),
             #littleTicks = TRUE,
             sizes=c(1.5,3,2,2) #3,
             #sizes=c(1.5,5,1,1,2)) #,5
  )
}

# Plot genes of interest --------------------------------------------------

ids1 <- c('XM_012962922.1','XP_017946799.1','NP_001107970.1','NP_001008449.1',
          'NP_001017034.2','NP_001027481.1','NP_001165369.1','NP_001016741.1',
          'NP_989197.1','NP_001016283.1','XP_002944388.2','NP_001027524.1','XP_012825712.1',
          'XP_002939354.2','XP_012814299.1','NM_001017208.2','NP_001008133.1','NP_001016735.1')
ids2 <- c('XP_012818376.1','XP_017946799.1','NP_001107970.1','NP_001008449.1',
          'XP_012823915.1','NP_001027481.1','NP_001165369.1','NP_001016741.1',
          'NP_989197.1', 'NP_001016283.1','XP_002944388.2','NP_001027524.1','XP_012825712.1',
          'XP_002939354.2','XP_012814299.1', 'NP_001017208.1','NP_001008133.1','NP_001016735.1')
name <- c( 'Bra','Tbx2','Foxb1','fzd2','Bmp4','Vangl2','Noggin','Fzd7',
           'Bmp7','Dkk1','Six3', 'Tbx3','Sp5', 'Wnt3a','Wnt5b','Wnt-8a',
           'Wnt-11b-1','Wnt-11b-2')

path = '../processing/#figures/screenshots/Screenshots_NoInput/'
for (i in 1:length(ids1)){
  flname <- paste0(path,'xtr', name[i],'.fx.jpg')
  print(flname)
  jpeg(flname,  width = 1080, height = 600)
  plot.species('Xtro', ids1[i], ids2[i], add=F, flank=50000, name[i])
  dev.off()
}


# Plot genes of interest manually -----------------------------------------
 
plot.species('Xtro', 'XM_012962922.1', 'XP_012818376.1',add=F,flank=15000, 'Bra') # 300
plot.species('Xtro', 'XP_017946799.1', 'XP_017946799.1',add=F,flank=15000, 'Tbx2')
plot.species('Xtro', 'NP_001107970.1','NP_001107970.1',add=F,flank=15000, 'Foxb1')
plot.species('Xtro', 'NP_001008449.1','NP_001008449.1',add=F,flank=15000, 'fzd2')
plot.species('Xtro', 'NP_001017034.2','XP_012823915.1',add=F,flank=15000, 'Bmp4')
plot.species('Xtro', 'NP_001027481.1','NP_001027481.1',add=F,flank=50000, 'Vangl2')
plot.species('Xtro', 'NP_001165369.1', 'NP_001165369.1',add=F,flank=15000, 'Noggin')
plot.species('Xtro', 'NP_001016741.1', 'NP_001016741.1',add=F,flank=15000, 'Fzd7')
plot.species('Xtro', 'NP_989197.1', 'NP_989197.1', add=F, flank=35000, 'Bmp7')
plot.species('Xtro', 'NP_001016283.1', 'NP_001016283.1', add=F, flank=15000, 'Dkk1')
plot.species('Xtro', 'XP_002944388.2', 'XP_002944388.2', add=F, flank=15000, 'Six3')
plot.species('Xtro', 'NP_001027524.1', 'NP_001027524.1', add=F, flank=15000, 'Tbx3')
plot.species('Xtro', 'XP_012825712.1','XP_012825712.1', add=F, flank=15000, 'Sp5')
plot.species('Xtro', 'XP_002934799.2','XP_002934799.2', add=F, flank=10000, 'Tll1')
plot.species('Xtro', 'NP_001006776.2', 'NP_001006776.2',add=F,flank=25000, 'Pax3')
plot.species('Xtro', 'XP_012823259.1', 'XP_012823259.1',add=F,flank=25000, 'Meis2')
plot.species('Xtro', 'XP_002939354.2', 'XP_002939354.2',add=F,flank=25000, 'Wnt3a')
plot.species('Xtro', 'XP_012814299.1', 'XP_012814299.1',add=F,flank=25000, 'Wnt5b')
plot.species('Xtro', 'NM_001017208.2', 'NP_001017208.1',add=F,flank=25000, 'Wnt-8a')
plot.species('Xtro', 'NP_001008133.1', 'NP_001008133.1',add=F,flank=25000, 'Wnt-11b-1')
plot.species('Xtro', 'NP_001016735.1', 'NP_001016735.1',add=F,flank=25000, 'Wnt-11b-2')

