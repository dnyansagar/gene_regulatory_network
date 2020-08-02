# rm(list = ls())

# Load packages -----------------------------------------------------------

pacman::p_load(Gviz, grid, mgcv,Hmisc, data.table,GenomicFeatures)

# Load Files --------------------------------------------------------------

setwd("/home/dnyansagar/Bra/processing/txdb")

txdbs <- list(
  spur=makeTxDbFromGFF('Transcriptome.strPur2.WHLids.nooverlap.gff3'))
peak.files <- list(
  spur=read.table('sea_urchin_Bra_Ma_366peaks_score_cutoff_0.5_Spur2_170228.targets1.txt',
                  header=T,sep='\t',fill=T))

# read all the bedgraph files in and convert them to data frames ----------

read.bg <- function(short.id) {
  df <- read.table(gzfile(paste0(short.id,'.bedGraph.gz')), 
                   col.names = c('chromosome', 'start', 'end', 'value'),
                   colClasses=c('character','numeric','numeric','numeric'))
  # convert to a data table afterward. so far fread doesn't like gzipped files
  setDT(df) 
  df
}
read.all <- function(species) {
  lapply(c('C0','C1', 'C2', 'C3'),function(x) read.bg(paste0(species,x)))
}
bedGraphData <- sapply(c('spur'),read.all,simplify = F,USE.NAMES = T)

# Functions ---------------------------------------------------------------

get.data.track <- function(bg,species,replicate,chrom, shade, tkname) {
  print(tkname)
  DataTrack(
    range = bg[[replicate]][bg[[replicate]]$chromosome == chrom,],
    name=tkname, type = "histogram", fill.histogram = shade,
    col.histogram = shade, window=-1, windowSize=1, 
    cex.axis=.8,
    cex.title=3, 
    #cex.title=1, 
    ylim =c(0,50),   genome = 'species')
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
  
  dtrack1 <- get.data.track(bg,species,1,chrom, 'black', 'Input1')
  dtrack2 <- get.data.track(bg,species,2,chrom, 'black', 'AB')
  dtrack3 <- get.data.track(bg,species,3,chrom, 'black', 'Input2')
  dtrack4 <- get.data.track(bg,species,4,chrom, 'black', 'AB')
  
  if (missing(gene)){ name=''} else {name=paste0(' (', gene, ')')}
  
  ptrack <- AnnotationTrack(start=peak$start,end=peak$end, chromosome = chrom,
                            lineheight=0, min.height=1, fill="red", 
                            name = "Peaks",
                            cex.title=2.5)
  gtrack <- GenomeAxisTrack(col="black")
  
  grtrack <- GeneRegionTrack( txdb, chromosome = chrom, #start = start, end = end,
                              fill = "navy", col="navy",
                              showId = TRUE,name = "Gene",
                              cex.title=2.5,
                              cex.group=2,
                              fontface.group=2,
                              background.panel="white",
                              fontcolor.group="black",
                              col.line="navy"
  )
  
  plotTracks(list(gtrack, dtrack2, dtrack4, ptrack, grtrack), #dtrack2,
             from = start, to = end, background.title='white',
             col.axis='black', col.title='black', add=add,
             main=paste0(capitalize(species), name),
             #littleTicks = TRUE,
             sizes=c(1.5, 3, 3, 2.5, 2.5) #,3
             #sizes=c(1.5,5,5,1,1,2))
  )
}

# Plot genes of interest --------------------------------------------------

ids <- c('WHL22.124093.0',  'WHL22.169409.0','WHL22.677570.0',
         'WHL22.600041.0', 'WHL22.636100.0', 'WHL22.242915.0',
         'WHL22.738092.0', 'WHL22.606237.0','WHL22.78013.0',  
         'WHL22.754673.3','WHL22.417635.0','WHL22.438846.0', 
         'WHL22.614286.0', 'WHL22.80978.1', 'WHL22.45238.0', 
         'WHL22.343932.0', 'WHL22.532435.2','WHL22.273041.0', 
         'WHL22.306197.0', 'WHL22.647633.0','WHL22.52081.0',
         'WHL22.735232.0','WHL22.596775.0')


gene <- c('Chrd','Pdx1', 'Neurog3', 
          'Brachyury','Prickle2', 'Prdm14', 
          'Pcdh10', 'Notch4', 'GATA4', 
          'Rgmb', 'Calm5', 'Bmp2', 
          'Tgif2', 'Notum',  'Fzd2', 
          'Fgf16', 'Otx2', 'CYP2C8', 
          'Tenm4', 'Chrna2', 'Wnt5',  
          'Wnt16', 'Wnt-9a')

path = '/home/dnyansagar/Bra/processing/#figures/screenshots/Screenshots_NoInput/22Sep/temp/'

for (i in 1:length(ids)){
  flname <- paste0(path,'spur', gene[i],'.jpg')
  jpeg(flname,  width = 1080, height = 700, pointsize = 12, quality = 100)
  plot.species('spur', ids[i], add=F, flank=20000, gene[i])
  dev.off()
}

# Plot genes of interest manually -----------------------------------------

plot.species('spur', 'WHL22.124093.0',add=F,flank=30000, 'Chrd')
plot.species('spur', 'WHL22.169409.0', add=F,flank=10000, 'Pdx1')
plot.species('spur', 'WHL22.677570.0', add=F,flank=10000, 'Neurog3')
plot.species('spur', 'WHL22.42488.0', add=F,flank=10000, 'fzd8')
plot.species('spur', 'WHL22.600041.0', add=F,flank=10000, 'Brachyury')
plot.species('spur', 'WHL22.636100.0', add=F,flank=10000, 'Prickle2')
plot.species('spur', 'WHL22.2236.1', add=F,flank=10000, 'Meis1')
plot.species('spur', 'WHL22.143854.0', add=F,flank=10000, 'Isl1')
plot.species('spur', 'WHL22.40221.0', add=F,flank=10000, 'POU3f4')
plot.species('spur', 'WHL22.104525.0', add=F,flank=10000, 'Sox2')
plot.species('spur', 'WHL22.242915.0', add=F,flank=10000, 'Prdm14')
plot.species('spur', 'WHL22.738092.0', add=F,flank=10000, 'Pcdh10')
plot.species('spur', 'WHL22.169355.0', add=F,flank=10000, 'Dach2')
plot.species('spur', 'WHL22.606237.0', add=F,flank=10000, 'Notch4')
plot.species('spur', 'WHL22.78013.0', add=F,flank=10000, 'GATA4')
plot.species('spur', 'WHL22.119881.0', add=F,flank=10000, 'MSX2')
plot.species('spur', 'WHL22.754673.3', add=F,flank=10000, 'Rgmb')
plot.species('spur', 'WHL22.323545.0', add=F,flank=10000, 'Fgfrl1')
plot.species('spur', 'WHL22.417635.0', add=F,flank=10000, 'Calm5')
plot.species('spur', 'WHL22.438846.0', add=F,flank=10000, 'Bmp2')
plot.species('spur', 'WHL22.614286.0', add=F,flank=10000, 'Tgif2')
plot.species('spur', 'WHL22.735232.0', add=F,flank=10000, 'Wnt16')
plot.species('spur', 'WHL22.80978.1', add=F,flank=10000, 'Notum')
plot.species('spur', 'WHL22.42488.0', add=F,flank=10000, 'Fzd8')
plot.species('spur', 'WHL22.128311.0', add=F,flank=10000, 'Ascl1')
plot.species('spur', 'WHL22.45238.0', add=F,flank=10000, 'Fzd2')
plot.species('spur', 'WHL22.343932.0', add=F,flank=10000, 'Fgf16')
plot.species('spur', 'WHL22.532435.2', add=F,flank=10000, 'Otx2')
plot.species('spur', 'WHL22.738139.0', add=F,flank=10000, 'POU4f1')
plot.species('spur', 'WHL22.738264.5', add=F,flank=10000, 'Spata18')
plot.species('spur', 'WHL22.273041.0', add=F,flank=10000, 'CYP2C8')
plot.species('spur', 'WHL22.306197.0', add=F,flank=10000, 'Tenm4')
plot.species('spur', 'WHL22.576126.0', add=F,flank=10000, 'TEXT')
plot.species('spur', 'WHL22.423696.0', add=F,flank=10000, 'text2')
plot.species('spur', 'WHL22.647633.0', add=F,flank=10000, 'Chrna2')
plot.species('spur', 'WHL22.52081.0', add=F,flank=50000, 'Wnt5')
plot.species('spur', 'WHL22.735232.0', add=F,flank=10000, 'Wnt16')
plot.species('spur', 'WHL22.596775.0', add=F,flank=10000, 'Wnt-9a')

# END ---------------------------------------------------------------------
