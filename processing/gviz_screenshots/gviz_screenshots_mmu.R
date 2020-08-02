# rm(list = ls())

# Load packages -----------------------------------------------------------

pacman::p_load(Gviz, grid, mgcv,Hmisc, data.table,GenomicFeatures)

# Load Files --------------------------------------------------------------

setwd('/home/dnyansagar/Bra/processing/txdb/')
txdbs <- list(
  mmu=makeTxDbFromGFF('/home/dnyansagar/Bra/processing/txdb/mm10_ensembl_compr.protein_coding.supergene.srt.gff3'))
peak.files <- list(
  mmu=read.table('/home/dnyansagar/Bra/processing/txdb/Mmu_T_ChIP.targets2.bed',
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
  lapply(c('C1', 'C2'),function(x) read.bg(paste0(species,x)))
}
bedGraphDataM <- sapply(c('mmu'),read.all,simplify = F,USE.NAMES = T)

# Functions ---------------------------------------------------------------

get.data.track <- function(bg,species,replicate,chrom, shade, tkname) {
  print(tkname)
  DataTrack(
    range = bg[[replicate]][bg[[replicate]]$chromosome == chrom,],
    name=tkname, type = "histogram", fill.histogram = shade,
    col.histogram = shade, window=-1, windowSize=1, cex.axis=.8,
    #cex.title=1, 
    cex.title=3, 
    ylim =c(0,300),   genome = 'species')
}

plot.species <- function(species,txid,flank=1000,add=T, gene) {
  print(gene)
  txdb <- txdbs[[species]]
  peak.file <- peak.files[[species]]
  bg <- bedGraphDataM[[species]]
  txs <- transcripts(txdb,filter=list(tx_name=txid))
  peak <- peak.file[peak.file$GeneID == txid | peak.file$GeneID.1 == txid,]
  start = min(min(peak$start),start(txs))-flank
  end = max(max(peak$end),end(txs))+flank
  chrom = peak$seqnames[1]
  options(ucscChromosomeNames = F)
  
  dtrack1 <- get.data.track(bg,species,1,chrom, 'black', 'Input')
  dtrack2 <- get.data.track(bg,species,2,chrom, 'black', 'AB')
  
  if (missing(gene)){ name=''} else {name=paste0(' (', gene, ')')}
  
  ptrack <- AnnotationTrack(start=peak$start,end=peak$end, col='red',
                            background.panel="white",
                            chromosome = chrom,lineheight=0,min.height=1, 
                            fill="red", cex.title=3,  name = "Peaks")
  
  gtrack <- GenomeAxisTrack(col="black")
  
  grtrack <- GeneRegionTrack( txdb, chromosome = chrom, start = start, end = end,
                              fill = "navy", col="navy",
                              showId = TRUE,name = "Gene",
                              cex.title=3,
                              cex.group=2,
                              fontface.group=2,
                              # stacking(txdb ,dense),
                              background.panel="white",
                              fontcolor.group="black",
                              col.line="navy"
  )
  # dtrack1,
  plotTracks(list(gtrack, 
                   dtrack2,
                  ptrack, grtrack), 
             from = start, to = end, background.title='white',
             col.axis='black', col.title='black', add=add,
             main=paste0(capitalize(species), name),
             #littleTicks = TRUE,
             sizes=c(1.5,5,2,2)
  )
}

save.plot <- function(flname, ht = 0, wd = 0, genename, geneid ){
  jpeg(flname,  width = wd, height = ht)
  plot.species('mmu','ENSMUSG00000038132Supergene',add=F, flank=10000,'rbm24')
  dev.off()
}

# Plot genes of interest --------------------------------------------------

# get screenshots of these genes

path = '../Bra/processing/#figures/screenshots/Screenshots_NoInput/'

ids <- c('ENSMUSG00000038132Supergene','ENSMUSG00000000127Supergene','ENSMUSG00000000263Supergene',
         'ENSMUSG00000050295Supergene', 'ENSMUSG00000038132Supergene','ENSMUSG00000030543Supergene',
         'ENSMUSG00000061524Supergene', 'ENSMUSG00000035799Supergene','ENSMUSG00000057098Supergene', 
         'ENSMUSG00000032035Supergene','ENSMUSG00000023991Supergene','ENSMUSG00000023781Supergene',
         'ENSMUSG00000029844Supergene','ENSMUSG00000060969Supergene','ENSMUSG00000047002Supergene',
         'ENSMUSG00000020647Supergene','ENSMUSG00000063972Supergene','ENSMUSG00000028023Supergene',
         'ENSMUSG00000051367Supergene','ENSMUSG00000021540Supergene','ENSMUSG00000045680Supergene',
         'ENSMUSG00000071757Supergene','ENSMUSG00000030544Supergene','ENSMUSG00000029673Supergene',
         'ENSMUSG00000044813Supergene', 'ENSMUSG00000050700Supergene','ENSMUSG00000044338Supergene',
         'ENSMUSG00000034227Supergene','ENSMUSG00000037025Supergene','ENSMUSG00000019960Supergene',
         'ENSMUSG00000062327Supergene', 'ENSMUSG00000028031Supergene','ENSMUSG00000029671Supergene', 
         'ENSMUSG00000030170Supergene' )

gene <- c('rbm24','unkown1', 'unkown2', 
          'Foxc1', 'Rbm24','Mesp2',
          'Zic2', 'Twist1','Ebf1',
          'Ets1','Foxp4', 'Hes7',
          'Hoxa1','Irx1','Msgn1',
          'Ncoa1','Nr6a1','Pitx2',
          'Six1','Smad5','Tcf21',
          'Zhx2', 'Mesp1','Auts2',
          'Shb','Emilin3','Aplnr',
          'Foxj1','Foxa2','Dusp6', 
          'T', 'Dkk2', 'Wnt16', 'Wnt5b')

for (i in 1:length(ids)){
  flname <- paste0(path,'mmu', gene[i],'fx.jpg')
  jpeg(flname,  width = 1080, height = 600, pointsize = 12, quality = 100)
  plot.species('mmu', ids[i], add=F, flank=25000, gene[i])
  dev.off()
}

# Plot genes of interest manually -----------------------------------------

plot.species('mmu','ENSMUSG00000050700Supergene',add=F, flank=50000,'Emilin3')
plot.species('mmu','ENSMUSG00000050295Supergene',add=F, flank=35000,'Foxc1')
plot.species('mmu','ENSMUSG00000020647Supergene',add=F, flank=25000,'Ncoa1')
plot.species('mmu','ENSMUSG00000062327Supergene',add=F, flank=10000,'T')
plot.species('mmu','ENSMUSG00000044813Supergene',add=F, flank=25000,'Shb')
plot.species('mmu','ENSMUSG00000032035Supergene',add=F, flank=25000,'Ets1')
plot.species('mmu','ENSMUSG00000028031Supergene',add=F, flank=25000,'Dkk2')
plot.species('mmu','ENSMUSG00000029671Supergene',add=F, flank=100000,'Wnt16')
plot.species('mmu','ENSMUSG00000030170Supergene',add=F, flank=100000,'Wnt5b')


# END ---------------------------------------------------------------------


