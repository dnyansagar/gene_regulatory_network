# rm(list = ls())

# Load packages -----------------------------------------------------------

pacman::p_load(Gviz, grid, mgcv,Hmisc, data.table,GenomicFeatures)

# Load Files --------------------------------------------------------------

setwd('/export/home/dnyansagar/Bra/BAM')
txdbs <- list(
  nve=makeTxDbFromGFF('/export/home/dnyansagar/Bra/processing/txdb/nveGenes.good.130208.longCDS.gff3'))
peak.files <- list(
  nve=read.table('/export/home/dnyansagar/Bra/processing/txdb/nematostella_Bra_peaks_score_cutoff_7_160108.targets3.txt',
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
  lapply(c('C0','C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'),
         function(x) read.bg(paste0(species,x)))
}
bedGraphData <- sapply(c('nve'),read.all,simplify = F,USE.NAMES = T)

# Functions ---------------------------------------------------------------

get.data.track <- function(bg,species,replicate,chrom, shade, tkname) {
  print(tkname)
  par(las=2)
  DataTrack(
    range = bg[[replicate]][bg[[replicate]]$chromosome == chrom,],
    name=tkname, type = "histogram", fill.histogram = shade, #srt = 45, 
    col.histogram = shade, window=-1, windowSize=1,  srt = 60,
    cex.title=1.5,
    rotation.title=1,
   ylim =c(0,10),   genome = 'species')
  #text(tkname, )
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
  dtrack1 <- get.data.track(bg,species,1,chrom, 'black', 'Input')
  dtrack2 <- get.data.track(bg,species,2,chrom, 'black', 'AB')
  dtrack3 <- get.data.track(bg,species,3,chrom, 'black', 'AB2')
  dtrack4 <- get.data.track(bg,species,4,chrom, 'darkgreen', 'H3K4me1')
  dtrack5 <- get.data.track(bg,species,5,chrom, 'purple3', 'p300')
  dtrack6 <- get.data.track(bg,species,6,chrom, 'magenta3', 'Pol2')
  dtrack7 <- get.data.track(bg,species,7,chrom, 'slategrey', 'H3K27me3')
  dtrack8 <- get.data.track(bg,species,8,chrom, 'turquoise', 'K27Ac1') 
  dtrack9 <- get.data.track(bg,species,9,chrom, 'navy', 'K27Ac') 
  dtrack10 <- get.data.track(bg,species,10,chrom, 'red', 'H3K36me3')
  if (missing(gene)){ name=''} else {name=paste0(' (', gene, ')')}
  ptrack <- AnnotationTrack(start=peak$start,end=peak$end,
                            chromosome = chrom,lineheight=0,min.height=1, 
                            fill="red", cex.title=1.5, rotation.title=1,
                            background.panel="white",name = "Peak")
  gtrack <- GenomeAxisTrack(col="black")
  ot1 <- OverlayTrack(trackList=list(dtrack9, dtrack8), name="OverlayTrack")
  grtrack <- GeneRegionTrack( txdb, chromosome = chrom, # start = start, end = end,
                              fill = "navy", col="navy",
                              showId = TRUE,name = "Gene",  geneSymbol=TRUE,
                              cex.title=1.5,
                              cex.group=1.2,
                              fontface.group=2,
                              rotation=1,
                              rotation.title=1,
                              background.panel="transparent",
                              fontcolor.group="black",
                              col.line="navy"
  )
  # dtrack1, 
  displayPars(grtrack) <- list(size=2)
  displayPars(ptrack) <- list(size=2)
  plotTracks(list(gtrack, 
                  dtrack2, dtrack3, dtrack5, dtrack6,  ot1, dtrack4, dtrack7, dtrack10, 
                  ptrack, grtrack), 
                  from = start, to = end, background.title='transparent',
             col.axis='black', col.title='black', add=add, 
             transformation = function(x) {x/mean(x)},
             main=paste0(capitalize(species), name), margin=50, innerMargin=8
             #littleTicks = TRUE,
             # sizes=c(1.5,3,2,2) #3,
             # sizes=c(1.5,3,3,2,2,2,2,2,2,2,2)
  )}

# Plot genes of interest --------------------------------------------------

ids_short <- c("NVE12349", "NVE5480", "NVE22735", "NVE7113", "NVE21288",
               "NVE10414", "NVE26097","NVE18184", "NVE2323","NVE8194",
               "NVE3111", "NVE4014","NVE12095", "NVE12960", "NVE17595",
               "NVE17746", "NVE963", "NVE4834", "NVE7485","NVE15735",
               "NVE20014", "NVE25316", "NVE9693", "NVE20664", "NVE18683",
               "NVE18678", "NVE18680", "NVE18895", "NVE13981", "NVE25227",
               "NVE19592", "NVE5293", "NVE5433", "NVE5790","NVE7022",
               "NVE11443", "NVE11781", "NVE13735", "NVE19002", "NVE19164",
               "NVE19917", "NVE20578", "NVE20632", "NVE21886", "NVE23698",
               "NVE24657", "NVE3498", "NVE4532", "NVE5908","NVE6048",
               "NVE6690", "NVE7213", "NVE7854", "NVE8062", "NVE8117", 
               "NVE8397","NVE8595", "NVE9269", "NVE9727", "NVE10272", 
               "NVE10336","NVE11414","NVE12324", "NVE13368","NVE18494",
               "NVE18693", "NVE19007", "NVE19228", "NVE19942", "NVE25610",
               'NVE3568', 'NVE1835', 'NVE7119','NVE19736','NVE21992',
               'NVE18028', 'NVE10444')
gene_short <- c("Bmper","NvAdmp-related","Chordin","Ctnnb1","Dickkopf",
                "Dlg2","Dvl2","Fzd4","Vangl2","Prickle2",
                "NvWnt11","NvWnt6","NvWntA","NvWnt1","NvWnt3",
                "NvWnt4","NvBMP2_4","Noggin2","Notum","Rgmb",
                "Tgfbr1l","ZSWIM5_6","Jag1","Hes1","Hes2",
                "HES-2-like","Hes3","Dll4","GLIS","Ptch1",
                "Talpid3","Fgf1l","Vegfc","Map3k9","Rasal2",
                "FGFa1","Mapk4","Igf2bp3","fgfr2i9","FGF8A",
                "Trib2","Rasip1","Rab18","Rapgef2","Map4k4",
                "Rab31","COUP Tf","Adrb2","QRFPR","Mc5r",
                "Zic1","Nova-1","Grid1","Npy2r","HMCN1",
                "Grik4","Slc17a6","Glra2","Olfr319","Dscaml1",
                "NPFFR1L2","Ahnak","Syt17","Sema3a","Ddc",
                "Syt11","Sema6d","Htr5b","Rims2","Nv-ZicA",
                'Bra', 'Fzd10', 'Fzd1', 'Fzd5_8','NvWnt2',
                'Ihh', 'Isl1'
)

path = '../Bra/processing/#figures/screenshots/Screenshots_NoInput/'
for (i in 1:length(ids_short)){
  flname <- paste0(path,'nve', gene_short[i],'.jpg')
  jpeg(flname,  width = 1080, height = 600) 
  plot.species('nve', ids_short[i], add=F, flank=20000, gene_short[i])
  dev.off()
}

# Plot genes of interest manually -----------------------------------------

plot.species('nve','NVE3568',add=F,flank=15000, 'Bra')
plot.species('nve','NVE22735',add=F,flank=10000, 'Chordin')
plot.species('nve','NVE23709',add=F,flank=10000, 'SoxB1')
plot.species('nve','NVE21766',add=F,flank=10000, 'Tbx2/3')
plot.species('nve','NVE3111',add=F,flank=10000, 'NvWnt11')
plot.species('nve','NVE4014',add=F,flank=10000, 'NvWnt6')
plot.species('nve','NVE12095',add=F,flank=10000, 'NvWntA')
plot.species('nve','NVE12960',add=F,flank=10000, 'NvWnt1')
plot.species('nve','NVE17595',add=F,flank=10000, 'NvWnt3')
plot.species('nve','NVE17746',add=F,flank=10000, 'NvWnt4')
plot.species('nve','NVE963',add=F,flank=10000, 'NvBMP2_4')
plot.species('nve','NVE21992',add=F,flank=10000, 'NvWnt2')
plot.species('nve',"NVE18683",add=F,flank=10000,"Hes2")
plot.species('nve','NVE7119',add=F,flank=10000,'Fzd1')
plot.species('nve',"NVE19164", add=F,flank=10000, "FGF8A")


# END ---------------------------------------------------------------------

