# rm(list = ls())

# Load packages -----------------------------------------------------------

pacman::p_load(ggpubr, magrittr,topGO,ggplot2,cowplot,grid,gridExtra, patchwork)


# Working directory -------------------------------------------------------

# Set working directory
setwd("~/evo/Bra/processing/GO")

# setwd("/Users/dnyansagar/mnt/Bra/processing/GO")

# GO ontology function 

go.analysis <- function(gene2go.file, braGenes.file, allGenes.file, ontology, identifier, highlight){
  gene2go <-    readMappings(gene2go.file)
  allgenes <-   readLines(allGenes.file)
  braGenes <-   readLines(braGenes.file)
  geneList <-   factor(as.integer(allgenes %in% braGenes))
  names(geneList) <- allgenes
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = gene2go, geneSel = braGenes)
  resultFis.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  allRes <- GenTable(GOdata, p.value = resultFis.elim, 
                     orderBy = "p.value", ranksOf = "classic", topNodes = 50, numChar = 100)
  allRes.filtered <- allRes[allRes$Expected >= 1 & 
                              as.numeric(sub(',','.', allRes$p.value, fixed = TRUE)) < 0.01 &
                              ! allRes$GO.ID %in% GO_ignore,]
  allRes.filtered$score <- -log10(as.numeric(sub(',','.', allRes.filtered$p.value, 
                                                fixed = TRUE)))
  allRes.filtered$Term <- factor(allRes.filtered$Term,
                          levels=allRes.filtered$Term[seq(from=length(allRes.filtered$Term),to=1)])
  levels(allRes.filtered$Term) <- sub("homophilic cell adhesion via plasma membrane adhesion molecules", "Cadherins",
                                      levels(allRes.filtered$Term))
  allRes.filtered$Enrichment <- allRes.filtered$Significant/allRes.filtered$Expected
  allRes.filtered$highlight  <- as.factor(ifelse(is.element(allRes.filtered$Term, highlight), 
                                       "Developmental", "Others"))
  return(na.omit(allRes.filtered))
}


# Different representations -----------------------------------------------


# Simple horizontal barplot 
limits = c(100, 600)
plotGO  <- function(df, title, clr){
  #df$colo <- ifelse(df$highlight=="Developmental", "**", "")
  p <- ggplot(df, aes(x=Term, y=score)) + geom_bar(stat='identity', fill=clr) + #
    theme(axis.text.y = element_text(size=16,face="plain"))+
    theme(legend.title=element_blank(), axis.title.y=element_blank(),
          axis.text.y = element_blank()) +
    theme(axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
    labs(title=title)+ coord_flip() + theme(legend.position="none") +
    #annotate("text", x = 1, y = 22, label = "(Genes)") +
    geom_text(aes(y=0.1, label=paste0(Term,' (', Significant,')'),
                  hjust = 0,  size = 4, colour = highlight)) + 
    #, colo, colour = highlight))+ 
    scale_color_manual(values = c("black", "grey50")) +
    ylab(expression(paste("-log" [10], italic(" p"), "-value")))
  
} 

# Lollipop chart
plotGO2 <- function(df , title, clr, limits, breaks){
  ggplot(df, aes(x=Term, y=score)) +
    geom_segment(aes(x=Term, xend=Term, y=0, yend=score),color=clr) +
    geom_point(aes(size= Significant), color = clr, alpha = 0.99) +
    scale_size_continuous(limits=limits , breaks=breaks)+
    coord_flip() + theme_light() + 
    labs(title=title,  y="Score") +
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    geom_text(aes(y=0.1, label=Term , hjust = 0,  size = 400), 
              show.legend = FALSE, nudge_x = 0.4) + 
    theme( legend.title=element_blank(), 
      axis.title.y=element_blank(), 
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank() 
    ) 
}

# Dotplot representation

plotGO3 <- function(df, title){
  require(hrbrthemes)
  ggplot(df, aes(x=Term, y=1)) +
    geom_point(aes(size= Significant,  color = score)) +
    coord_flip() + theme_cowplot()+  
    labs(title=title) + 
    scale_color_gradient(low="blue", high="red")+
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + 
    theme(text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())
}

# Other info
GO_ignore <- c('GO:0008150','GO:0009987')
h1 <- read.table('highlight.txt', header = F, sep = '\n')
highlight <- as.vector(h1$V1)

#############################
##### MOUSE GO ANALYSIS #####
#############################
mmu_gene2go <- 'mmu_goSlim.txt'
mmu_allgenes <- "mmu_all_genes.list"

# Mouse Perturbation ------------------------------------------------------

mmu_bragenes <- "mmu_bra_de_dec.list"
mmu.de.bra <- go.analysis(mmu_gene2go, mmu_bragenes, mmu_allgenes, "BP", "mmu_de", highlight)

mmu_de_plot    <- plotGO2(mmu.de.bra[order(mmu.de.bra$score),],"M. musculus", "#e41a1c",
                          c(0,160), c(40, 80, 120, 160))

mmu_de_plot3    <- plotGO3(mmu.de.bra[order(mmu.de.bra$score),],"M. musculus")

str(mmu.de.bra$Significant)
pretty(range(mmu.de.bra$Significant))
mmu_up <- "mmu_upDec.list"
mmu.de.up <- go.analysis(mmu_gene2go, mmu_up, mmu_allgenes, "BP", "mmu_up", highlight)
mmu.de.up.plot   <- plotGO2(mmu.de.up[order(mmu.de.up$score),],"GO: M. musculus upregulated genes", "#e41a1c")

mmu_down <- "mmu_DownDec.list"
mmu.de.down <- go.analysis(mmu_gene2go, mmu_down, mmu_allgenes, "BP", "mmu_down", highlight)
mmu.de.down.plot <- plotGO2(mmu.de.down[order(mmu.de.down$score),], "GO: M. musculus downregulated genes", "#e41a1c")


# Mouse ChIP targets ------------------------------------------------------

mmu_chip <- '~/evo/Bra/processing/02_select_target/mmutarget2018-12-31.txt'
mmu.chip.bra <- go.analysis(mmu_gene2go, mmu_chip, mmu_allgenes, "BP", "mmu_chip", highlight)
mmu.chip.plot  <- plotGO(mmu.chip.bra[order(mmu.chip.bra$score),], "M. musculus", "#e41a1c")
mmu.chip.plot2 <- plotGO2(mmu.chip.bra[order(mmu.chip.bra$score),], "M. musculus", "#e41a1c",
                          c(100, 800), c(100, 200, 400, 600, 800))
mmu.chip.plot2
mmu.chip.plot3  <- plotGO3(mmu.chip.bra[order(mmu.chip.bra$score),], "M. musculus")

###############################
##### XENOPUS GO ANALYSIS #####
###############################

xtr_gene2go <- "xtr9_slim.txt"
xtr_allgenes <- "xtr_all_genes2.list"

# Xenopus perturbation ----------------------------------------------------

xtr_bra_new <- 'xtr_bra_de.list'
xtr.new.bra <- go.analysis(xtr_gene2go, xtr_bra_new, xtr_allgenes, "BP", "xtr_de", highlight)
xtr.de.fil <- xtr.new.bra[xtr.new.bra$Significant >= 10, ]
#xtr_de_plot   <- plotGO(xtr.de.fil[order(xtr.de.fil$score),],"GO: X.tropicalis DE genes", "orange3")
xtr_de_plot   <- plotGO2(xtr.de.fil[order(xtr.de.fil$score),],"X.tropicalis", "#377eb8",
                         c(0,160), c(40, 80, 120, 160))

xtr_de_plot3   <- plotGO3(xtr.new.bra[order(xtr.new.bra$score),],"X.tropicalis")
xtr_de_plot3

xtr_bra_up <- 'xtr_up.list'
xtr.new.up <- go.analysis(xtr_gene2go, xtr_bra_up, xtr_allgenes, "BP", "xtr_up", highlight)
xtr_de_up_plot     <- plotGO2(xtr.new.up[order(xtr.new.up$score),], "GO: X. tropicalis upregulated genes", "#377eb8")

xtr_bra_down <- 'xtr_down.list'
xtr.new.down <- go.analysis(xtr_gene2go, xtr_bra_down, xtr_allgenes, "BP", "xtr_down", highlight)
xtr_de_down_plot   <- plotGO2(xtr.new.down[order(xtr.new.down$score),], "GO: X. tropicalis downregulated genes", "#377eb8")


# Xenopus ChIP targets ----------------------------------------------------

xtr_chipGenes <- '~/evo/Bra/processing/02_select_target/xtrtarget2018-12-31.txt'
xtr.chip.bra <- go.analysis(xtr_gene2go, xtr_chipGenes, xtr_allgenes, "BP", "xtr_chip", highlight)
xtr.chip.bra_slim <- go.analysis(xtr_gene2go, xtr_chipGenes, xtr_allgenes, "BP", "xtr_chip", highlight)
xtr.chip.fil <- xtr.chip.bra[xtr.chip.bra$Significant >= 10, ]
xtr_chip_plot_slim <- plotGO(xtr.chip.bra_slim[order(xtr.chip.bra_slim$score),],"X. tropicalis" , "#377eb8")
xtr_chip_plot_slim <- plotGO2(xtr.chip.bra_slim[order(xtr.chip.bra_slim$score),], "X. tropicalis" , "#377eb8",
                              c(0,200), c(50, 100, 150, 200))
xtr_chip_plot  <- plotGO2(xtr.chip.bra[order(xtr.chip.bra$score),],"X. tropicalis" , "#377eb8")

xtr_chip_plot3  <- plotGO3(xtr.chip.bra[order(xtr.chip.bra$score),],"X. tropicalis")


####################################
##### NEMATOSTELLA GO ANALYSIS #####
####################################
nve_gene2go  <- "nve_gene2go.tab"
nve_allgenes <- "nve_all_genes.list"

# Nematostella perturbation -----------------------------------------------

nve_braGenes <- "nve_bra_de.list"
nve.de.bra <- go.analysis(nve_gene2go,nve_braGenes, nve_allgenes, "BP", "nve_de",highlight )
nve_de_plot   <- plotGO2(nve.de.bra[order(nve.de.bra$score),], "N. vectensis", "#984ea3",
                         c(0,160), c(40, 80, 120, 160))
nve_de_plot3   <- plotGO3(nve.de.bra[order(nve.de.bra$score),], "N. vectensis")

nve.upGenes <- 'nve_up.list'
nve.de.up <- go.analysis(nve_gene2go,nve.upGenes, nve_allgenes, "BP", "nve_up", highlight )
nve_de_up_plot     <- plotGO2(nve.de.up[order(nve.de.up$score),], "GO: N. vectensis upregulated genes", "#984ea3")

nve.downGenes <- 'nve_down.list'
nve.de.down <- go.analysis(nve_gene2go,nve.downGenes, nve_allgenes, "BP", "nve_down", highlight )
nve_de_down_plot   <- plotGO2(nve.de.down[order(nve.de.down$score),], "GO: N. vectensis downregulated genes", "#984ea3")


# Nematostella ChIP targets -----------------------------------------------


nve_chipGenes <- '~/evo/Bra/processing/02_select_target/nvetarget2018-12-31.txt'
nve.chip.bra <- go.analysis(nve_gene2go, nve_chipGenes, nve_allgenes, "BP", "nve_chip", highlight)
nve_chip_plot <- plotGO2(nve.chip.bra[order(nve.chip.bra$score),], "N. vectensis", "#984ea3", 
                         c(0,160), c(40, 80, 120, 160))
nve_chip_plot1 <- plotGO(nve.chip.bra[order(nve.chip.bra$score),], "N. vectensis", "#984ea3")
nve_chip_plot
nve_chip_plot3 <- plotGO3(nve.chip.bra[order(nve.chip.bra$score),], "N. vectensis")

##
t_sox_all <- "nve_t_sox_all.txt"
t_sox_all  <- "t_sox_both_list.txt"
t_sox_all <- "stboth.txt"
nve_t_sox_all <- go.analysis(nve_gene2go, t_sox_all, nve_allgenes, "BP", "nve_bra_sox_both", highlight)
nve_tsoxall_plot <- plotGO2(nve_t_sox_all[order(nve_t_sox_all$score),], "Bra and Sox motifs", "#984ea3")
nve_tsoxall_plot
##
sox_only <- "sox_only_list.txt"
sox_only <- "sonly.txt"
nve_sox_only <- go.analysis(nve_gene2go, sox_only, nve_allgenes, "BP", "sox_only", highlight)
nve_sox_only_plot <- plotGO2(nve_sox_only[order(nve_sox_only$score),], "Sox motif Only", "#984ea3")
##
t_only <- "t_only_list.txt"
t_only <- "tonly.txt"
nve_t_only <- go.analysis(nve_gene2go, t_only, nve_allgenes, "BP", "Bra_only", highlight)
nve_t_only_plot <- plotGO2(nve_t_only[order(nve_t_only$score),], "Bra motif Only", "#984ea3")
nve_t_only_plot

####################################
###### SEA URCHIN GO ANALYSIS ######
####################################
spu_gene2go <- "spur_gene2go.tab"
spu_allgenes <- "spu_all_genes.list"

# Sea Urchin perturbation -------------------------------------------------

spu_de <- "spur_de.list"
spu.bra.de <- go.analysis(spu_gene2go, spu_de, spu_allgenes, "BP", "spu_de", highlight)
spu.de.plot <- plotGO2(spu.bra.de[order(spu.bra.de$score),], "S. purpuratus", "#4daf4a",
                       c(0,100), c(20, 40, 60, 80, 100))

spu.de.plot3 <- plotGO3(spu.bra.de[order(spu.bra.de$score),], "S. purpuratus")

spu_up <- "spur_up.list"
spu.bra.up <- go.analysis(spu_gene2go, spu_up, spu_allgenes, "BP", "spu_up", highlight)
spu.up.plot <- plotGO2(spu.bra.up, "GO: S. purpuratus upregulated genes", "#4daf4a")

spu_down <- "spur_down.list"
spu.bra.down <- go.analysis(spu_gene2go, spu_down, spu_allgenes, "BP", "spu_down", highlight)
spu.down.plot <- plotGO2(spu.bra.down, "GO: S. purpuratus downregulated genes", "#4daf4a")


# Sea urchin ChIP targets -------------------------------------------------


spu.chip <- 'spu.chip.first.list'
spu.bra.chip <- go.analysis(spu_gene2go, spu.chip, spu_allgenes, "BP", "spu_chip", highlight)
spu_chip_plot <- plotGO2(spu.bra.chip[order(spu.bra.chip$score),], "S. purpuratus", "#4daf4a", 
                         c(0,80), c(20,40,60,80))
spu_chip_plot3 <- plotGO3(spu.bra.chip[order(spu.bra.chip$score),], "S. purpuratus")
spu_chip_plot3


# CAPSASPORA GO ANALYSIS --------------------------------------------------


####################################
###### CAPSASPORA GO ANALYSIS ######
####################################
cow_gene2go <- "cowc_gene2go2.tab"
cow_allgenes <- "cow_list"
cow.braGenes <- "bra_ATAC.txt"
cow.bra.ch <- go.analysis(cow_gene2go, cow.braGenes, cow_allgenes, "BP", "cow_chip", highlight)
cow.bra.chip.plot <- plotGO2(cow.bra.ch, "C. owczarzaki *", "#ff7f00",
                             c(0,60), c(20, 40, 60))
cow.bra.chip.plot3 <- plotGO3(cow.bra.ch, "C. owczarzaki *")



####################################
########## Arrange plots ###########
####################################
grpI <- 'groupI.list'
cow_grpI <- go.analysis(cow_gene2go, grpI, cow_allgenes, "BP", "grpI_chip", highlight)
grpI_plot <- plotGO2(cow_grpI, "Group I genes", "#ff7f00")
# Save it in 1080 X 600 resolution 


# Grid arrange ------------------------------------------------------------


grid.arrange(mmu.chip.plot, xtr_chip_plot_slim,  
             nrow = 1, ncol=2, top = textGrob("GO term analysis of ChIP targets",gp=gpar(fontsize=22,font=1)))
grid.arrange(mmu.chip.plot, xtr_chip_plot_slim, spu_chip_plot,nve_chip_plot, 
             nrow = 1, ncol=2 )#, top = textGrob("Gene Ontology (GO) term analysis",gp=gpar(fontsize=22,font=1)))
grid.arrange(cow.bra.chip.plot, nve_chip_plot, spu_chip_plot, xtr_chip_plot_slim, mmu.chip.plot, nrow = 2, ncol=3)


grid.arrange(mmu_de_plot,xtr_de_plot, spu.de.plot, nve_de_plot,  nrow = 2, ncol=2, 
             top = textGrob("Gene Ontology (GO) term analysis",gp=gpar(fontsize=22,font=1)))



grid.arrange(nve_chip_plot,nve_de_plot,  nrow = 2, ncol=1, 
             top = textGrob("Gene Ontology (GO) term analysis",gp=gpar(fontsize=22,font=1)))
grid.arrange(xtr_chip_plot, xtr_de_plot , nrow = 2, ncol=1, 
             top = textGrob("Gene Ontology (GO) term analysis",gp=gpar(fontsize=22,font=1)))
grid.arrange(nve_chip_plot, nve_de_plot,
             xtr_chip_plot, xtr_de_plot,
             mmu.chip.plot ,mmu_de_plot , nrow = 3, ncol=2, 
             top = textGrob("Gene Ontology (GO) term analysis",gp=gpar(fontsize=22,font=1)))

grid.arrange(nve_chip_plot, nve_de_plot,xtr_chip_plot, xtr_de_plot,nrow = 2, ncol=2, 
             top = textGrob("Gene Ontology (GO) term analysis",gp=gpar(fontsize=22,font=1)))


mmu.chip.bra$species <- 'mmu'
xtr.chip.bra_slim$species <- 'xtr'
nve.chip.bra$species <- 'nve'
spu.bra.chip$species <- 'spu'

ab <- merge(merge(merge(mmu.chip.bra, 
            xtr.chip.bra_slim, by="GO.ID",  all=TRUE),
            spu.bra.chip, by="GO.ID",  all=TRUE),
            nve.chip.bra, by="GO.ID",  all=TRUE)

ab.melt <- reshape2::melt(ab, id.var = 'GO.ID')

ggplot(ab.melt, aes(x=Term, y=score)) +
  geom_segment(aes(x=Term, xend=Term, y=0, yend=score,  color=highlight)) +
  geom_point(aes(size= Significant, color = highlight),  alpha = 0.9) +
  coord_flip() + theme_light() + 
  labs(title=title,  y="Score") +
  scale_color_manual(values = c("red", "grey10")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(y=0.1, label=Term , hjust = 0,  size = 99),nudge_x = 0.4)+ 
  theme(
    legend.title=element_blank(), 
    axis.title.y=element_blank(), 
    axis.text.y = element_blank()) + 
  #theme(legend.position="none") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()) + facet_grid(. ~ variable)


# Patchwork arrange -------------------------------------------------------


# ChIP plots 
mmu.chip.plot2 |xtr_chip_plot_slim / spu_chip_plot |nve_chip_plot | cow.bra.chip.plot
cow.bra.chip.plot /nve_chip_plot | spu_chip_plot /xtr_chip_plot_slim |mmu.chip.plot2 
mmu.chip.plot2 | xtr_chip_plot_slim / spu_chip_plot |nve_chip_plot /cow.bra.chip.plot

mmu.chip.plot3 | xtr_chip_plot3 / spu_chip_plot3 |nve_chip_plot3 /cow.bra.chip.plot3
# ChIP patchwork
(mmu.chip.plot3 | xtr_chip_plot3) /  ( nve_chip_plot3 | spu_chip_plot3/cow.bra.chip.plot3) +
  plot_layout(widths = c(2, 1))

# DE plots 
nve_de_plot3 / spu.de.plot3 | xtr_de_plot3 |mmu_de_plot3

# DE patchwork
(mmu_de_plot3 |xtr_de_plot3)/ (spu.de.plot3 |nve_de_plot3)

mmu.chip.plot3 |mmu_de_plot3
xtr_chip_plot3 |xtr_de_plot3
spu_chip_plot3 |spu.de.plot3
nve_chip_plot3 |nve_de_plot3

scale_color_ipsum()
library(ggplot2)
xtr_de_plot3   <- plotGO3(xtr.de.fil[order(xtr.de.fil$score),],"X.tropicalis")
mmu_de_plot3    <- plotGO3(mmu.de.bra[order(mmu.de.bra$score),],"M. musculus")
nve_de_plot3   <- plotGO3(nve.de.bra[order(nve.de.bra$score),], "N. vectensis")
spu.de.plot3 <- plotGO3(spu.bra.de[order(spu.bra.de$score),], "S. purpuratus")

#####################################################################################
#####################################################################################