# rm(list = ls())


# Load libraries ----------------------------------------------------------

pacman::p_load(pacman, UpSetR, ggplot2, ggthemes, 
               tidyverse, patchwork, hrbrthemes)


# Functions ---------------------------------------------------------------


co <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#999999")

make.table <- function(filename){
  # Read table
  dtf <- read.table(filename, sep='\t', header=FALSE, row.names=1)
  # Assign column names
  names(dtf) = c('bHLH', 'hmg', 'pax' , 'fox' , 'homeo', 'tbox1', 'tbox2', 
                 'bHLH_d', 'hmg_d', 'pax_d' , 'fox_d' , 'homeo_d', 'tbox1_d', 'tbox2_d', 'extra')
  # Merge T-box1(Half) and T-box2(Full)  for upset intersection
  #print(names(dtf))
  dtf$tbox = ifelse(dtf['tbox1'] == 1 | dtf['tbox2'] == 1, 1, 0)
  # Select T-Box disatnce 
  dtf$tbox_d <- ifelse((is.na(dtf$tbox1_d) & !is.na(dtf$tbox2_d)), dtf$tbox2_d, 
                       ifelse((!is.na(dtf$tbox1_d) & is.na(dtf$tbox2_d)), dtf$tbox1_d, 
                              ifelse((dtf$tbox1_d < dtf$tbox2_d ), dtf$tbox1_d , 
                                     ifelse((dtf$tbox2_d < dtf$tbox1_d) , dtf$tbox2_d, NA))))
  return(dtf)
}
my.boxplot <- function(dataframe, order, col){
  library("ggsci")
  library("reshape2")
  dataframe %>% 
    select(order) %>% stack %>%
    ggplot(aes(x = factor(ind, levels = order), y = values)) +
    geom_boxplot(alpha=0.9, outlier.alpha = 0, width=0.5, fill=col) +  
    labs( x = 'Motif family', y = 'Distance from the summit (bp)')+
    theme_ipsum() + coord_flip() +
    scale_x_discrete(labels = gsub("_d","",order)) +
    theme(text = element_text(size=24),
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24) ,
          axis.title.x = element_text(size = 24),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = -0.45, vjust=2.12)
          )
    # + # +
    # scale_y_continuous(position = "right") +
    # scale_x_discrete(position = "left")
}
my.boxplot2 <- function(dataframe){
  dataframe %>% 
    select(matches('_')) %>% 
    Filter(function(x) !all(is.na(x)), .) %>%
    select(-c("tbox1_d", "tbox2_d"))%>%
    select(sort(current_vars())) %>%
    dplyr::rename_all(
      ~stringr::str_replace_all(., "_d", ""))%>%
    stack() %>%
    ggplot(aes(x = ind, y = values)) +
    geom_boxplot(alpha=0.9, outlier.alpha = 0, width=0.5, fill="royalblue") +  
    labs(x = 'Motif family', y = 'Distance from the summit')+
    theme_bw() +
    theme(text = element_text(size=24),
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24) ,
          #axis.title.x = element_text(color = "navy"),
          axis.title.y = element_text(color = "white")) +
    coord_flip() +
    scale_y_continuous(position = "right") +
    scale_x_discrete(position = "left")
}
my.density <- function(dataframe, order, col){
  library("ggsci")
  library("reshape2")
  dataframe %>% 
    select(order) %>% 
    dplyr::rename_all(
      ~stringr::str_replace_all(., "_d", "")) %>%
    melt() %>% 
    ggplot(aes(x = value, colour=variable)) +
    geom_density(size=1.25) + 
    scale_color_manual(values=col) + theme_ipsum() +
    labs( x = "Distance from Summit (bp)\n",
          y = "Density",
         color = "Motif family\n") +
    theme(legend.text = element_text(size = 24),
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          title = element_text(size = 24))
}

# N. vectensis
u_nve = make.table('~/evo/Bra/motif_analysis/Final_meme-chip-Nve/plots/Upset.txt')
u_nve[,"tbox"] <- sapply(u_nve[,"tbox"], as.numeric)

# X. tropicalis
u_xtr = make.table('~/evo/Bra/motif_analysis/Final_meme-chip-Xtr_new_genome/plots/Upset.txt')
u_xtr[,"tbox"] <- sapply(u_xtr[,"tbox"], as.numeric)

# M. musculus
u_mmu =  make.table('~/evo/Bra/motif_analysis/Final_meme-chip-Mmu_new_paper/plots/Upset.txt')
u_mmu[,"tbox"] <- sapply(u_mmu[,"tbox"], as.numeric)

# S. purpuratus
u_spu =  make.table('~/evo/Bra/motif_analysis/Final_meme-chip-Spu366LowThreshold/plots/Upset.txt')
u_spu[,"tbox"] <- sapply(u_spu[,"tbox"], as.numeric)

u_spu =  make.table('~/evo/Bra/motif_analysis/Spur5//plots/Upset.txt')
u_spu[,"tbox"] <- sapply(u_spu[,"tbox"], as.numeric)


# Generate Plots ----------------------------------------------------------

myboxplot1 <- my.boxplot(u_nve, c('tbox_d','hmg_d','fox_d' ,'homeo_d','bHLH_d','pax_d'),
                         co[c(1, 2, 3, 4, 5, 6)])
myboxplot2 <- my.boxplot(u_xtr, c('tbox_d','homeo_d','pax_d'), 
                         co[c(1, 4, 6)])
myboxplot3 <- my.boxplot(u_mmu, c('tbox_d','hmg_d','homeo_d','fox_d' ,'bHLH_d','pax_d'),
                         co[c(1, 2, 4, 3, 5, 6)])
myboxplot4 <- my.boxplot(u_spu, c('tbox_d','bHLH_d'), 
                         co[c(1, 5)]) 


mydensityplot1 <- my.density(u_nve, c('tbox_d','hmg_d','fox_d' ,'homeo_d','bHLH_d','pax_d'),
                             co)
mydensityplot2 <- my.density(u_xtr, c('tbox_d','homeo_d','pax_d'), 
                             co[c(1, 4, 6)])
mydensityplot3 <- my.density(u_mmu, c('tbox_d','hmg_d','homeo_d','fox_d' ,'bHLH_d','pax_d'),
                             co[c(1, 2, 4, 3, 5, 6)])
mydensityplot4 <- my.density(u_spu, c('tbox_d','bHLH_d'), 
                             co[c(1, 5)])
# Patchwork arrange -------------------------------------------------------


( mydensityplot1 / myboxplot1)
( mydensityplot2 / myboxplot2)
( mydensityplot3 / myboxplot3)
( mydensityplot4 / myboxplot4)


# Upset Plots -------------------------------------------------------------

upset(u_mmu, sets= c('tbox','hmg','homeo','fox' ,'bHLH','pax'),
      sets.bar.color = "royalblue", main.bar.color = "royalblue",
      order.by = "degree", decreasing = "F", point.size =2,
      scale.intersections="identity", text.scale=c(2,2,2,2,2,1.5),
      sets.x.label = "ChIP sequences with motifs",
      mainbar.y.label = "Motif Intersections",
      queries = list(list(query = intersects, params = list('tbox', 'hmg'),   color = "orange", active = T),
                     list(query = intersects, params = list('tbox', 'fox'),  color = "orange", active = T),
                     list(query = intersects, params = list('tbox','homeo'), color = "orange", active = T),
                     list(query = intersects, params = list('tbox','pax'),   color = "orange", active = T),
                     list(query = intersects, params = list('tbox','bHLH'),  color = "orange", active = T)))

upset(u_xtr, sets = c( "tbox" , "pax","homeo"),
      sets.bar.color = "royalblue", main.bar.color = "royalblue",
      order.by = "degree", decreasing = "F", point.size =2,
      scale.intersections="identity", text.scale=c(2,3,3,2,2.5,3),
      sets.x.label = "ChIP sequences with motifs",
      mainbar.y.label = "Motif Intersections",
      queries = list(list(query = intersects, params = list('tbox', 'homeo'),   color = "orange", active = T),
                     list(query = intersects, params = list('tbox', 'pax'),  color = "orange", active = T),
                     list(query = intersects, params = list('tbox', 'pax', 'homeo'),  color = "orange", active = T)))

upset(u_spu, sets = c('bHLH','tbox'),
      sets.bar.color = "royalblue", main.bar.color = "royalblue",
      order.by = "degree", decreasing = "F", point.size =4,
      scale.intersections="identity", text.scale=c(4,4,4,3,4,4),
      sets.x.label = "ChIP sequences with motifs",
      mainbar.y.label = "Motif Intersections", 
      queries = list(list(query = intersects, params = list('tbox','bHLH'),  color = "orange", active = T)))

upset(u_nve, sets=c('tbox','hmg' ,'fox','homeo','pax'), #' 'bHLH'
      sets.bar.color = "royalblue", main.bar.color = "royalblue",
      order.by = "degree", decreasing = "F", point.size =2,
      scale.intersections="identity", text.scale=c(2,2,2,2,2,1.5),
      sets.x.label = "ChIP sequences with motifs",
      mainbar.y.label = "Motif Intersections",
      queries = list(list(query = intersects, params = list('tbox', 'hmg'),   color = "orange", active = T),
                      #list(query = intersects, params = list('tbox', 'fox'),  color = "orange", active = T),
                      list(query = intersects, params = list('tbox','homeo'), color = "orange", active = T)
                      #list(query = intersects, params = list('tbox','pax'),   color = "orange", active = T),
                      #list(query = intersects, params = list('tbox','bHLH'),  color = "orange", active = T)
                     ))



# END ---------------------------------------------------------------------
