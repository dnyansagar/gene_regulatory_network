
# Load packages -----------------------------------------------------------


pacman::p_load(UpSetR,tidyverse, data.table)


# Load data ---------------------------------------------------------------


uset <- read.table('/Users/dnyansagar/evo/Bra/processing/03_make_tables/upset2020-07-30.txt', 
                   header = F, sep = '\t')
names(uset) <- c('OMA','C. owczarzaki', 'N. vectensis', 'S. purpuratus', 
                'C. intestinalis', 'X. tropicalis', 'M. musculus')

oma_annot <- read.table('/Users/dnyansagar/evo/Bra/processing/03_make_tables/OMA_ANNO_FILT.txt',
                        header = F, sep = '\t')
names(oma_annot) <- c("OMA","MMU", "XTR", "SPU", "NVE")


# Upset plot --------------------------------------------------------------

upset(uset, 
      sets = c('C. owczarzaki', 'N. vectensis', 'S. purpuratus', 
                'C. intestinalis', 'X. tropicalis', 'M. musculus'), 
      nintersects = NA, keep.order = T, order.by = "freq" , 
      sets.bar.color = "black", main.bar.color = "black",
      mainbar.y.label = "Brachyury target Intersection", 
      queries = list(list(query = intersects,
                          params = list('C. owczarzaki', 'N. vectensis',
                                        'S. purpuratus', 'C. intestinalis',
                                        'X. tropicalis', 'M. musculus'),
                          color = "red",
                          active = T),
                     list(query = intersects,
                          params = list('S. purpuratus', 'C. intestinalis',
                                        'X. tropicalis', 'M. musculus'),
                          color = "limegreen",
                          active = T),
                     list(query = intersects,
                          params = list('C. intestinalis', 'X. tropicalis',
                                        'M. musculus'),
                          color = "blue",
                          active = T),
                     list(query = intersects,
                          params = list('X. tropicalis', 'M. musculus'),
                          color = "coral",
                          active = T),
                     list(query = intersects,
                          params = list('S. purpuratus', 'N. vectensis'),
                          color = "darkorchid",
                          active = T)))
grid.text("Brachyury target intersection",x = 0.65, y=0.95, gp=gpar(fontsize=20))



# Intersection annotation -------------------------------------------------


f <- function(data, annot, ...) {
  take <- select(data, one_of(...)) %>% colnames
  notake <- select(data,-c("OMA"), -one_of(...)) %>% colnames
  data %>%
  filter_at(vars(take),  ~.==1) %>%
  filter_at(vars(notake),  ~.!=1) %>% 
    select(OMA) ->tempdata
  annot %>% 
    filter(OMA %in% tempdata$OMA) %>% 
    select(OMA, MMU, XTR, SPU)
}


vertebrates <- f(uset, oma_annot, "X. tropicalis", "M. musculus")
fwrite(vertebrates, file="~/Documents/projects/Bra/vertebrates.tsv", sep="\t")

chordates   <- f(uset, oma_annot, "C. intestinalis", "X. tropicalis", "M. musculus" )
fwrite(chordates, file="~/Documents/projects/Bra/chordates.tsv", sep="\t")

bilaterans  <- f(uset, oma_annot,"S. purpuratus","C. intestinalis", "X. tropicalis", "M. musculus")
fwrite(bilaterans, file="~/Documents/projects/Bra/bilaterans.tsv", sep="\t")

non_chordates <- f(uset, oma_annot,"N. vectensis", "S. purpuratus")
fwrite(non_chordates, file="~/Documents/projects/Bra/non_chordates.tsv", sep="\t")


non_chordate <- uset %>% 
  filter_at(vars(c("N. vectensis", "S. purpuratus")), ~.==1) %>% select(OMA)
  filter_at(vars(-c("N. vectensis", "S. purpuratus")), ~.!=1) %>% 
  select(OMA)

oma_annot %>% 
  filter(OMA %in% non_chordate$OMA) %>% 
  select(OMA, MMU, XTR, SPU) %>% View()


# END ---------------------------------------------------------------------


