getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions")

librarian::shelf(tidyverse,RColorBrewer)

#get master dataset
master <- read.csv("Master dataset/WL_paleo_masterdataset_20Nov2024.csv")
glimpse(master)

diatom <- read.csv("WL_diatoms_5percabund.csv")
glimpse(diatom)

#filter master dataset to just diatoms in diatom list
diatoms_abund <- master %>% select(lake,year_loess,c(diatom$name)) %>%
  #pivot to get diatom names as variable
  pivot_longer(cols=3:28,names_to="diatom",values_to="rel_abund")
glimpse(diatoms_abund)

#get ecology notes
ecology <- read.csv("WL_diatoms_5percabund_ecology.csv") %>% select(diatom,mixing)
glimpse(ecology)
#add colors groups for P and B
#how many of each?
length(which(ecology$mixing=="P"))
length(which(ecology$mixing=="B"))

color_P <- grDevices::colorRampPalette(c("#0e5f08","#7bdf79"))(12)
color_B <- colorRampPalette(c("#ff5207","#ffcc00"))(14)
color_palettes <- c(color_B,color_P)
#order diatoms by B/P
ecology_colors <- ecology %>% arrange(mixing) %>% mutate(color=color_palettes) 
glimpse(ecology_colors)

#merge with diatoms_abund
diatoms_abund_ecology <- diatoms_abund %>% left_join(ecology_colors) %>%
  arrange(lake,year_loess,mixing)
glimpse(diatoms_abund_ecology)

#stacked area plot
ggplot(diatoms_abund_ecology,aes(y=rel_abund,x=year_loess,fill=diatom))+
  geom_area()+
  scale_fill_manual(values=diatoms_abund_ecology$color)+
  coord_flip()+
  facet_wrap(~lake,scales="free")+
  theme_bw(base_size=15)
