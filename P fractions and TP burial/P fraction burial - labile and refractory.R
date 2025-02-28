rm(list=ls())

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/P fractions and TP burial")

librarian::shelf(tidyverse)

master_dat <- read.csv("WL_paleo_masterdataset_20Nov2024.csv")
glimpse(master_dat)

pfrac <- master_dat %>% select(lake,year_loess,dmar_loess,depth,ex_p:recalcitrant_o_p)
glimpse(pfrac)

#calculate labile and refractory based on Edlund 2017
pfrac$labile <- pfrac$ex_p+pfrac$na_oh_p
pfrac$refractory <- pfrac$mineral_p+pfrac$labile_o_p #mineral_p is HCl_P, labile_o_p is organic P?

pfrac %>% select(lake,year_loess,labile,refractory) %>% pivot_longer(cols=c(labile,refractory)) %>%
  ggplot(aes(x=year_loess,y=value,fill=name))+
  geom_area(alpha=0.6,size=0.5, color="black")+
  scale_fill_manual(values=c("yellow","darkgreen"))+
  theme_bw()+
  facet_wrap(~lake,scales="free_y")
