rm(list=ls())

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/P fractions and TP burial")

librarian::shelf(tidyverse)

master_dat <- read.csv("WL_paleo_masterdataset_20Nov2024.csv")
glimpse(master_dat)

pfrac <- master_dat %>% select(lake,year_loess,dmar_loess,depth,ex_p:recalcitrant_o_p)
glimpse(pfrac)

#calculate labile and refractory based on Edlund 2017
pfrac$labile <- pfrac$ex_p+pfrac$na_oh_p
pfrac$refractory <- pfrac$mineral_p+pfrac$recalcitrant_o_p #mineral_p is HCl_P, labile_o_p is organic P?

pfrac %>% select(lake,depth,year_loess,labile,refractory,bd_p) %>% pivot_longer(cols=c(labile,refractory,bd_p)) %>%
  ggplot(aes(x=depth,y=value,fill=name))+
  geom_area(alpha=0.6,size=0.5, color="black")+
  scale_fill_manual(values=c("hotpink","yellow","darkgreen"))+
  theme_bw()+
  geom_vline(xintercept=5)+
  scale_x_reverse()+
  coord_flip()+
  facet_wrap(~lake,scales="free_y")

#ratio between BD-P+labile : refractory
#calculate this ratio for the top 5 cm, look at the difference between ratio between mixing regimes

pfrac$labile_refrac_ratio <- (pfrac$bd_p+pfrac$labile_o_p+pfrac$ex_p)/pfrac$refractory
pfrac$bd_TP_ratio <- (pfrac$bd_p)/pfrac$tp_results_mg_p_g
pfrac$bd_labile_ratio <- pfrac$bd_p/(pfrac$ex_p+pfrac$labile_o_p)

glimpse(pfrac)

pfrac$lake <- factor(pfrac$lake,
                        levels=c("burnt","smoke","elbow","etwin",
                                 "flame","wtwin",
                                 "dunnigan","finger"))

#(pfrac$bd_p+pfrac$labile_o_p+pfrac$ex_p)/pfrac$refractory
pfrac %>% filter(depth<5) %>% 
  group_by(lake) %>%
  summarize(labile_refrac_ratio = mean(labile_refrac_ratio,na.rm=T)) %>%
ggplot(aes(x=lake,y=labile_refrac_ratio))+
  geom_point(size=3)+
  theme_bw(base_size=16)

#change in productivity
