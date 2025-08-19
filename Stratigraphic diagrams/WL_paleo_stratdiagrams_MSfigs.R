## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse,rioja,tidypaleo)

getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Read in data ----
## ----------------------------------- ##

#master dat
master_dat <- read.csv("Master dataset/WL_paleo_masterdataset_20Nov2024.csv")
glimpse(master_dat)

## ----------------------------------- ##
# Select data for strat plots ----
## ----------------------------------- ##

#select variables for algae strat diagram
paleo_prod <- master_dat %>% select(c(lake,year_loess,
                                      chl_b,fuco,b_car,lutein,diato,allo,myxo,cantha,echine,
                                      aul_subarctica,aul_ambigua,pse_elliptica,ast_formosa,tab_flocculosa,fra_crotonensis,sta_venter))
glimpse(paleo_prod)
unique(paleo_prod$lake)

## ------------------------------- ##
## --- strat plots with ggplot --- ##
## ------------------------------- ##

theme_set(theme_paleo())

paleo_prod %>%
  mutate(lake = factor(lake,levels=c("dunnigan","finger","flame","wtwin","burnt","smoke","elbow","etwin"),
                       labels=c("Dunnigan","Finger","Flame","W. Twin","Burnt","Smoke","Elbow","E. Twin"))) %>%
  filter(!if_all(chl_b:sta_venter,is.na)) %>% #remove rows where all paleo data is NA
  pivot_longer(cols=chl_b:sta_venter) %>%
  mutate(name = factor(name,
                       levels=c("chl_b","fuco","b_car","lut_zea","diato","allo","myxo","cantha","echine",
                                "aul_subarctica","aul_ambigua","pse_elliptica","ast_formosa","tab_flocculosa",
                                "fra_crotonensis","sta_venter"))) %>%
  filter(!is.na(name)) %>% #remove rows where "name" is NA
  ggplot(aes(y=value,x=year_loess))+
  coord_flip()+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  facet_geochem_gridh(vars(name),grouping=vars(lake))+
  theme_classic()+
  theme(strip.text.x=element_text(angle=45,size=10),
        strip.text.y=element_text(size=10),
        strip.background=element_blank(),
        axis.text.x=element_text(angle=90),
        axis.title=element_blank())


#Need to plot lakes separately to control scales

## ------------------------------ ##
## --- strat plots with rioja --- ##
## ------------------------------ ##

#finger stratigraph
finger <- paleo_prod %>%
  filter(lake=="finger")
colnames(finger)
finger_plot <- strat.plot(finger[,3:18],yvar=finger$year_loess,
                          y.tks=seq(1860,2020,20),
                          plot.poly=T,plot.bar=T,col.bar="black",
                          srt.xlabel=30,title="Finger")
#finger CONISS
#distance matrix  
dist.mat <- vegdist(finger[,2:9],method="euclidian", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
#coniss cluster
chclust.obj <- chclust(dist.mat,method="coniss")
#find optimal number of clusters
bstick(chclust.obj)

cutree(chclust.obj,k=3)
#finger=3
addClustZone(finger_plot,chclust.obj,nZone=3,lwd=1.5,lty="dashed",col="gray25")

## ------------------------ ##
## --- cowplot together --- ##
## ------------------------ ##

cowplot::plot_grid(dunnigan_paleo,finger_paleo,
                   flame_paleo,wtwin_paleo,
                   burnt_paleo,smoke_paleo,
                   elbow_paleo,etwin_paleo,
                   nrow=4,labels="AUTO")

