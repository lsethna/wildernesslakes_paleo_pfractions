## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, readxl,rioja)

getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Read in data ----
## ----------------------------------- ##
#master data
master_dat <- read.csv("raw_data/WL_paleo_masterdataset_20Nov2024.csv") %>% select(!X)
glimpse(master_dat)
#Organic C variable; keep only lake, date,and org_C variables
org_c <- read.csv("raw_data/deriv_orgC.csv") %>% select(lake,year_loess,orgC_burial_g_cm2_yr,deriv_orgC)
glimpse(org_c)
#pfrac data
pfrac <- read.csv("raw_data/pfrac.csv") %>% pivot_wider(id_cols=c(lake,bottom_interval_cm)) %>%
  rename(depth=bottom_interval_cm)
glimpse(pfrac)

## merge together
master_dat_v2 <- master_dat %>% left_join(org_c) %>% left_join(pfrac)

## ----------------------------------- ##
# Select data for strat plots ----
## ----------------------------------- ##


#Burnt Lake plot
burnt_paleo <- paleo %>%
  filter(lake=="burnt")
burnt.paleo <- strat.plot(burnt_paleo[1:27,7:18],yvar=burnt_paleo$year_loess[1:27], #burnt[rows,col]
                          y.tks=seq(1860,2020,20),
                          plot.poly=T,plot.bar=T,col.bar="black",
                          srt.xlabel=45,title="Burnt")
smoke_paleo <- paleo %>%
  filter(lake=="smoke")
smoke.paleo <- strat.plot(smoke_paleo[1:27,7:18],yvar=smoke_paleo$year_loess[1:27], #burnt[rows,col]
                          y.tks=seq(1860,2020,20),
                          plot.poly=T,plot.bar=T,col.bar="black",
                          srt.xlabel=45,title="Smoke")
etwin_paleo <- paleo %>%
  filter(lake=="etwin")
etwin.paleo <- strat.plot(etwin_paleo[1:71,7:18],yvar=etwin_paleo$year_loess[1:71], #burnt[rows,col]
                          y.tks=seq(1860,2020,20),
                          plot.poly=T,plot.bar=T,col.bar="black",
                          srt.xlabel=45,title="etwin")
#plot using rioja?
librarian::shelf(vegan)

#Burnt
burnt <- wrc %>%
  filter(lake=="burnt")
glimpse(burnt)
burnt.plot <- strat.plot(burnt[1:32,3:17],yvar=burnt$year_loess[1:32], #burnt[rows,col]
                         y.tks=seq(1860,2020,20),
                         plot.poly=T,plot.bar=T,col.bar="black",
                         srt.xlabel=45,title="Burnt")
#flame
flame <- wrc %>%
  filter(lake=="flame")
glimpse(flame)
flame.plot <- strat.plot(flame[1:21,3:11],yvar=flame$year_loess[1:21],
                         y.tks=seq(1850,2015,15),
                         plot.poly=T,plot.bar=T,col.bar="black",
                         srt.xlabel=45,title="Flame")
