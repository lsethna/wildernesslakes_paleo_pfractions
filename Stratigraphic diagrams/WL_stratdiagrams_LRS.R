## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, googledrive,readxl,rioja)

getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Download data ----
## ----------------------------------- ##
# Identify URLs of Drive folders with needed data
master_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/13rrJv6QRCx4q342zEjsmbiO1J0xZuyzk")

# Identify needed data in the Drive
wanted_files <- googledrive::drive_ls(path = master_url) %>%
  # Filter to only needed files
  dplyr::filter(name %in% c("WL_paleo_masterdataset_25Sept2024.csv",
                            "interpolated_master_dat_15Oct24.csv"))
# Check those files
wanted_files
# Create a folder to download data into
dir.create(path = file.path("raw_data"), showWarnings = F)
# Download that data
purrr::walk2(.x = wanted_files$name, .y = wanted_files$id,
             .f = ~ googledrive::drive_download(file = .y, overwrite = T,
                                                path = file.path("raw_data", .x)))

## ----------------------------------- ##
# Read in data ----
## ----------------------------------- ##
master_dat <- read.csv("raw_data/WL_paleo_masterdataset_25Sept2024.csv") %>% select(!X)
#check
glimpse(master_dat)
#interpolated master dat
int_master_dat <- read.csv("raw_data/interpolated_master_dat_15Oct24.csv")
glimpse(int_master_dat)

## ----------------------------------- ##
# Select data for strat plots ----
## ----------------------------------- ##

#Worne isotope proposal
worne <- master_dat %>% select(c(lake,year_loess,tp_results_mg_p_g,sio2_wt_percent,chl_a,echine,cantha,fra_crotonensis,ast_formosa))
#WRC presentation
wrc <- master_dat %>% select(c(lake,year_loess,percent_organic,
                               tp_results_mg_p_g,bd_p,
                               sio2_wt_percent,
                               chl_a,echine,cantha,
                               fra_crotonensis,ast_formosa))
#plot using rioja?
librarian::shelf(vegan)

#Burnt
burnt <- wrc %>%
  filter(lake=="burnt")
glimpse(burnt)
burnt.plot <- strat.plot(burnt[1:27,3:11],yvar=burnt$year_loess[1:27],
                         y.tks=seq(1860,2020,20),
                         plot.poly=T,plot.bar=T,col.bar="black",
                         srt.xlabel=45,title="Burnt")
#wtwin
flame <- wrc %>%
  filter(lake=="flame")
glimpse(flame)
flame.plot <- strat.plot(flame[1:21,3:11],yvar=flame$year_loess[1:21],
                         y.tks=seq(1850,2015,15),
                         plot.poly=T,plot.bar=T,col.bar="black",
                         srt.xlabel=45,title="Flame")

## ------------------------------------------------------------- ##
## ---- Selected variables from PCoA in strat plot + CONISS ---- ##
## ------------------------------------------------------------- ##

#finger stratigraph
finger <- int_master_dat %>%
  filter(lake=="finger") %>%
  select(c(year_loess,percent_organic,
           tp_results_mg_p_g,bd_p,
           concentration_si_o2,
           chl_a,ast_formosa,allo,fuco))
glimpse(finger)
finger <- strat.plot(finger[,2:9],yvar=finger$year_loess,
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
#finger=3
addClustZone(finger,chclust.obj,nZone=3,lwd=1.5,lty="dashed",col="gray25")

#wtwin stratigraph
wtwin <- int_master_dat %>%
  filter(lake=="wtwin") %>%
  select(c(year_loess,percent_organic,
           tp_results_mg_p_g,bd_p,
           concentration_si_o2,
           chl_a,ast_formosa,allo,fuco))
glimpse(wtwin)
wtwin <- strat.plot(wtwin[,2:9],yvar=wtwin$year_loess,
                     y.tks=seq(1860,2020,20),
                     plot.poly=T,plot.bar=T,col.bar="black",
                     srt.xlabel=30,title="wtwin")
#wtwin CONISS
#distance matrix  
dist.mat <- vegdist(wtwin[,2:9],method="euclidian", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
#coniss cluster
chclust.obj <- chclust(dist.mat,method="coniss")
#find optimal number of clusters
bstick(chclust.obj)
#wtwin=2
addClustZone(wtwin,chclust.obj,nZone=2,lwd=1.5,lty="dashed",col="gray25")

#burnt stratigraph
burnt <- int_master_dat %>%
  filter(lake=="burnt") %>%
  select(c(year_loess,percent_organic,
           tp_results_mg_p_g,bd_p,
           concentration_si_o2,
           tab_flocculosa,fra_crotonensis,
           chl_a,echine,allo,diato,lutein))
glimpse(burnt)
burnt.plot <- strat.plot(burnt[,2:12],yvar=burnt$year_loess,
                    y.tks=seq(1860,2020,20),
                    plot.poly=T,plot.bar=T,col.bar="black",
                    srt.xlabel=30,las=2,title="Burnt")
#burnt CONISS
#distance matrix  
dist.mat <- vegdist(burnt[,2:12],method="euclidian", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
#coniss cluster
chclust.obj <- chclust(dist.mat,method="coniss")
#find optimal number of clusters
bstick(chclust.obj)
#burnt=2
addClustZone(burnt.plot,chclust.obj,nZone=2,lwd=1.5,lty="dashed",col="gray25")
#clear plots between sites
dev.off()
