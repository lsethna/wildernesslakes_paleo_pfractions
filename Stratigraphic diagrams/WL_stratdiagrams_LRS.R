## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, googledrive,readxl,rioja)

getwd()
setwd("C:/Users/16512/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Download data ----
## ----------------------------------- ##
# Identify URLs of Drive folders with needed data
master_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/13rrJv6QRCx4q342zEjsmbiO1J0xZuyzk")

# Identify needed data in the Drive
wanted_files <- googledrive::drive_ls(path = master_url) %>%
  # Filter to only needed files
  dplyr::filter(name %in% c("WL_paleo_masterdataset_20Nov2024.csv",
                            "interpolated_master_dat_10Dec24.csv"))
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
master_dat <- read.csv("raw_data/WL_paleo_masterdataset_20Nov2024.csv") #%>% select(!X)
#check
glimpse(master_dat)
#interpolated master dat
int_master_dat <- read.csv("raw_data/interpolated_master_dat_10Dec24.csv")
glimpse(int_master_dat)

## ----------------------------------- ##
# Select data for strat plots ----
## ----------------------------------- ##

#Worne isotope proposal
worne <- master_dat %>% select(c(lake,year_loess,tp_results_mg_p_g,sio2_wt_percent,chl_a,echine,cantha,fra_crotonensis,ast_formosa))
#WRC presentation
wrc <- master_dat %>% select(c(lake,year_loess,percent_organic, #selecting col
                               tp_results_mg_p_g,bd_p,
                               sio2_wt_percent,
                               chl_a,echine,cantha,
                               fra_crotonensis,ast_formosa))

#Paleo paper diagrams
master_dat <- master_dat %>% mutate(recalcitrant_o_p_flux=recalcitrant_o_p*dmar_loess,
                               labile_o_p_flux=labile_o_p*dmar_loess,
                                tp_results_flux=tp_results_mg_p_g*dmar_loess)
  

paleo <- master_dat %>% select(c(lake,year_loess,percent_organic,
                                 recalcitrant_o_p_flux,labile_o_p_flux,tp_results_flux,
                                 ast_formosa, aul_ambigua, fra_crotonensis, lin_intermedia, tab_flocculosa, aul_ambigua, sta_venter,
                                 sio2_wt_percent,
                                 echine,cantha,myxo,chl_a))


custom_colors <- c("skyblue", "salmon", "palegreen", "orange", "purple", "pink", "yellow", "cyan", "red", "blue", "brown", "green", "gold")

color_palette <- c(
  # Dark Gray
  "#4B4B4B",  # Dark Gray
  
  # Blue Shades
  "#AEC6FF",  # Light Blue
  "#5A9BD4",  # Medium Blue
  "#1F4E79",  # Dark Blue
  
  
  "#FFE4B5",  # Very Light Orange (Lightened)
  "#FFD700",  # Light Orange (Golden)
  "#FFA500",  # Medium Orange (Original)
  "#FF7F32",  # Darker Orange (Brighter)
  "#FF6600",  # Dark Orange (Pure Orange)
  "#FF4500",  # Darkest Orange (Vibrant Red-Orange)
  
  #Light gray
  "#d3d3d3", #silica
  
  # Green Shades
  "#D3E4B3",  # Very Light Green (new)
  "#B0E57C",  # Light Green
  "#66B032",  # Medium Green
  "#2E6E1D",  # Dark Green
  
  # Purple Shades
  "#D7BDE2",  # Light Purple
  "#A569BD",  # Medium Purple
  "#6C3483",  # Dark Purple
  
  # Brown Shades
  "#D2B48C",  # Light Brown (Tan)
  "#A0522D",  # Medium Brown
  "#4B2E00"   # Dark Brown
)

# Display the color palette
color_palette


#Burnt Lake plot
burnt_paleo <- paleo %>%
  filter(lake=="burnt")
burnt.paleo <- strat.plot(burnt_paleo[1:27,3:17],yvar=burnt_paleo$year_loess[1:27], #burnt[rows,col]
                          y.tks=seq(1860,2020,20),
                          plot.poly=T,plot.bar=T,col.bar="black",
                          col.poly = color_palette,
                          srt.xlabel=45,title="Burnt")



















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

