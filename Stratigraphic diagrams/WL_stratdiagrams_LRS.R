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
  dplyr::filter(name %in% c("WL_paleo_masterdataset_25Sept2024.csv"))
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

## ----------------------------------- ##
# Select data for strat plots ----
## ----------------------------------- ##

#Worne isotope proposal
worne <- master_dat %>% select(c(lake,year_loess,tp_results_mg_p_g,sio2_wt_percent,chl_a,echine,cantha,fra_crotonensis,ast_formosa))

#plot using rioja?
librarian::shelf(rioja,vegan)

#E. Twin
etwin <- worne %>%
  filter(lake=="etwin")
glimpse(etwin)
etwin.plot <- strat.plot(etwin[,3:9],yvar=etwin$year_loess,
                         y.tks=seq(1860,2020,20),
                         plot.poly=T,plot.bar=T,col.bar="black",
                         srt.xlabel=45,title="Burnt")

#clear plots between sites
dev.off()