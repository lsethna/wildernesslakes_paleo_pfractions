
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readxl)


######Load in data######

#pull data from googledrive
library(googledrive)
install.packages("googlesheets4")
library(googlesheets4)


#authenticates getting your drive
googledrive::drive_auth()

#load data in from your google drive
#link to folder
bsi <- googledrive::as_id("https://drive.google.com/drive/folders/1T1OwrQpbi3JoLwBleRM6tquK_7rX6rYP")

# Identify needed data in the folder
wanted_files <- googledrive::drive_ls(path = bsi) %>%
  # Filter to only needed files
  dplyr::filter(name %in% c("WL_BSi_all.xlsx"))

# Create a folder to download data into
dir.create(path = file.path("raw_data"), showWarnings = F)

# Download that data
purrr::walk2(.x = wanted_files$name, .y = wanted_files$id,
             .f = ~ googledrive::drive_download(file = .y, overwrite = T,
                                                path = file.path("raw_data", .x)))

#import downloaded file from local drive
getwd()
bsi_dat <- read_excel("raw_data/WL_BSi_all.xlsx")
glimpse(bsi_dat)

#######Didn't work, keep trying to make it work######
library(readr)
bsi <- read_csv("BSi_for_R.csv")
glimpse(bsi)



#Make both lists so you can plot them together using cowplot
#make plot
wt_perc <- 
bsi %>%
  ggplot(aes(y=wt_percent,x=date))+
  geom_line()+
  coord_flip()+
  theme_bw()+
  ylab("wt % BSi")+
  facet_grid(rows=vars(lake),scales="free")+
  theme(strip.placement="outside")


#make plot
flux_si <- 
bsi %>%
  ggplot(aes(y=flux_mg,x=date))+
  geom_line()+
  coord_flip()+
  theme_bw()+
  facet_grid(rows=vars(lake),scales="free")+
  ylab("Flux SiO2(mg/cm2 yr)")+
  theme(strip.placement="outside")


library(cowplot)

cowplot::plot_grid(wt_perc,flux_si)

#Just change "lake" and make new graph

bsi %>% filter(lake=="dunnigan") %>% 
  ggplot(aes(y=flux_mg,x=date))+
  geom_line()+
  coord_flip()+
  theme_bw()+
  facet_grid(rows=vars(lake),scales="free")+
  ylab("Flux SiO2(mg/cm2 yr)")+
  theme(axis.title.y = element_blank())+
  theme(strip.placement="outside")
 


#DMAR vs Date
library(readxl)
WL_dating <- read_excel("H:/.shortcut-targets-by-id/1r-S5jPa9KzL2cJXTVseFiknRKoMTF35v/Wilderness Lakes/Sediment cores/Dating files/Final WL Dating files/WL dating.xlsx")
view()
glimpse(WL_dating)


#Note lakes are capitalized
WL_dating %>%
  filter(Lake == "Flame") %>%
  ggplot(aes(x=Sediment_dmar, y=Date_Base))+
  geom_point()+
  geom_path()+
  xlab("Sedimentation Rate (g/cm2 yr)")+
  theme(axis.title.y = element_blank())+
  facet_wrap(~Lake)


  
  
