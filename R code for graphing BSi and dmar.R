
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readxl)
library(praise)

#for the feels
praise()

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
bsi_dat %>%
  #pretty up lake names and order
  mutate(lake = case_when(lake=="burnt"~"Burnt",
                          lake=="dunnigan"~"Dunnigan",
                          lake=="elbow"~"Elbow",
                          lake=="etwin"~"East Twin",
                          lake=="finger"~"Finger",
                          lake=="flame"~"Flame",
                          lake=="smoke"~"Smoke",
                          lake=="wtwin"~"West Twin")) %>%
  mutate(lake = factor(lake, 
                       levels=c("Dunnigan","Finger","Burnt","Smoke","Elbow","East Twin","Flame","West Twin"),
                       labels=c("(a) Dunnigan","(b) Finger",
                                      "(c) Burnt","(d) Smoke",
                                      "(e) Elbow","(f) East Twin",
                                      "(g) Flame","(h) West Twin"))) %>%
  #plot
  ggplot(aes(x=wt_percent,y=date))+
  geom_path(color="gray40",linewidth=0.75)+
  theme_classic()+
  xlab("Sediment weight % BSi")+
  scale_y_continuous(limits=c(1830,2030),breaks=seq(1830,2030,50))+
  facet_wrap(~lake,scales="free_x",axes="all_y",nrow=4)+
  theme(strip.background=element_blank(),
        strip.text=element_text(hjust=0,size=11),
        axis.title.y=element_blank())

wt_perc

#make plot
flux_si <- 
bsi_dat %>%
  ggplot(aes(y=flux_mg,x=date))+
  geom_line()+
  coord_flip()+
  theme_classic()+
  facet_grid(rows=vars(lake),scales="free")+
  ylab("Flux SiO2(mg/cm2 yr)")+
  facet_wrap(~lake,scales="free",nrow=8)+
  theme(strip.background=element_blank(),
        axis.title.y=element_blank())

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


  
  
