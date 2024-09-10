
library(janitor)
library(ggplot2)
library(tidyverse)
library(dplyr)


######Load in data######

#pull data from googledrive
library(googledrive)
install.packages("googlesheets4")
library(googlesheets4)


#authenticates getting your drive
googledrive::drive_auth()

#load data in from your google drive
bsi <- read_sheet("https://drive.google.com/drive/folders/1T1OwrQpbi3JoLwBleRM6tquK_7rX6rYP")

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


  
  
