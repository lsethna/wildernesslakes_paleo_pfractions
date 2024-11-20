## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, googledrive,readxl)

getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Download data ----
## ----------------------------------- ##

# Identify URLs of Drive folders with needed data
main_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1r-S5jPa9KzL2cJXTVseFiknRKoMTF35v") #main WL folder
bsi_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1T1OwrQpbi3JoLwBleRM6tquK_7rX6rYP") #contains WL_BSi_all.xlsx
dates_dmar_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1qzXj9PzjJPqfbgQhzEZGg8Yr2Xr9knv4") #contains sections_interp_year_dmar_30July2024.csv
loi_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1WhQLgIDeu5N1SvcdrYztCWfjBdCb-HR_") #contains WL_LOI_allcores.xlsx
pigment_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1y-zoTHbMW9eBUMPLkWLC_P6Qvcw-n-yL") #raw pigment data, all cores WL_pigments_allcores_14Oct24.csv
pfracs_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/11lWhi9zlPjE2PDnLIAPNbDnbVXrGAxQo") #Pfrac_mass_focuscorrect google sheet
diatoms_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/111Cev1JyY5i4EFqkmUcvtsrBZHt_Z_f8") #wilddiatom_rawdat.csv
isotopes_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1_NSK63Qv61FwfsvhzQi-e1bogk0Fbd5E") #WL_isotopes_rawdat.xlsx

# Identify needed data in the Drive
wanted_files <- googledrive::drive_ls(path = bsi_url) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = dates_dmar_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = loi_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = pigment_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = pfracs_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = diatoms_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = isotopes_url)) %>%
  # Filter to only needed files
  dplyr::filter(name %in% c("WL_BSi_all.xlsx",
                            "sections_interp_year_dmar_30July2024.csv",
                            "WL_LOI_allcores.xlsx",
                            "WL_allpigments_20Nov24.csv",
                            "Pfrac_mass_focuscorrect", #saved as Google sheet, no file extension needed
                            "wilddiatom_rawdat.csv",
                            "WL_isotopes_rawdat.xlsx"
                            ))
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
bsi <- read_excel("raw_data/WL_BSi_all.xlsx")
dates_dmar <- read.csv("raw_data/sections_interp_year_dmar_30July2024.csv")
loi <- read_excel("raw_data/WL_LOI_allcores.xlsx") %>% janitor::clean_names()
pfracs <- read_excel("raw_data/Pfrac_mass_focuscorrect.xlsx")
pigments <- read.csv("raw_data/WL_allpigments_20Nov24.csv")
diatoms <- read.csv("raw_data/wilddiatom_rawdat.csv")
isotopes <- read_excel("raw_data/WL_isotopes_rawdat.xlsx")

#check the data
glimpse(bsi)
glimpse(dates_dmar)
glimpse(loi)
glimpse(pfracs)
glimpse(pigments)
glimpse(diatoms)
glimpse(isotopes)

## ----------------------------------- ##
# Match columns and lake names ----
## ----------------------------------- ##
unique(bsi$lake) #master lake name convention: lowercase, no spaces, abbreviate East and West as e and w
unique(dates_dmar$lake)
unique(loi$lake)
unique(pfracs$lake)
unique(pigments$lake)
unique(diatoms$Lake)
unique(isotopes$Lake)

#format column names to match
#lake, depth, analyte before name...
bsi_v2 <- bsi %>% select(!c(date,dmar,flux_g)) %>% janitor::clean_names() %>%
  rename(flux_sio2_mg=flux_mg,
         sio2_wt_percent=wt_percent) 
glimpse(bsi_v2)

dates_dmar_v2 <- 
  dates_dmar %>% select(lake,year_loess,dmar_loess,base) %>%
  rename(depth=base)%>%
  mutate(lake=case_when(lake=="Burnt"~"burnt",
                        lake=="Dunnigan"~"dunnigan",
                        lake=="East Twin"~"etwin",
                        lake=="Elbow"~"elbow",
                        lake=="Finger"~"finger",
                        lake=="Flame"~"flame",
                        lake=="Smoke"~"smoke",
                        lake=="West Twin"~"wtwin"))

loi_v2 <- loi %>% select(!c(top)) %>%
  rename(depth=base)%>%
  mutate(lake=case_when(lake=="Burnt"~"burnt",
                        lake=="Dunnigan"~"dunnigan",
                        lake=="E Twin"~"etwin",
                        lake=="Elbow"~"elbow",
                        lake=="Finger"~"finger",
                        lake=="Flame"~"flame",
                        lake=="Smoke"~"smoke",
                        lake=="W Twin"~"wtwin"))

pfracs_v2 <- pfracs %>% select(!c(`...1`,id,notes,DMAR,year)) %>%
  rename(depth=depth_base)

pigments_v2 <- pigments %>% select(!c(X)) %>%
  mutate(lake=case_when(lake=="BURNT"~"burnt",
                        lake=="Dunnigan"~"dunnigan",
                        lake=="East Twin"~"etwin",
                        lake=="Elbow"~"elbow",
                        lake=="FINGER"~"finger",
                        lake=="Flame"~"flame",
                        lake=="Smoke"~"smoke",
                        lake=="West Twin"~"wtwin"))

diatoms_v2 <- diatoms %>% select(!c(X,Year)) %>% #use dates from loess model
  dplyr::rename(depth=botepth,
                lake=Lake) %>%
  mutate(lake=case_when(lake=="Burnt"~"burnt",
                        lake=="Dunnigan"~"dunnigan",
                        lake=="ETwin"~"etwin",
                        lake=="Elbow"~"elbow",
                        lake=="Finger"~"finger",
                        lake=="Flame"~"flame",
                        lake=="Smoke"~"smoke",
                        lake=="WTwin"~"wtwin"))

isotopes_v2 <- isotopes %>% select (!c(Year)) %>%
  dplyr::rename(lake=Lake,
                depth=`Depth (cm)`,
                d15N_perc_air=`δ15N ‰ AIR`,
                d13C_org_VPDB=`δ¹³Corg (VPDB)`,
                perc_TOC=`%TOC`,
                perc_TN=`%TN`,
                TOC_TN_ratio=`TOC/TN`) %>%
  mutate(lake=case_when(lake=="Burnt"~"burnt",
                        lake=="Dunn"~"dunnigan",
                        lake=="ETW"~"etwin",
                        lake=="Elbow"~"elbow",
                        lake=="Finger"~"finger",
                        lake=="Flame"~"flame",
                        lake=="Smoke"~"smoke",
                        lake=="WTW"~"wtwin"))

#calculate relative abundance
diatoms_v3 <- diatoms_v2 %>% mutate(tot_count=rowSums(across(`Ach..microcephala`:`Uln..ulna`), na.rm=T)) %>%
  mutate_at(vars(-c(lake,depth,tot_count)),funs(./tot_count*100)) %>%
  mutate(tot_abund=rowSums(across(`Ach..microcephala`:`Uln..ulna`),na.rm=T)) %>%
  janitor::clean_names()
glimpse(diatoms_v3)

## ----------------------------------- ##
# Merge together ----
## ----------------------------------- ##
master_v1 <- dates_dmar_v2 %>% 
  full_join(loi_v2) %>%
  full_join(bsi_v2) %>%
  full_join(pfracs_v2) %>%
  full_join(pigments_v2) %>%
  full_join(diatoms_v3) %>%
  full_join(isotopes_v2)
glimpse(master_v1)

#plot to check funky stuff
#just date v. sample to check for doubles?
ggplot(master_v1,aes(x=year_loess,y=depth))+
  geom_point()+
  facet_wrap(~lake)

getwd()
write.csv(master_v1,file="Master dataset/WL_paleo_masterdataset_20Nov2024.csv")
