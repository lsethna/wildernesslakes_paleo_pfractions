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
pigment_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1y-zoTHbMW9eBUMPLkWLC_P6Qvcw-n-yL") #raw pigment data, all cores WL_pigments_allcores.csv
pfracs_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/11lWhi9zlPjE2PDnLIAPNbDnbVXrGAxQo") #Pfrac_mass_focuscorrect google sheet

# Identify needed data in the Drive
wanted_files <- googledrive::drive_ls(path = bsi_url) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = dates_dmar_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = loi_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = pigment_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = pfracs_url)) %>%
  # Filter to only needed files
  dplyr::filter(name %in% c("WL_BSi_all.xlsx",
                            "sections_interp_year_dmar_30July2024.csv",
                            "WL_LOI_allcores.xlsx",
                            "WL_pigments_allcores.csv",
                            "Pfrac_mass_focuscorrect" #saved as Google sheet, no file extension needed
                            ))

# Check those files
wanted_files

# Create a folder to download data into
dir.create(path = file.path("raw_data"), showWarnings = F)

# Download that data
purrr::walk2(.x = wanted_files$name, .y = wanted_files$id,
             .f = ~ googledrive::drive_download(file = .y, overwrite = T,
                                                path = file.path("raw_data", .x)))
