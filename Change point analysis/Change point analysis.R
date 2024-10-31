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
#interpolated master dat
int_master_dat <- read.csv("raw_data/interpolated_master_dat_15Oct24.csv")
glimpse(int_master_dat)

## write a loop ##

lakes <- unique(master_dat$lake)
colnames(master_dat) #what variables do we have?
variable <- c("percent_organic","tp_results_mg_p_g","bd_p","sio2_wt_percent") #list the variables you want to analyze
lakes
 
for (i in 1:length(lakes)) {
  d <- subset(master_dat,lake==lakes[i])
  
  for (j in 1:length(variable)) {
     #change point analysis for each variable
    var_d <- subset(d, variable==variable[j])
    change_point_results[[j]] <- analysis_output
  }
  
  #do some analysis
  
  output_list[[i]] <- change_point_results
}



























