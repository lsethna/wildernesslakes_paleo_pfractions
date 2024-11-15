## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, googledrive,readxl,rioja)

getwd()
#change this to match local GitHub folder
setwd("C:/Users/16512/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

#read in data from folder
library(readr)
WL_paleo<- read_csv("WL_paleo_masterdataset_14Oct2024.csv")


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



wrc <- master_dat %>% select(c(lake,year_loess,percent_organic, #selecting col
                               tp_results_mg_p_g,bd_p,
                               sio2_wt_percent,
                               chl_a,echine,cantha,
                               fra_crotonensis,ast_formosa))

#check normality of data to determine if you need to transform it:

ggplot(master_dat, aes(sample = percent_organic)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = percent_inorg)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = sio2_wt_percent)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = tp_results_mg_p_g)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = labile_o_p)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = recalcitrant_o_p)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = echine)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = cantha)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")

ggplot(master_dat, aes(sample = myxo)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ lake, scales = "free")


#make new data set with col of interest and filter out NA

cpt.meanvar(int_master_dat$percent_organic,penalty="MBIC",pen.value=0,method="AMOC",Q=5,test.stat="Normal",
            class=TRUE,param.estimates=TRUE,shape=1,minseglen=2)
## write a loop ##

lakes <- unique(master_dat$lake)
colnames(master_dat) #what variables do we have?
#list the variables you want to analyze
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



#Change point analysis
library(changepoint)
library(changepoint.np)

data(HC1)





















