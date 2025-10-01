## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, readxl,rioja,vegan)

getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Read in data ----
## ----------------------------------- ##
#interpolated master data
master_dat <- read.csv("raw_data/interpolated_master_dat_4Sep25.csv") %>% select(!X)
glimpse(master_dat)
#Organic C variable; keep only lake, date,and org_C variables
org_c <- read.csv("raw_data/deriv_orgC.csv") %>% select(lake,year_loess,orgC_burial_g_cm2_yr,deriv_orgC)
glimpse(org_c)
#pfrac data
pfrac <- read.csv("raw_data/pfrac.csv") %>% pivot_wider(id_cols=c(lake,bottom_interval_cm)) %>%
  rename(depth=bottom_interval_cm)
glimpse(pfrac)

## merge together
master_dat_v2 <- master_dat %>% left_join(org_c) %>% left_join(pfrac)
glimpse(master_dat_v2)

## get variables from PCoA
imp_variables <- read.csv("Paleo PCoA/sig_variables_pcoa_by_lake.csv")
#get codes to change imp_variables back to match var names in master data
variable_code <- read.csv("Paleo PCoA/pcoa_variable_abbr.csv")

## --------------------------------------------------------------------------- ##
## ---------------------------- strat plots ---------------------------------- ##
## --------------------------------------------------------------------------- ##

lakes = unique(master_dat_v2$lake)
pretty_lake_names = c("Burnt","Dunnigan","East Twin","Elbow","Finger","Flame","Smoke","West Twin")

for (i in 1:length(lakes)){
# Select data for strat plots
vars_code <- imp_variables %>% filter(lake==lakes[i]) %>% select(var) %>% pull()
#get full name from code df
vars <- variable_code %>% filter(code %in% vars_code) %>% select(var) %>% pull()

lake_dat <- master_dat_v2 %>% filter(lake==lakes[i]) %>%
  select(year_loess,all_of(vars),deriv_orgC,labile,ironbound)

#CONISS
#get column numbers to skip
col_names_to_skip <- c("year_loess","deriv_orgC","labile","ironbound")
col_n_to_skip <- which(names(lake_dat) %in% col_names_to_skip)
#distance matrix  
dist.mat <- vegdist(lake_dat[,-c(col_n_to_skip)],method="euclidian", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
#coniss cluster
chclust.obj <- chclust(dist.mat,method="coniss")

#create a vector of colors based on variable type
var_colors <- imp_variables %>% filter(lake==lakes[i]) %>% select(var_type)
var_colors <- var_colors %>% mutate(color = case_when(var_type=="diatoms"~"#ffa500",
                                                      var_type=="pigment"~"#006400"))
strat_colors <- c(#colors for bio vars,
  var_colors$color,
  rep("hotpink",3) #for the three geochem vars
  )

#will also need to pretty up the lake names for the titles

#run this to save plot
file = paste0("Bio proxy CONISS/",lakes[i],"_strat_coniss.pdf")

pdf(
  file=file,
  width=11,height=8.5
)

strat.plot(lake_dat[,-1],yvar=lake_dat$year_loess, #[rows,col]
           y.tks=seq(plyr::round_any(min(lake_dat$year_loess),10),2020,10),
           plot.poly=T,plot.bar=F,col.poly=strat_colors,
           srt.xlabel=45,title=pretty_lake_names[i],
           las=2,mgp=c(3,1,0.25),xSpace=0.01,
           clust=chclust.obj)

dev.off() #clear before rerunning next lake

}
