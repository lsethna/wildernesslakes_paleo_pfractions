## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, readxl,rioja,vegan)

getwd()
setwd("/Users/lsethna/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

#create date stamp for file exports
date <- format(Sys.Date(),format="%d%b%y")

## ----------------------------------- ##
# Read in data ----
## ----------------------------------- ##
#interpolated master data -- only diatoms and pigments
master_dat <- read.csv("raw_data/interpolated_diatom_pigment_dat_16Nov26.csv") %>% select(!X)
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

## create chla/pheo
master_dat_v2$chla_pheoa <- master_dat_v2$chl_a/master_dat_v2$pheo_a

## get variables from PCoA
imp_variables <- read.csv("Paleo PCoA/sig_variables_pcoa_by_lake_17Jul26.csv")
glimpse(imp_variables)
unique(imp_variables$var)

#get codes to change imp_variables back to match var names in master data
variable_code <- read.csv("Paleo PCoA/pcoa_variable_abbr.csv") %>% select(!X)
glimpse(variable_code)
#add codes for org C and P fractions
variable_code <- variable_code %>%
  add_row(var = c("chla_pheoa","deriv_orgC","labile","ironbound"),
          code = c("chl-a:pheo-a","dC/dt","Labile P","Fe-bound P"))


###change master data to have pretty column names
name_map <- setNames(variable_code$code, variable_code$var)
# rename columns in the dataset
master_dat_v3 <- master_dat_v2 %>%
  rename_with(~ name_map[.x], .cols = names(name_map))
colnames(master_dat_v3)

write.csv(master_dat_v3,file=paste0("Bio proxy CONISS/interp_paleo_var_for_CONISS_",date,".csv"))

## --------------------------------------------------------------------------- ##
## ---------------------------- strat plots ---------------------------------- ##
## --------------------------------------------------------------------------- ##

lakes = unique(master_dat_v2$lake)
pretty_lake_names = c("(c) Burnt","(a) Dunnigan","(f) East Twin","(e) Elbow","(b) Finger","(g) Flame","(d) Smoke","(h) West Twin")
#list to store zone boundary years
zone_boundary_years <- list()
# list to store plots
strat_plot_list <- list()

for (i in 1:length(lakes)){
# Select data for strat plots
vars_code <- imp_variables %>% filter(lake==lakes[i]) %>% select(var) %>% pull()

lake_dat <- master_dat_v3 %>% filter(lake==lakes[i]) %>%
  select(year_loess,all_of(vars_code),`chl-a:pheo-a`,`dC/dt`,`Labile P`,`Fe-bound P`)

#### CONISS ###
#get column numbers to skip
col_names_to_skip <- c("year_loess","chl-a:pheo-a","dC/dt","Labile P","Fe-bound P")
col_n_to_skip <- which(names(lake_dat) %in% col_names_to_skip)
#distance matrix  
dist.mat <- vegdist(lake_dat[,-c(col_n_to_skip)],method="euclidian", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=T)
#coniss cluster
chclust.obj <- rioja::chclust(dist.mat,method="coniss")
#optimal number of clusters
n_cluster <- rioja::bstick(chclust.obj)
n_cluster_lines <- n_cluster[which(n_cluster$dispersion<n_cluster$bstick, arr.ind=TRUE)[1],] %>% select(nGroups) %>% pull()

### pull out the years of the significant zone boundaries ###
#assign each sample to a zone, based on the optimal number of clusters
zone_membership <- cutree(chclust.obj, k = n_cluster_lines)

#lake_dat is in the same row order as zone_membership (both derived from the same
#filtered/ordered data), so a change in zone id between consecutive rows marks a boundary
boundary_idx <- which(diff(zone_membership) != 0)

#take the midpoint year between the two samples flanking each boundary as the "date of change"
boundary_years <- (lake_dat$year_loess[boundary_idx] + lake_dat$year_loess[boundary_idx + 1]) / 2

#save to the list, keyed by lake name
zone_boundary_years[[i]] <- data.frame(lake = lakes[i],
                                       boundary_years = boundary_years)

#create a vector of colors based on variable type
var_colors <- imp_variables %>% filter(lake==lakes[i]) %>% select(type)
var_colors <- var_colors %>% mutate(color = case_when(type=="diatom"~"#919191",
                                                      type=="pigment"~"#d6d6d6"))
strat_colors <- c(#colors for bio vars,
  var_colors$color,
  rep("#f0f0f0",4) #for the preservation index + three geochem vars
  )

strat_plot_list[[i]] <- local({
  dat_i <- lake_dat
  cols_i <- strat_colors
  title_i <- pretty_lake_names[i]
  clust_i <- chclust.obj
  n_lines_i <- n_cluster_lines
  
  # Generate clean year breaks within the dated range of this specific lake
  y_breaks <- pretty(dat_i$year_loess, n = 6)
  # # Keep only breaks above the min of the core so axis doesn't overhang
  y_breaks <- y_breaks[y_breaks >= min(dat_i$year_loess, na.rm = TRUE)]

  cowplot::ggdraw(function() {
    par(oma = c(0, 0, 0, 0))
    
    sp <- strat.plot(
      dat_i[, -1],
      yvar = dat_i$year_loess,
      y.tks = y_breaks,
      plot.poly = TRUE,
      plot.bar = TRUE,
      col.poly = cols_i,
      srt.xlabel = 45,
      title = "",
      # Shorten the stratigraphic body vertically so text/labels get more room
      yTop = 0.64,
      yBottom = 0.14,
      xLeft = 0.14,
      xRight = 0.98,
      las = 2,
      # Tighten label-to-axis distance: c(title, tick_labels, axis_line)
      mgp = c(2.2, 0.5, 0),
      xSpace = 0.02,
      # Larger text sizes relative to a 8.5 x 9.5 in page
      cex.ylabel = 0.75,
      cex.xlabel = 0.75,
      # Smaller x-axis scale numbers and shorter tick marks
      cex.axis = 0.6,
      tcl = -0.25
    )
    addClustZone(sp, clust_i, nZone = n_lines_i, lwd = 1.2, lty = 2, col = "grey25")
  }) +
    cowplot::draw_plot_label(
      label = title_i,
      x = 0.14,          # Matches xLeft
      y = 0.96,          # Sits right above the rotated column headers
      hjust = 0,
      vjust = 1,
      size = 11,
      fontface = "bold"
    )
})

}


## create facet plot
# Reorder list so panels appear (a) through (h) left-to-right, top-to-bottom
# Reorder panels (a) through (h)
panel_order <- order(pretty_lake_names)
ordered_plots <- strat_plot_list[panel_order]

# Use 'scale' (0.90 to 0.98) to adjust whitespace between all 8 panels
panels_grid <- cowplot::plot_grid(
  plotlist = ordered_plots,
  ncol = 2,
  nrow = 4,
  scale = 0.94
)

# Save to PDF 8.5 x 10.5 inch dimensions
cowplot::save_plot(
  filename = paste0("Bio proxy CONISS/Bio proxy CONISS - v4/All_Lakes_strat_coniss_", date, ".pdf"),
  plot = panels_grid,
  base_width = 8.5,
  base_height = 10.5
)

### export zone boundary years
zone_boundary_years_df <- dplyr::bind_rows(zone_boundary_years)
write.csv(zone_boundary_years_df,file=paste0("Bio proxy CONISS/zone_boundary_years_CONISS_",date,".csv"))
