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

for (i in 1:length(lakes)) {
  
  # 1. Select data for strat plots
  vars_code <- imp_variables %>% filter(lake == lakes[i]) %>% select(var) %>% pull()
  
  lake_dat <- master_dat_v3 %>% 
    filter(lake == lakes[i]) %>%
    select(year_loess, all_of(vars_code), `chl-a:pheo-a`, `dC/dt`, `Fe-bound P`, `Labile P`)
  
  # 2. CONISS clustering
  col_names_to_skip <- c("year_loess", "chl-a:pheo-a", "dC/dt", "Fe-bound P", "Labile P")
  col_n_to_skip <- which(names(lake_dat) %in% col_names_to_skip)
  
  dist.mat <- vegdist(lake_dat[, -c(col_n_to_skip)], method = "euclidian",
                      binary = FALSE, diag = FALSE, upper = FALSE, na.rm = TRUE)
  chclust.obj <- rioja::chclust(dist.mat, method = "coniss")
  
  n_cluster <- rioja::bstick(chclust.obj)
  n_cluster_lines <- n_cluster[which(n_cluster$dispersion < n_cluster$bstick, arr.ind = TRUE)[1], ] %>%
    select(nGroups) %>% pull()
  
  zone_membership <- cutree(chclust.obj, k = n_cluster_lines)
  boundary_idx <- which(diff(zone_membership) != 0)
  boundary_years <- (lake_dat$year_loess[boundary_idx] + lake_dat$year_loess[boundary_idx + 1]) / 2
  
  zone_boundary_years[[i]] <- data.frame(lake = lakes[i],
                                         boundary_years = boundary_years)
  
  # 3. Map variable types and fill colors
  var_types <- imp_variables %>% 
    filter(lake == lakes[i]) %>% 
    select(var, type) %>%
    bind_rows(tibble(
      var = c("chl-a:pheo-a", "dC/dt", "Fe-bound P", "Labile P"),
      type = rep("geochem", 4)
    ))
  
  # Preserve exact column order from lake_dat
  var_order <- colnames(lake_dat)[-1]
  
  # 4. Pivot longer for ggplot
  lake_long <- lake_dat %>%
    pivot_longer(
      cols = -year_loess,
      names_to = "var",
      values_to = "value"
    ) %>%
    left_join(var_types, by = "var") %>%
    mutate(
      var = factor(var, levels = var_order),
      fill_col = case_when(
        type == "diatom"  ~ "#919191",
        type == "pigment" ~ "#d6d6d6",
        type == "geochem" ~ "#ffffff"
      )
    ) %>%
    arrange(var, year_loess)
  
  # 5. Control y-axis breaks & bounds (min inside core, top capped at 2020)
  y_min <- min(lake_dat$year_loess, na.rm = TRUE)
  y_breaks <- pretty(c(y_min, 2020), n = 5)
  y_breaks <- y_breaks[y_breaks >= y_min & y_breaks <= 2020]
  y_top_limit <- 2020
  
  # 6. Build ggplot stratigraphic diagram
  p <- ggplot(lake_long, aes(y = year_loess)) +
    # Interpolated area polygon (orientation = "y" fills horizontally from 0 to value)
    geom_ribbon(
      aes(xmin = pmin(0, value, na.rm = TRUE),
          xmax = pmax(0, value, na.rm = TRUE),
          fill = fill_col),
      orientation = "y",
      color = "black",
      linewidth = 0.25,
      na.rm = TRUE
    ) +
    # Horizontal sample bars
    geom_segment(
      aes(x = 0, xend = value, yend = year_loess),
      color = "black",
      linewidth = 0.25,
      na.rm = TRUE
    ) +
    # CONISS zone boundary lines
    geom_hline(
      yintercept = boundary_years,
      linetype = "dashed",
      color = "grey25",
      linewidth = 0.45
    ) +
    scale_fill_identity() +
    # Clean x-axis: upper and lower bounds on column
    scale_x_continuous(
      breaks = function(limits) {
        # Keep 0 and the largest pretty break inside the data range
        b <- pretty(limits, n = 2)
        b <- b[b >= limits[1] & b <= limits[2]]
        if (length(b) > 2) b[c(1, length(b))] else b
      },
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    # Exact y-axis control: bottom locked to y_min, top extends to highest y_break
    scale_y_continuous(
      limits = c(y_min, y_top_limit),
      breaks = y_breaks,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    facet_grid(~ var, scales = "free_x") +
    labs(
      title = pretty_lake_names[i],
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      # Title close to the rotated facet labels
      plot.title = element_text(size = 10, face = "plain", hjust = 0, margin = margin(b = 2)),
      # Rotated 45-degree variable names at top of each column
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.clip = "off",
      strip.text.x = element_text(
        angle = 45,
        hjust = 0,
        vjust = 0,
        size = 7.5,
        margin = margin(b = 2, l = 2)
      ),
      # Rotated 45-degree x-axis numbers & short ticks
      axis.text.x = element_text(angle=90, hjust=1, vjust=1, size=6.5, color="black"),
      axis.text.y = element_text(size=7.5, color="black"),
      axis.ticks.length = unit(2, "pt"),
      axis.line = element_line(linewidth=0.3),
      axis.ticks = element_line(linewidth=0.3),
      # Spacing between stratigraphic columns within the lake
      panel.spacing.x = unit(5, "pt"),
      # Outer margin around each lake panel (now actually responds!)
      plot.margin = margin(t=6, r=8, b=6, l=6, unit="pt")
    ) +
    coord_cartesian(clip="off")
  
  strat_plot_list[[i]] <- p
}


## create facet plot
# Reorder list so panels appear (a) through (h) left-to-right, top-to-bottom
panel_order <- order(pretty_lake_names)
ordered_plots <- strat_plot_list[panel_order]

# Use 'scale' (0.90 to 0.98) to adjust whitespace between all 8 panels
panels_grid <- cowplot::plot_grid(
  plotlist = ordered_plots,
  ncol=2,
  nrow=4
)

# Save to PDF (base dimensions in inches)
cowplot::save_plot(
  filename = paste0("Bio proxy CONISS/Bio proxy CONISS - v4/All_Lakes_strat_coniss_", date, ".pdf"),
  plot = panels_grid,
  base_width = 6.5,
  base_height = 8.5
)

### export zone boundary years
zone_boundary_years_df <- dplyr::bind_rows(zone_boundary_years)
write.csv(zone_boundary_years_df,file=paste0("Bio proxy CONISS/zone_boundary_years_CONISS_",date,".csv"))
