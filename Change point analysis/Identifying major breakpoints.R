librarian::shelf(tidyverse,grid,vegan,cowplot,zoo)

library(readr)
WL_paleo_master <- read_csv("GitHub/wildernesslakes_paleo_pfractions/Master dataset/WL_paleo_masterdataset_20Nov2024.csv")
glimpse(WL_paleo_master)

#pull Relevant data 
cluster <- WL_paleo_master %>% select(c(lake,year_loess,percent_organic,
                               tp_results_mg_p_g,
                               TOC_TN_ratio,
                               chl_a,
                               phaeo_b,
                               lutein,#indicates increase in green algae and tresstrial plants
                               diato, #increase in diatoms
                               echine,cantha,myxo, #cyanos
                               sio2_wt_percent,
                               ast_formosa, aul_ambigua,
                               aul_granulata,
                               aul_subarctica,
                               fra_crotonensis,
                               lin_comta,
                               nav_cryptocephala,
                               nit_perminuta,
                               pse_elliptica,
                               sel_nigri,
                               sel_saugerresii,
                               sta_construens,
                               sta_pinnata,
                               sta_venter,
                               tab_fenestrata,
                               tab_flocculosa))


#assign mixing_regime and trophic status based off of established contemporary data
cluster <- cluster %>%
  mutate(mixing_regime = case_when(
    lake == "flame" ~ "dimictic",
    lake == "wtwin" ~ "dimictic",
    lake == "finger" ~ "well_mixed",
    lake == "dunnigan" ~ "well_mixed",
    TRUE ~ "polymictic"),
    trophic_status=case_when(
    lake == "dunnigan" ~ "oligotrophic",
    lake == "wtwin" ~ "oligotrophic",
    TRUE ~ "mesotrophic"))  # Default case if none match

#####Just diatoms! do calculations
diatoms <- cluster %>% select(c(lake,year_loess,
                                        ast_formosa, aul_ambigua,
                                        aul_granulata,
                                        aul_subarctica,
                                        fra_crotonensis,
                                        lin_comta,
                                        nav_cryptocephala,
                                        nit_perminuta,
                                        pse_elliptica,
                                        sel_nigri,
                                        sel_saugerresii,
                                        sta_construens,
                                        sta_pinnata,
                                        sta_venter,
                                        tab_fenestrata,
                                        tab_flocculosa,
                                mixing_regime,trophic_status))
                             
                              
diatoms <- diatoms %>% mutate(
  strat_ind=ast_formosa+fra_crotonensis+lin_comta+tab_fenestrata+tab_flocculosa,
  dimictic_ind=aul_ambigua+aul_granulata+aul_subarctica,
  polymictic_ind=pse_elliptica+sta_construens+sta_pinnata+sta_venter,
  ratio_strat_dimictic=strat_ind/dimictic_ind,
  ratio_strat_polymictic=strat_ind/polymictic_ind)

#For some reason dunnigan is missing data? Remove and filter out NA
diatoms_clean <- diatoms %>%filter(lake!= "dunnigan")%>%
  filter(!is.na(ratio_strat_dimictic))


#Strat / Dimictic
diatoms %>% na.omit() %>% 
ggplot(aes(x=year_loess, y=ratio_strat_dimictic))+
  geom_line()+
  facet_wrap(~lake)

#Strat / Polymictic
diatoms %>% na.omit() %>% 
  ggplot(aes(x=year_loess, y=ratio_strat_polymictic))+
  geom_line()+
  facet_wrap(~lake, scales = "free")



library(strucchange) ###Don't know enough about this but there is lots of potential
#########Experimenting with strucchange#########
#Break point  analysis 
lm1 <- lm(ratio_strat_dimictic ~ year_loess+lake, data = diatoms_clean) #Create a linear model
coef1<- coef(lm1)
summary(lm1)


re.test <- efp(ratio_strat_dimictic ~ year_loess, data = diatoms_clean, type = "RE")
plot(re.test)

bp.lake <- breakpoints(ratio_strat_dimictic ~ year_loess, data = diatoms_clean, h = 0.1)
summary(bp.lake)
lines(bp.lake, breaks = 2)

plot(bp.lake)
breakpoints(bp.lake) #choose breakpoints around 2
######

#"segmented": Best for situations where you are specifically looking to identify 
#points where the relationship between variables changes from one linear trend to another

library(segmented)


###############Trying a loop, ratio_strat_dimictic#############
# Initialize an empty list to store results
brakepoint_results_list <- list()

# Loop through each lake
for (lake_name in unique(diatoms_clean$lake)) {
  
  # Subset the data for the current lake
  lake_data <- subset(diatoms_clean, lake == lake_name)
  
  # Fit the model
  lake_lm <- lm(ratio_strat_dimictic ~ year_loess, data = lake_data)
  segmented_model <- segmented(lake_lm)
  
  # Extract statistics
  breakpoint <- segmented_model$psi[, "Est."]
  ci <- confint(segmented_model)
  davies_pvalue <- davies.test(lake_lm)$p.value
  
  # Store results in a data frame
  lake_results <- data.frame(
    lake = lake_name,
    Breakpoint = round(breakpoint, 2),
    CI_Lower = round(ci[1, 1], 2),
    CI_Upper = round(ci[1, 2], 2),
    Davies_Test_pvalue = round(davies_pvalue, 4)
  )
  
  # Append to the results list
  brakepoint_results_list[[lake_name]] <- lake_results
}
# Combine all lake results into a single table
strat_dimictic_results_table <- do.call(rbind, brakepoint_results_list)

#Add in Trophic status and mixing regime
strat_dimictic_results_table <- strat_dimictic_results_table %>%
  mutate(contemporay_mixing_regime = case_when(
    lake == "flame" ~ "dimictic",
    lake == "wtwin" ~ "dimictic",
    lake == "finger" ~ "well_mixed",
    lake == "dunnigan" ~ "well_mixed",
    TRUE ~ "polymictic"),
    trophic_status=case_when(
      lake == "dunnigan" ~ "oligotrophic",
      lake == "wtwin" ~ "oligotrophic",
      TRUE ~ "mesotrophic"))  # Default case if none match


# Print the final table
print(strat_dimictic_results_table)

#Make a graph of the data
strat_dimictic_results_table %>%
  ggplot(aes(x = lake, y = Breakpoint, shape = contemporay_mixing_regime, color = -log10(Davies_Test_pvalue))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  theme_minimal() +
  scale_color_gradient(low = "pink", high = "red") + 
  labs(
    title = "Stratified:dimictic diatom indicator species Breakpoints",
    x = "Lake",
    y = "Breakpoint Year",
    shape = "contemporay_mixing_regime",
    color="-log10 (Davies_Test_pvalue)"
  ) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#############Trying a loop, ratio_strat_polymictic #############
brakepoint_results_list2 <- list()
# Loop through each lake
for (lake_name in unique(diatoms_clean$lake)) {
  
  # Subset the data for the current lake
  lake_data <- subset(diatoms_clean, lake == lake_name)
  
  # Fit the model
  lake_lm <- lm(ratio_strat_polymictic ~ year_loess, data = lake_data)
  segmented_model <- segmented(lake_lm)
  
  # Extract statistics
  breakpoint <- segmented_model$psi[, "Est."]
  ci <- confint(segmented_model)
  davies_pvalue <- davies.test(lake_lm)$p.value
  
  # Store results in a data frame
  lake_results <- data.frame(
    lake = lake_name,
    Breakpoint = round(breakpoint, 2),
    CI_Lower = round(ci[1, 1], 2),
    CI_Upper = round(ci[1, 2], 2),
    Davies_Test_pvalue = round(davies_pvalue, 4)
  )
  
  # Append to the results list
  brakepoint_results_list2[[lake_name]] <- lake_results
}
# Combine all lake results into a single table
strat_polymictic_results_table <- do.call(rbind, brakepoint_results_list2)


#Add in Trophic status and mixing regime
strat_polymictic_results_table <- strat_polymictic_results_table %>%
  mutate(contemporay_mixing_regime = case_when(
    lake == "flame" ~ "dimictic",
    lake == "wtwin" ~ "dimictic",
    lake == "finger" ~ "well_mixed",
    lake == "dunnigan" ~ "well_mixed",
    TRUE ~ "polymictic"),
    trophic_status=case_when(
      lake == "dunnigan" ~ "oligotrophic",
      lake == "wtwin" ~ "oligotrophic",
      TRUE ~ "mesotrophic"))  # Default case if none match


# Print the final table
print(strat_polymictic_results_table)

#Make a graph of the data
strat_polymictic_results_table %>%
ggplot(aes(x = lake, y = Breakpoint, shape = contemporay_mixing_regime, color = -log10(Davies_Test_pvalue))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  theme_minimal() +
  labs(
    title = "stratified:polymictic diatom indicator species Breakpoints",
    x = "Lake",
    y = "Breakpoint Year",
    shape = "contemporay_mixing_regime",
    color="-log10 (Davies_Test_pvalue)"
  ) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################################################














diatoms %>% na.omit() %>%  
  ggplot(aes(x=year_loess, y=ratio_strat_polymictic))+
  geom_line()+
  facet_wrap(~lake)


diatoms <- diatoms %>%
  mutate(mixing_indicator = case_when(
    species == "ast_formosa" ~ "stratification",
    species == "aul_ambigua" ~ "dimictic",
    species == "aul_granulata" ~ "dimictic",
    species == "aul_subarctica" ~ "dimictic",
    species == "fra_crotonensis" ~ "stratification",
    species == "lin_comta" ~ "stratification",
    species == "tab_fenestrata" ~ "stratification", #listed as can be strat indication
    species == "tab_flocculosa" ~ "stratification", #listed as can be strat indication
    species == "pse_elliptica" ~ "polymictic",
    species == "sta_construens" ~ "polymictic",
    species == "sta_pinnata" ~ "polymictic",
    species == "sta_venter" ~ "polymictic",
    TRUE ~ NA))
   



diatoms <- diatoms %>% pivot_longer(cols = c(3:18),
                                    names_to = "species",
                                    values_to = "count" )

#Calculate ratios















