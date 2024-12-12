####Making fancy graphs

librarian::shelf(tidyverse, googledrive,readxl,rioja,cowplot, patchwork)

master_dat <- read_csv("~/GitHub/wildernesslakes_paleo_pfractions/Master dataset/WL_paleo_masterdataset_20Nov2024.csv")

#add mixing regime
master_dat <- master_dat %>%
  mutate(mixing_regime = case_when(
    lake == "flame" ~ "dimictic",
    lake == "wtwin" ~ "dimictic",
    lake == "dunnigan" ~ "well_mixed",
    lake == "finger" ~ "well_mixed",
    TRUE ~ "polymictic"))  # Default case if none match

#convert things to flux
master_dat <- master_dat %>% mutate(recalcitrant_o_p_flux=recalcitrant_o_p*dmar_loess,
                                    labile_o_p_flux=labile_o_p*dmar_loess,
                                    tp_results_flux=tp_results_mg_p_g*dmar_loess)


#increasing productivity graph
###########################################################################
#increasing productivity graph
 productivity<- master_dat %>% select(c(lake,year_loess,percent_organic,
                                 recalcitrant_o_p_flux,labile_o_p_flux,
                                 ,chl_a, mixing_regime))
 
 #This looks not great
 productivity %>%
   pivot_longer(cols = c(percent_organic, recalcitrant_o_p_flux, labile_o_p_flux, chl_a),
                names_to = "variable", values_to = "value")%>%
   group_by(variable) %>%
   mutate(scaled_value = scale(value))%>%
   ggplot(aes(x = year_loess, y = scaled_value, color = variable)) +
   geom_line() +
   facet_wrap(~ lake, drop = TRUE)   # Rows for mixing regime, columns for lake

 
 
 color_palette <- c(
   # Dark Gray
   "#4B4B4B",  # Dark Gray
   
   # Blue Shades
   "#AEC6FF",  # Light Blue
   "#5A9BD4",  # Medium Blue
   "#2E6E1D",  # Dark Green
   
   
   "#FFE4B5",  # Very Light Orange (Lightened) tycho-planktonic taxa
   "#FFD700",  # Light Orange (Golden) tycho-planktonic taxa
   "#FFA500",  # Medium Orange (Original) tycho-planktonic taxa
   "#FF7F32",  # Darker Orange (Brighter)
   "#FF6600",  # Dark Orange (Pure Orange)
   "#FF4500",  # Darkest Orange (Vibrant Red-Orange)
   
   #Light gray
   "#d3d3d3", #silica
   
   # Green Shades
   "#D3E4B3",  # Very Light Green (new)
   "#B0E57C",  # Light Green
   "#66B032",  # Medium Green
   "#2E6E1D",  # Dark Green
   
   # Purple Shades
   "#D7BDE2",  # Light Purple
   "#A569BD",  # Medium Purple
   "#6C3483",  # Dark Purple
   
   # Brown Shades
   "#D2B48C",  # Light Brown (Tan)
   "#A0522D",  # Medium Brown
   "#4B2E00"   # Dark Brown
 )
 
 # Display the color palette
 color_palette
 
 
 #Burnt Lake plot
 burnt_productivity <- productivity %>%
   filter(lake=="burnt")
 burnt.productivity <- par(oma = c(0, 0, 2, 0)) #Leave space for the title
 burnt.productivity <-strat.plot(burnt_productivity[1:27,3:6],yvar=burnt_productivity$year_loess[1:27], #burnt[rows,col]
                           y.tks=seq(1860,2020,20),
                           plot.poly=T,plot.bar=T,col.bar="black",
                           col.poly = color_palette,
                           srt.xlabel=45,title="Burnt")

 
 
 
 #Elbow Lake plot
 elbow_productivity <- productivity %>%
   filter(lake=="elbow")
 elbow.productivity <- par(oma = c(0, 0, 2, 0)) #Leave space for the title, can't pipe
 elbow.productivity <-strat.plot(elbow_productivity[1:27,3:6],yvar=elbow_productivity$year_loess[1:27], #elbow[rows,col]
                           y.tks=seq(1860,2020,20),
                           plot.poly=T,plot.bar=T,col.bar="black",
                           col.poly = color_palette,
                           srt.xlabel=45,title="Elbow")
 
 
 #finger Lake plot
 finger_productivity <- productivity %>%
   filter(lake=="finger")
 finger.productivity <- par(oma = c(0, 0, 2, 0))
 finger.productivity <-strat.plot(finger_productivity[1:27,3:6],yvar=finger_productivity$year_loess[1:27], #finger[rows,col]
                            y.tks=seq(1860,2020,20),
                            plot.poly=T,plot.bar=T,col.bar="black",
                            col.poly = color_palette,
                            srt.xlabel=45,title="Finger")
 
 #Flame lake plot
 flame_productivity <- productivity %>%
   filter(lake=="flame")
 flame.productivity <- par(oma = c(0, 0, 2, 0))
 flame.productivity <-strat.plot(flame_productivity[1:27,3:6],yvar=flame_productivity$year_loess[1:27], #flame[rows,col]
                           y.tks=seq(1860,2020,20),
                           plot.poly=T,plot.bar=T,col.bar="black",
                           col.poly = color_palette,
                           srt.xlabel=45,title="Flame")
 
 
 
 #wtwin lake plot
 wtwin_productivity <- productivity %>%
   filter(lake=="wtwin")
 wtwin.productivity <- par(oma = c(0, 0, 2, 0))
 wtwin.productivity <-strat.plot(wtwin_productivity[1:27,3:6],yvar=wtwin_productivity$year_loess[1:27], #wtwin[rows,col]
                           y.tks=seq(1860,2020,20),
                           plot.poly=T,plot.bar=T,col.bar="black",
                           col.poly = color_palette,
                           srt.xlabel=45,title="West Twin")
 
 
 #smoke lake plot
 smoke_productivity <- productivity %>%
   filter(lake=="smoke")
 smoke.productivity <- par(oma = c(0, 0, 2, 0))
 smoke.productivity <- strat.plot(smoke_productivity[1:27,3:6],yvar=smoke_productivity$year_loess[1:27], #smoke[rows,col]
                           y.tks=seq(1860,2020,20),
                           plot.poly=T,plot.bar=T,col.bar="black",
                           col.poly = color_palette,
                           srt.xlabel=45,title="Smoke")
 
 #etwin lake plot
 etwin_productivity <- productivity %>%
   filter(lake=="etwin")
 etwin.productivity <- par(oma = c(0, 0, 2, 0))
 etwin.productivity <-strat.plot(etwin_productivity[1:27,3:6],yvar=etwin_productivity$year_loess[1:27], #etwin[rows,col]
                           y.tks=seq(1860,2020,20),
                           plot.poly=T,plot.bar=T,col.bar="black",
                           col.poly = color_palette,
                           srt.xlabel=45,title="East Twin")
 
 #dunnigan lake plot
 dunnigan_productivity <- productivity %>%
   filter(lake=="dunnigan")
 dunnigan.productivity <- par(oma = c(0, 0, 2, 0))
 dunnigan.productivity <-strat.plot(dunnigan_productivity[1:27,3:6],yvar=dunnigan_productivity$year_loess[1:27], #dunnigan[rows,col]
                              y.tks=seq(1860,2020,20),
                              plot.poly=T,plot.bar=T,col.bar="black",
                              col.poly = color_palette,
                              srt.xlabel=45,title="Dunnigan")
 
 
 library(gridExtra)
 
 # Assuming flame.productivity, wtwin.productivity, etc. are lists of stratigraphic plots
 # Each plot is stored as an individual plot object in the list
 
 # Example: Combine Dimictic Lakes plots (flame and wtwin)
 grid.arrange(
   flame.productivity[[2]],  # First stratigraphic plot from flame
   wtwin.productivity[[2]],  # First stratigraphic plot from wtwin
   ncol = 2,  # Arrange in 2 columns
   top = "Dimictic Lakes"  # Title for the top
 )
 #could not get this code to work
 
 
 ###########################################################################
 #Pfractions graph!
 #########################################################################
 #increasing productivity graph
 pfrac<- master_dat %>% select(c(lake,year_loess,
                                 ex_p,bd_p,na_oh_p,
                                recalcitrant_o_p,labile_o_p,tp_results_mg_p_g,
                                         mixing_regime))
#convert to %
 pfrac <- pfrac %>% mutate(loosely_bound_percent=(ex_p/tp_results_mg_p_g)*100,
                  iron_bound_percent=(bd_p/tp_results_mg_p_g)*100,
                  aluminum_bound_percent=(na_oh_p/tp_results_mg_p_g)*100,
                  labile_o_percent=(labile_o_p/tp_results_mg_p_g)*100,
                  recalcitrant_o_p=(recalcitrant_o_p/tp_results_mg_p_g)*100)


 #means of data
 summary_table <- pfrac %>%
 group_by(lake,mixing_regime) %>%
   summarize(
     iron_bound_percent_mean = mean(iron_bound_percent, na.rm = TRUE),
     aluminum_bound_percent_mean = mean(aluminum_bound_percent, na.rm = TRUE),
     .groups = "drop"  # This removes the grouping after summarization
   )
 
 #summary table
 library(knitr)
 
 # Use kable to create a more readable table
 kable(summary_table, caption = "Fe-P and Al-P Means")
 
 library(writexl)
 
 write_xlsx(summary_table, "C:/Users/16512/Documents/GitHub/wildernesslakes_paleo_pfractions/Graphs/mean2_summary_table.xlsx")
 
 
 
 pfrac %>% 
   pivot_longer(cols = c(loosely_bound_percent, iron_bound_percent, aluminum_bound_percent, labile_o_percent,recalcitrant_o_p),
                names_to = "pfraction", values_to = "percentage")%>%
   na.omit() %>% 
   ggplot(aes(x=percentage, y=year_loess, color=pfraction))+
   geom_point()+
   geom_path()+
   facet_wrap(~lake, scales="free")
   
    
#Nitrogen Graphs
###########################################################################
 
 nitrogen<- master_dat %>% select(c(lake,year_loess,TOC_TN_ratio,
                                    tab_flocculosa, ast_formosa, fra_crotonensis,
                                    mixing_regime))
 
 long_data <- nitrogen %>%
   pivot_longer(cols = c(tab_flocculosa, ast_formosa, fra_crotonensis),
                names_to = "species", values_to = "abundance")
 
 # Plotting the data diatoms
 long_data %>% 
 ggplot(aes(x = year_loess, y = abundance, color = species, group = species)) +
   geom_point() +
   facet_wrap(~lake, scales = "free") +  # Facet by mixing regime
   labs(title = "Species Abundance Over Time by Mixing Regime",
        x = "Year", y = "Abundance") +
   theme(legend.position = "bottom")
 
 #Plotting the data C:N
 nitrogen %>% 
 ggplot(aes(x = year_loess, y =TOC_TN_ratio , color = lake)) +
   geom_point() +
   geom_path()+
   facet_wrap(~mixing_regime) +  # Facet by mixing regime
   labs(title = "C:N Ratio",
        x = "Year", y = "C:N") +
   theme(legend.position = "bottom")
 
 
 
 
 
 
 