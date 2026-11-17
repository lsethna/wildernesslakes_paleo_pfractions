## create stacked WL drivers and responses for Wilson-Jackson paleo MS

rm(list=ls())
librarian::shelf(tidyverse)

getwd()
setwd("/Users/lsethna/Documents/GitHub/wildernesslakes_paleo_pfractions")

#read in data
#PCoA scores
pcoa <- read.csv("Paleo PCoA/PCoA_scores_by_lake_24Jul26.csv")
#climate
climate <- read.csv("raw_data/tower_annual_temperatures.csv")
#deriv org C
carbon <- read.csv("raw_data/deriv_orgC.csv")
#p fractions
pfrac <- read.csv("raw_data/pfrac.csv")

## climate plot
glimpse(climate)
climate_plot <- 
ggplot(climate,aes(x=Year,y=(Annual_Average_Temperature_F-32)*(5/9))) +
  geom_line(color="gray20") +
  scale_x_continuous(limits=c(1825,2025),breaks=seq(1825,2025,25)) +
  scale_y_continuous(breaks=seq(0,6,1),name="Avg annual \ntemp. (°C)") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none")

#pcoa scores plot
glimpse(pcoa)
pcoa_plot <- 
pcoa %>%
mutate(lake = factor(lake,levels=c("dunnigan","finger","burnt","smoke","elbow","etwin","flame","wtwin"),
                     labels=c("Dunnigan","Finger","Burnt","Smoke","Elbow","East Twin","Flame","West Twin"))
) %>%
  ggplot(aes(x=year,y=X1,color=lake))+
  #geom_vline(xintercept = 1900,lty="dashed")+
  #geom_vline(xintercept = 1970,lty="dashed")+
  geom_smooth(se=F,span=0.7)+ #span - smaller numbers produce wigglier lines, larger numbers produce smoother lines
  scale_color_manual(values=c("#980B83", #Dunnigan
                              "#f464df", #Finger
                              "#122135", #Burnt
                              "#3968A6", #Smoke
                              "#7199CE", #East Twin
                              "#B8CCE7", #Elbow
                              "#ffa600", #Flame
                              "#FFD586" #West Twin
  )) +
  labs(y="PCoA Axis 1",color="Lake") +
  scale_x_continuous(limits=c(1825,2025),breaks=seq(1825,2025,25)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none")

#carbon plot
glimpse(carbon)
carbon_plot <- 
carbon %>%
  mutate(lake = factor(lake,levels=c("dunnigan","finger","burnt","smoke","elbow","etwin","flame","wtwin"),
                       labels=c("Dunnigan","Finger","Burnt","Smoke","Elbow","East Twin","Flame","West Twin"))
  ) %>%
ggplot(aes(x=year_loess,y=deriv_orgC,color=lake)) +
  geom_line() +
  scale_color_manual(values=c("#980B83", #Dunnigan
                              "#f464df", #Finger
                              "#122135", #Burnt
                              "#3968A6", #Smoke
                              "#7199CE", #East Twin
                              "#B8CCE7", #Elbow
                              "#ffa600", #Flame
                              "#FFD586" #West Twin
  )) +
  labs(y="dC/dt",color="Lake") +
  scale_y_continuous(limits=c(-0.03,0.03))+
  scale_x_continuous(limits=c(1825,2025),breaks=seq(1825,2025,25)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none") 

#pfrac plot
glimpse(pfrac)
pfrac_plot <- 
pfrac %>%
  filter(name %in% c("labile","ironbound")) %>%
  rename(depth=bottom_interval_cm)%>%
  left_join(carbon) %>%
  mutate(lake = factor(lake,levels=c("dunnigan","finger","burnt","smoke","elbow","etwin","flame","wtwin"),
                       labels=c("Dunnigan","Finger","Burnt","Smoke","Elbow","East Twin","Flame","West Twin"))
  ) %>%
  ggplot(aes(x=year_loess,y=value,color=lake,lty=name)) +
  geom_line() +
  scale_color_manual(values=c("#980B83", #Dunnigan
                              "#f464df", #Finger
                              "#122135", #Burnt
                              "#3968A6", #Smoke
                              "#7199CE", #East Twin
                              "#B8CCE7", #Elbow
                              "#ffa600", #Flame
                              "#FFD586" #West Twin
  )) +
  labs(y="P fraction (mg P L \u207b\u00b9)") +
  scale_x_continuous(limits=c(1825,2025),breaks=seq(1825,2025,25)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) 

## cowplot
 cowplot::plot_grid(climate_plot,pcoa_plot,carbon_plot,pfrac_plot,
                    ncol=1,align="v",rel_heights = c(1,1,1,1.78))

 
 #get land use times
 land_use <- carbon %>%
   select(contains("landuse"),lake) %>%
   distinct()
