rm(list=ls())

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/Paleo PCoA")

librarian::shelf(tidyverse)

pcoa_loadings <- read.csv("WLPCoAScores.csv")
glimpse(pcoa_loadings)

ggplot(pcoa_loadings, aes(x=year,y=Dim1,color=lake))+
  geom_line(linewidth=1) +
  scale_color_manual(values=c("#d77932", #burnt
                              "#8b5a2b", #dunnigan
                              "#9a77b5", #elbow
                              "#6d91c1", #etwin
                              "#62497a", #finger
                              "#9b3d3d", #flame
                              "#e3b448", #smoke
                              "#3b5b7c" #wtwin
  )) +
  theme_classic()
