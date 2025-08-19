setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/Paleo PCoA")

librarian::shelf(tidyverse)

pcoa_loadings <- read.csv("WLPCoAScores.csv")
glimpse(pcoa_loadings)

