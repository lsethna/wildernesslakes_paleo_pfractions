## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##
rm(list=ls())

#install.packages(librarian)
librarian::shelf(tidyverse, googledrive,readxl,rioja)

getwd()
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions") #change this to match local GitHub folder

## ----------------------------------- ##
# Read in data ----
## ----------------------------------- ##

#interpolated master dat
int_master_dat <- read.csv("raw_data/interpolated_master_dat_10Dec24.csv")
glimpse(int_master_dat)

## ----------------------------------- ##
# Select data for strat plots ----
## ----------------------------------- ##

#select variables for algae strat diagram
paleo_prod <- int_master_dat %>% select(c(lake,year_loess,
                                          chl_b,fuco,b_car,lutein,diato,allo,myxo,cantha,echine,
                                          aul_subarctica,aul_ambigua,pse_elliptica,ast_formosa,tab_flocculosa,fra_crotonensis,sta_venter))
glimpse(paleo_prod)
unique(paleo_prod$lake)

## ------------------------------- ##
## --- strat plots with ggplot --- ##
## ------------------------------- ##

paleo_prod %>%
  mutate(lake = factor(lake,levels=c("dunnigan","finger","flame","wtwin","burnt","smoke","elbow","etwin"),
                       labels=c("Dunnigan","Finger","Flame","W. Twin","Burnt","Smoke","Elbow","E. Twin"))) %>%
  pivot_longer(cols=chl_b:sta_venter) %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area()+
  coord_flip()+
  facet_grid(rows=vars(lake),cols=vars(name))+
  theme_classic()

#Need to plot lakes separately to control scales

## ------------------------------ ##
## --- strat plots with rioja --- ##
## ------------------------------ ##

paleo_prod_long <- 
  paleo_prod %>%
  pivot_longer(cols=chl_b:sta_venter) %>%
  mutate(name = factor(name,
                       levels=c("chl_b","fuco","b_car","lut_zea","diato","allo","myxo","cantha","echine",
                                "aul_subarctica","aul_ambigua","pse_elliptica","ast_formosa","tab_flocculosa",
                                "fra_crotonensis","sta_venter")))
#plot theme for all 
plot_theme <- 
  theme_classic()+
  theme(strip.background=element_blank(),
        strip.text=element_text(angle=45),
        axis.text.x=element_text(angle=90),
        panel.spacing=unit(1, "lines"))

#Dunnigan
dunnigan_paleo <- paleo_prod_long %>%
  filter(lake=="dunnigan") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="Dunnigan (well-mixed)")+
  plot_theme

#Finger
finger_paleo <- paleo_prod_long %>%
  filter(lake=="finger") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="Finger (well-mixed)")+
  plot_theme

#Flame
flame_paleo <- paleo_prod_long %>%
  filter(lake=="flame") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="Flame (dimictic)")+
  plot_theme

#West Twin
wtwin_paleo <- paleo_prod_long %>%
  filter(lake=="wtwin") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="West Twin (dimictic)")+
  plot_theme

#Burnt
burnt_paleo <- paleo_prod_long %>%
  filter(lake=="burnt") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="Burnt (polymictic)")+
  plot_theme
  
#Smoke
smoke_paleo <- paleo_prod_long %>%
  filter(lake=="smoke") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="Smoke (polymictic)")+
  plot_theme

#Elbow
elbow_paleo <- paleo_prod_long %>%
  filter(lake=="elbow") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="Elbow (polymictic)")+
  plot_theme

#East Twin
etwin_paleo <- paleo_prod_long %>%
  filter(lake=="etwin") %>%
  ggplot(aes(y=value,x=year_loess))+
  geom_area(color="black",fill="gray")+
  geom_bar(stat="identity",width=0.5)+
  coord_flip()+
  scale_x_continuous(breaks=seq(1860,2020,20))+
  facet_wrap(~name,ncol=16,scales="free_x")+
  labs(x="Year",title="East Twin (polymictic)")+
  plot_theme

## ------------------------ ##
## --- cowplot together --- ##
## ------------------------ ##

cowplot::plot_grid(dunnigan_paleo,finger_paleo,
                   flame_paleo,wtwin_paleo,
                   burnt_paleo,smoke_paleo,
                   elbow_paleo,etwin_paleo,
                   nrow=4,labels="AUTO")

