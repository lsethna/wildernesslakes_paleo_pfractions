librarian::shelf(tidyverse,grid,vegan,cowplot,zoo)

rm(list=ls())

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/Paleo PCoA")

#read in master data
master_dat <- read.csv("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/raw_data/WL_paleo_masterdataset_20Nov2024.csv") %>%
  select(!X)
glimpse(master_dat)

## ------------------------------------------------------------ ##
## ---- filter variables: pigments>0, diatoms>2% abundance ---- ##
## ------------------------------------------------------------ ##
vars <- data.frame(
  column_n=1:length(master_dat),
  vars=colnames(master_dat))
vars
#pigments 25:56
#diatoms 57:270
sel_diat <- master_dat %>% select(colnames(master_dat[,c(56:269)])) %>% 
  pivot_longer(cols=ach_microcephala:uln_ulna) %>% filter(value>5) %>% #filter to 5% relative abundance (might want to only filter to 2% and at least in two samples)
  select(name) %>% distinct()
  
##
sel_p <- master_dat %>% select(colnames(master_dat[,24:55])) %>% 
  summarize(across(chlide_a:car_z,max,na.rm=T))
sel_pig <- sel_p %>% 
  select(!c("sudan","chl_c1","chl_c2","fuco","diadino","chl_b","chl_a","chl_ap")) %>% #remove unstable or non-interpretable pigments
  pivot_longer(cols=chlide_a:car_z) %>% filter(value>0) #remove any non-detect pigments

#filter master data based on selected diatoms and pigments
master_v2 <- master_dat %>% select(lake,year_loess,depth, #id info
                                   9:11, #loi; perc inorg,org,calc
                                   #16:23, #P and P fractions
                                   14,sel_diat$name, #BSi flux, select diatoms
                                   sel_pig$name, #select pigments
                                   #274:276 #TOC, TN, TOC:TN

)  
glimpse(master_v2)

## ------------------------------- ##
## ---- interpolate variables ---- ##
## ------------------------------- ##
master_v3_interp <- master_v2 %>%
  group_by(lake) %>%
  mutate_at(vars(percent_organic:b_car),funs(zoo::na.approx(.,method="constant",rule=2))) %>% #rule=2 means extend nearest values to leading/trailing NAs
  ungroup() %>%
  drop_na() #remove any NA values - beyond dated intervals

#quick check to make sure all dates and dates are included
ggplot(master_v3_interp,aes(x=year_loess,y=depth))+
  geom_line()+facet_wrap(~lake)

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions")
write.csv(master_v3_interp,file="raw_data/interpolated_master_dat_4Sep25.csv")

## -------------- ##
## ---- PCoA ---- ##
## -------------- ##

colnames(master_v3_interp)
#diatoms=8:33
#pigments=34:45
#geochem=4:7

#standardize diatoms
diat.stand.hell <- vegan::decostand(master_v3_interp[,8:33], method="hellinger")
#diatom distance matrix
diat.dist <- dist(diat.stand.hell, method="euclidean")
#normalize the diatom vectors to allow for comparison across proxies
diat.dist <- diat.dist/max(diat.dist)

#standardize pigments
pig.stand.hell <- decostand(master_v3_interp[,34:45], method="hellinger")
#pigment distance matrix
pig.dist <- dist(pig.stand.hell, method="euclidean")
pig.dist <- pig.dist/max(pig.dist)

#standardize geochem
# geochem.stand.hell <- decostand(master_v3_interp[,c(4:7)], method="hellinger")
# #geochemment distance matrix
# geochem.dist <- dist(geochem.stand.hell, method="euclidean")
# geochem.dist <- geochem.dist/max(geochem.dist)

#combine distance matrices
comb.dist <- diat.dist+pig.dist#+geochem.dist
#and scale
comb.mds <- data.frame(cmdscale(d=comb.dist, k=2))

plot(comb.mds)

#get scores
site.comb <- data.frame(scores(comb.mds))

#add sample IDs
site.comb$lake <- master_v3_interp$lake
site.comb$year <- master_v3_interp$year_loess

#combined vectors for all vars
comb.vec <- envfit(comb.mds, env=master_v3_interp[,-3],na.rm=T) #include all paleo vars except depth
comb.vec <- data.frame(cbind(comb.vec$vectors$arrows, comb.vec$vectors$r, comb.vec$vectors$pvals))
comb.vec.sig <- subset(comb.vec, comb.vec$V4 <= 0.05 & comb.vec$V3 >= 0.6) #only saves vectors that are significant (p<0.05) and are >60% correlated with ordination
#just diatom vectors
diat.vec <- envfit(comb.mds,env=master_v3_interp[,8:33],na.rm=T)
diat.vec <- data.frame(cbind(diat.vec$vectors$arrows, diat.vec$vectors$r, diat.vec$vectors$pvals))
diat.vec.sig <- subset(diat.vec, diat.vec$V4 <= 0.05 & diat.vec$V3 >= 0.6)
#just pigment vectors
pig.vec <- envfit(comb.mds,env=master_v3_interp[,34:45],na.rm=T)
pig.vec <- data.frame(cbind(pig.vec$vectors$arrows, pig.vec$vectors$r, pig.vec$vectors$pvals))
pig.vec.sig <- subset(pig.vec, pig.vec$V4 <= 0.05 & pig.vec$V3 >= 0.6)
#just geochem vectors
# geochem.vec <- envfit(comb.mds,env=master_v3_interp[,c(4:7)],na.rm=T)
# geochem.vec <- data.frame(cbind(geochem.vec$vectors$arrows, geochem.vec$vectors$r, geochem.vec$vectors$pvals))
# geochem.vec.sig <- subset(geochem.vec, geochem.vec$V4 <= 0.05 & geochem.vec$V3 >= 0.6)

# Scale arrows by r2 to show correlation strength
pig.vec.sig.scaled <- pig.vec.sig %>%
  mutate(X1 = (X1*V3)/2,
         X2 = (X2*V3)/2)
diat.vec.sig.scaled <- diat.vec.sig %>%
  mutate(X1 = (X1*V3)/2,
         X2 = (X2*V3)/2)
# geochem.vec.sig.scaled <- geochem.vec.sig %>%
#   mutate(X1 = (X1*V3)/2,
#          X2 = (X2*V3)/2)

#plot PCoA paths for all lakes
ggplot(aes(X1, X2),data=site.comb) +  
  #lake trajectories 
  geom_path(aes(color=lake), size=0.75, 
            arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), 
            lineend="round") + 
  #pigment vectors
  geom_segment(data=pig.vec.sig.scaled,
               aes(x=0,y=0,xend=X1,yend=X2), 
               color="darkgreen",arrow=arrow(length=unit(0.15, "inches"))) + 
  ggrepel::geom_label_repel(data=pig.vec.sig.scaled,
            aes(x=X1*1.25, y=X2*1.25),label=rownames(pig.vec.sig.scaled), 
            color="darkgreen",label.size=NA,fill=NA) +  
  #diatom vectors
  geom_segment(data=diat.vec.sig.scaled,
               aes(x=0,y=0,xend=X1,yend=X2), 
               color="orange",arrow=arrow(length=unit(0.15, "inches"))) + 
  ggrepel::geom_label_repel(data=diat.vec.sig.scaled,
            aes(x=X1*1.25, y=X2*1.25),label=rownames(diat.vec.sig.scaled), 
            color="orange", label.size=NA,fill=NA) +
  #geochem vectors
  # geom_segment(data=geochem.vec.sig.scaled,
  #              aes(x=0,y=0,xend=X1,yend=X2), 
  #              color="hotpink",arrow=arrow(length=unit(0.15, "inches"))) + 
  # ggrepel::geom_label_repel(data=geochem.vec.sig.scaled,
  #           aes(x=X1*1.25, y=X2*1.25),label=rownames(geochem.vec.sig.scaled), 
  #           color="hotpink", label.size=NA,fill=NA) +
  #general aesthetics
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Lake") +
  scale_color_manual(values=c("#C1D9B7", #burnt
                              "#89CFF0", #dunnigan
                              "#007CAA", #elbow
                              "#665191", #etwin
                              "#d45087", #finger
                              "#7A871E", #flame
                              "#104210", #smoke
                              "#E97451" #wtwin
  )) +
  theme_classic(base_size=14) 

## ------------------------------------- ##
## ---- PCoA for each mixing regime ---- ##
## ------------------------------------- ##

#add mix regime to master data
master_v3_interp <- master_v3_interp %>%
  mutate(mix.regime = case_when(lake=="burnt"~"Polymictic",
                                lake=="etwin"~"Polymictic",
                                lake=="elbow"~"Polymictic",
                                lake=="smoke"~"Polymictic",
                                lake=="dunnigan"~"Well mixed",
                                lake=="finger"~"Well mixed",
                                lake=="wtwin"~"Dimictic",
                                lake=="flame"~"Dimictic"))

#set up variables to run PcoA for each mixing regime in a loop
mixing_regime = unique(master_v3_interp$mix.regime)
pcoa.scores_mixing.regime <- list() #list to save scores output
pcoa.variable.vectors <- list() #list to save variable vectors
pcoa.plots <- list()

for (i in 1:length(mixing_regime)) {
  #create df for each mixing regime
  mix.master.dat <- subset(master_v3_interp,master_v3_interp$mix.regime==mixing_regime[i])
  
  #standardize diatoms and pigments... 
  diat.stand.hell.mix <- decostand(mix.master.dat[,8:33], method="hellinger")
  pig.stand.hell.mix <- decostand(mix.master.dat[,34:45], method="hellinger")
  #...to create distance matrices
  diat.dist.mix <- dist(diat.stand.hell.mix, method="euclidean")
  pig.dist.mix <- dist(pig.stand.hell.mix, method="euclidean")
  
  #combine distance matrices...
  comb.dist.mix <- diat.dist.mix+pig.dist.mix
  #...and scale
  comb.mds.mix <- cmdscale(d=comb.dist.mix, k=2)
  
  #get scores
  comb.scores.mix <- data.frame(scores(comb.mds.mix))
  #add sample IDs
  comb.scores.mix$lake <- mix.master.dat$lake
  comb.scores.mix$year <- mix.master.dat$year_loess
  comb.scores.mix$mix_regime <- mixing_regime[i]
  
  #save to list
  pcoa.scores_mixing.regime[[i]] <- comb.scores.mix
  
  #calculate vectors for diatoms and pigments
  #diatom vectors
  diat.vec.mix <- envfit(comb.mds.mix,env=mix.master.dat[,8:33],na.rm=T)
  diat.vec.mix <- data.frame(cbind(diat.vec.mix$vectors$arrows, 
                                   diat.vec.mix$vectors$r, 
                                   diat.vec.mix$vectors$pvals),
                             var_type="diatoms")
  #pigments
  pig.vec.mix <- envfit(comb.mds.mix,env=mix.master.dat[,34:45],na.rm=T)
  pig.vec.mix <- data.frame(cbind(pig.vec.mix$vectors$arrows, 
                                  pig.vec.mix$vectors$r, 
                                  pig.vec.mix$vectors$pvals),
                                  var_type="pigment")
  #save these to list
  vec.mix <- rbind(diat.vec.mix,pig.vec.mix) %>% mutate(mix_regime=mixing_regime[1])
  pcoa.variable.vectors[[i]] <- vec.mix
  
  #keep only significant vectors for plotting
  diat.vec.mix.sig <- subset(diat.vec.mix, diat.vec.mix$V4 <= 0.05 & diat.vec.mix$V3 >= 0.6)
  pig.vec.mix.sig <- subset(pig.vec.mix, pig.vec.mix$V4 <= 0.05 & pig.vec.mix$V3 >= 0.6)
  
  # Scale arrows by r2 to show correlation strength
  diat.vec.mix.sig.scaled <- diat.vec.mix.sig %>%
    mutate(Dim1 = (Dim1*V3)/2,
           Dim2 = (Dim2*V3)/2)
  pig.vec.mix.sig.scaled <- pig.vec.mix.sig %>%
    mutate(Dim1 = (Dim1*V3)/2,
           Dim2 = (Dim2*V3)/2)
  
  #create PCoA biplot
  #colors
  comb.scores.mix <- comb.scores.mix %>% mutate(color = case_when(lake=="burnt" ~ "#C1D9B7",
                                                                  lake=="dunnigan" ~ "#89CFF0",
                                                                  lake=="elbow" ~ "#007CAA",
                                                                  lake=="etwin" ~ "#665191",
                                                                  lake=="finger" ~"#d45087",
                                                                  lake=="flame" ~"#7A871E",
                                                                  lake=="smoke" ~ "#104210",
                                                                  lake=="wtwin" ~"#E97451") )

  p <- 
  ggplot(aes(Dim1, Dim2),data=comb.scores.mix) +  
    #lake trajectories 
    geom_path(aes(color=lake), size=0.75, 
              arrow=arrow(type="closed", ends="first",length=unit(0.05,"inches")), lineend="round") + 
    #pigment vectors
    geom_segment(data=pig.vec.mix.sig.scaled,
                 aes(x=0,y=0,xend=Dim1,yend=Dim2), 
                 color="darkgreen",arrow=arrow(length=unit(0.15, "inches"))) + 
    geom_text(data=pig.vec.mix.sig.scaled,
              aes(x=Dim1*1.2, y=Dim2*1.2),label=rownames(pig.vec.mix.sig.scaled), 
              color="darkgreen", size=4,vjust=-0.5,
              position=position_jitter(width=0.1,height=0.1)) +  
    #diatom vectors
    geom_segment(data=diat.vec.mix.sig.scaled,
                 aes(x=0,y=0,xend=Dim1,yend=Dim2), 
                 color="orange",arrow=arrow(length=unit(0.15, "inches"))) + 
    geom_text(data=diat.vec.mix.sig.scaled,
              aes(x=Dim1*1.2, y=Dim2*1.2),label=rownames(diat.vec.mix.sig.scaled), 
              color="orange", size=4,vjust=-0.5,
              position=position_jitter(width=0.1,height=0.1)) +
    #general aesthetics
    labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Lake") +
    ggtitle(mixing_regime[i]) +
    scale_color_manual(values=unique(comb.scores.mix$color)) +
    theme_bw(base_size=14) +
    theme(legend.position="bottom")

  pcoa.plots[[i]] <- p
}

## ----------------------------------------------------------- ##
## ---- PCoA for one lake representing each mixing regime ---- ##
## ----------------------------------------------------------- ##

#create new dfs for one lake in each mixing regime
burnt.master.dat <- subset(poly.master.dat,poly.master.dat$lake=="burnt")
finger.master.dat <- subset(mixed.master.dat,mixed.master.dat$lake=="finger")
wtwin.master.dat <- subset(di.master.dat,di.master.dat$lake=="wtwin")

#check
glimpse(wtwin.master.dat)
colnames(finger.master.dat)

#standardize geochemistry
geochem.stand.hell.burnt <- decostand(burnt.master.dat[,c(4,8,13)],method="hellinger")
geochem.stand.hell.finger <- decostand(finger.master.dat[,c(4,8,13)],method="hellinger")
geochem.stand.hell.wtwin <- decostand(wtwin.master.dat[,c(4,8,13)],method="hellinger")
#geochemistry distance matrix
geochem.dist.burnt <- dist(geochem.stand.hell.burnt, method="euclidean")
geochem.dist.finger <- dist(geochem.stand.hell.finger, method="euclidean")
geochem.dist.wtwin <- dist(geochem.stand.hell.wtwin, method="euclidean")

geochem.dist.burnt <- geochem.dist.burnt/max(geochem.dist.burnt)
geochem.dist.finger <- geochem.dist.finger/max(geochem.dist.finger)
geochem.dist.wtwin <- geochem.dist.wtwin/max(geochem.dist.wtwin)

#standardize diatoms
diat.stand.hell.burnt <- decostand(burnt.master.dat[,16:31], method="hellinger")
diat.stand.hell.finger <- decostand(finger.master.dat[,16:31], method="hellinger")
diat.stand.hell.wtwin <- decostand(wtwin.master.dat[,16:31], method="hellinger")
#diatom distance matrix
diat.dist.burnt <- dist(diat.stand.hell.burnt, method="euclidean")
diat.dist.finger <- dist(diat.stand.hell.finger, method="euclidean")
diat.dist.wtwin <- dist(diat.stand.hell.wtwin, method="euclidean")

diat.dist.burnt <- diat.dist.burnt/max(diat.dist.burnt)
diat.dist.finger <- diat.dist.finger/max(diat.dist.finger)
diat.dist.wtwin <- diat.dist.wtwin/max(diat.dist.wtwin)

#standardize pigments
pig.stand.hell.burnt <- decostand(burnt.master.dat[,32:51], method="hellinger")
pig.stand.hell.finger <- decostand(finger.master.dat[,32:51], method="hellinger")
pig.stand.hell.wtwin <- decostand(wtwin.master.dat[,32:51], method="hellinger")
#pigment distance matrix
pig.dist.burnt <- dist(pig.stand.hell.burnt, method="euclidean")
pig.dist.finger <- dist(pig.stand.hell.finger, method="euclidean")
pig.dist.wtwin <- dist(pig.stand.hell.wtwin, method="euclidean")

pig.dist.burnt <- pig.dist.burnt/max(pig.dist.burnt)
pig.dist.finger <- pig.dist.finger/max(pig.dist.finger)
pig.dist.wtwin <- pig.dist.wtwin/max(pig.dist.wtwin)

#combine distance matrices
burnt.comb.dist <- geochem.dist.burnt+diat.dist.burnt+pig.dist.burnt
finger.comb.dist <- geochem.dist.finger+diat.dist.finger+pig.dist.finger
wtwin.comb.dist <- geochem.dist.wtwin+diat.dist.wtwin+pig.dist.wtwin

#and scale
burnt.comb.mds <- cmdscale(d=burnt.comb.dist, k=2)
finger.comb.mds <- cmdscale(d=finger.comb.dist, k=2)
wtwin.comb.mds <- cmdscale(d=wtwin.comb.dist, k=2)
#plot
plot(burnt.comb.mds)
plot(finger.comb.mds)
plot(wtwin.comb.mds)
#get scores
burnt.comb.scores <- data.frame(scores(burnt.comb.mds))
finger.comb.scores <- data.frame(scores(finger.comb.mds))
wtwin.comb.scores <- data.frame(scores(wtwin.comb.mds))

#add site IDs
burnt.comb.scores$lake <- burnt.master.dat$lake
burnt.comb.scores$year <- burnt.master.dat$year_loess
finger.comb.scores$lake <- finger.master.dat$lake
finger.comb.scores$year <- finger.master.dat$year_loess
wtwin.comb.scores$lake <- wtwin.master.dat$lake
wtwin.comb.scores$year <- wtwin.master.dat$year_loess

## geochem, diatom, and pigment vectors for each mixing regime ##
#burnt geochem vectors
burnt.geo.vec <- envfit(burnt.comb.mds,env=burnt.master.dat[,c(4,8,13)],na.rm=T)
burnt.geo.vec <- data.frame(cbind(burnt.geo.vec$vectors$arrows, 
                                  burnt.geo.vec$vectors$r, 
                                  burnt.geo.vec$vectors$pvals))
burnt.geo.vec.sig <- subset(burnt.geo.vec, burnt.geo.vec$V4 <= 0.001 & burnt.geo.vec$V3 >= 0.5)
#finger geochem vectors
finger.geo.vec <- envfit(finger.comb.mds,env=finger.master.dat[,c(4,8,13)],na.rm=T)
finger.geo.vec <- data.frame(cbind(finger.geo.vec$vectors$arrows, 
                                   finger.geo.vec$vectors$r, 
                                   finger.geo.vec$vectors$pvals))
finger.geo.vec.sig <- subset(finger.geo.vec, finger.geo.vec$V4 <= 0.001 & finger.geo.vec$V3 >= 0.5)
#wtwin geochem vectors
wtwin.geo.vec <- envfit(wtwin.comb.mds,env=wtwin.master.dat[,c(4,8,13)],na.rm=T)
wtwin.geo.vec <- data.frame(cbind(wtwin.geo.vec$vectors$arrows, 
                                  wtwin.geo.vec$vectors$r, 
                                  wtwin.geo.vec$vectors$pvals))
wtwin.geo.vec.sig <- subset(wtwin.geo.vec, wtwin.geo.vec$V4 <= 0.001 & wtwin.geo.vec$V3 >= 0.5)

#Burnt diatom vectors
burnt.diat.vec <- envfit(burnt.comb.mds,env=burnt.master.dat[,16:31],na.rm=T)
burnt.diat.vec <- data.frame(cbind(burnt.diat.vec$vectors$arrows, 
                                   burnt.diat.vec$vectors$r, 
                                   burnt.diat.vec$vectors$pvals))
burnt.diat.vec.sig <- subset(burnt.diat.vec, burnt.diat.vec$V4 <= 0.001 & burnt.diat.vec$V3 >= 0.5)
#finger mixed diatom vectors
finger.diat.vec <- envfit(finger.comb.mds,env=finger.master.dat[,16:31],na.rm=T)
finger.diat.vec <- data.frame(cbind(finger.diat.vec$vectors$arrows, 
                                   finger.diat.vec$vectors$r, 
                                   finger.diat.vec$vectors$pvals))
finger.diat.vec.sig <- subset(finger.diat.vec, finger.diat.vec$V4 <= 0.001 & finger.diat.vec$V3 >= 0.5)
#wtwin diatom vectors
wtwin.diat.vec <- envfit(wtwin.comb.mds,env=wtwin.master.dat[,16:31],na.rm=T)
wtwin.diat.vec <- data.frame(cbind(wtwin.diat.vec$vectors$arrows, 
                                wtwin.diat.vec$vectors$r, 
                                wtwin.diat.vec$vectors$pvals))
wtwin.diat.vec.sig <- subset(wtwin.diat.vec, wtwin.diat.vec$V4 <= 0.001 & wtwin.diat.vec$V3 >= 0.5)
#burnt pigment vectors
burnt.pig.vec <- envfit(burnt.comb.mds,env=burnt.master.dat[,32:51],na.rm=T)
burnt.pig.vec <- data.frame(cbind(burnt.pig.vec$vectors$arrows, 
                                 burnt.pig.vec$vectors$r, 
                                 burnt.pig.vec$vectors$pvals))
burnt.pig.vec.sig <- subset(burnt.pig.vec, burnt.pig.vec$V4 <= 0.001 & burnt.pig.vec$V3 >= 0.5)
#finger pigment vectors
finger.pig.vec <- envfit(finger.comb.mds,env=finger.master.dat[,32:51],na.rm=T)
finger.pig.vec <- data.frame(cbind(finger.pig.vec$vectors$arrows, 
                                  finger.pig.vec$vectors$r, 
                                  finger.pig.vec$vectors$pvals))
finger.pig.vec.sig <- subset(finger.pig.vec, finger.pig.vec$V4 <= 0.001 & finger.pig.vec$V3 >= 0.5)
#wtwin pigment vectors
wtwin.pig.vec <- envfit(wtwin.comb.mds,env=wtwin.master.dat[,32:51],na.rm=T)
wtwin.pig.vec <- data.frame(cbind(wtwin.pig.vec$vectors$arrows, 
                               wtwin.pig.vec$vectors$r, 
                               wtwin.pig.vec$vectors$pvals))
wtwin.pig.vec.sig <- subset(wtwin.pig.vec, wtwin.pig.vec$V4 <= 0.001 & wtwin.pig.vec$V3 >= 0.5)

#create PCoA plots
burnt.pcoa <- 
ggplot(aes(Dim1, Dim2), data=burnt.comb.scores) +  
  geom_segment(aes(x=rep(0, nrow(burnt.geo.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(burnt.geo.vec.sig)),yend=Dim2),
               colour="hotpink", data=burnt.geo.vec.sig) +
  geom_segment(aes(x=rep(0, nrow(burnt.diat.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(burnt.diat.vec.sig)),yend=Dim2),
               colour="orange", data=burnt.diat.vec.sig) +
  geom_segment(aes(x=rep(0, nrow(burnt.pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(burnt.pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=burnt.pig.vec.sig) + 
  geom_text(aes(label=round(burnt.comb.scores$year,0)),size=3,position=position_jitter(width=0.1,height=0.1)) + 
  geom_path(aes(Dim1, Dim2), arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), lineend="round") + 
  geom_text(aes(x=Dim1, y=Dim2), data=burnt.geo.vec.sig, label=rownames(burnt.geo.vec.sig), colour="hotpink", size=4) +
  geom_text(aes(x=Dim1, y=Dim2), data=burnt.diat.vec.sig, label=rownames(burnt.diat.vec.sig), colour="orange", size=4) +
  geom_text(aes(x=Dim1, y=Dim2), data=burnt.pig.vec.sig, label=rownames(burnt.pig.vec.sig), colour="darkgreen", size=4) + 
  scale_x_continuous(limits=c(-2,1.5),breaks=seq(-2,1.5,0.5))+
  scale_y_continuous(limits=c(-1.5,2))+
  ggtitle("Burnt Lake - Polymictic")+
  theme_bw(base_size=14)
finger.pcoa <- 
ggplot(aes(Dim1, Dim2), data=finger.comb.scores) +  
  geom_segment(aes(x=rep(0, nrow(finger.geo.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(finger.geo.vec.sig)),yend=Dim2),
               colour="hotpink", data=finger.geo.vec.sig) +
  geom_segment(aes(x=rep(0, nrow(finger.diat.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(finger.diat.vec.sig)),yend=Dim2),
               colour="orange", data=finger.diat.vec.sig) +
  geom_segment(aes(x=rep(0, nrow(finger.pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(finger.pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=finger.pig.vec.sig) + 
  geom_text(aes(label=round(finger.comb.scores$year,0)),size=3,position=position_jitter(width=0.1,height=0.1)) + 
  geom_path(aes(Dim1, Dim2), arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), lineend="round") + 
  geom_text(aes(x=Dim1, y=Dim2), data=finger.geo.vec.sig, label=rownames(finger.geo.vec.sig), colour="hotpink", size=4) +
  geom_text(aes(x=Dim1, y=Dim2), data=finger.diat.vec.sig, label=rownames(finger.diat.vec.sig), colour="orange", size=4) +
  geom_text(aes(x=Dim1, y=Dim2), data=finger.pig.vec.sig, label=rownames(finger.pig.vec.sig), colour="darkgreen", size=4) + 
  scale_x_continuous(limits=c(-2,1.5),breaks=seq(-2,1.5,0.5))+
  scale_y_continuous(limits=c(-1.5,2))+
  ggtitle("Finger Lake - Well mixed") +
  theme_bw(base_size=14)
wtwin.pcoa <- 
ggplot(aes(Dim1, Dim2), data=wtwin.comb.scores) +  
  geom_segment(aes(x=rep(0, nrow(wtwin.geo.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(wtwin.geo.vec.sig)),yend=Dim2),
               colour="hotpink", data=wtwin.geo.vec.sig) +
  geom_segment(aes(x=rep(0, nrow(wtwin.diat.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(wtwin.diat.vec.sig)),yend=Dim2),
               colour="orange", data=wtwin.diat.vec.sig) +
  geom_segment(aes(x=rep(0, nrow(wtwin.pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(wtwin.pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=wtwin.pig.vec.sig) + 
  geom_text(aes(label=round(wtwin.comb.scores$year,0)),size=3,position=position_jitter(width=0.1,height=0.1)) + 
  geom_path(aes(Dim1, Dim2), arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), lineend="round") + 
  geom_text(aes(x=Dim1, y=Dim2), data=wtwin.geo.vec.sig, label=rownames(wtwin.geo.vec.sig), colour="hotpink", size=4) +
  geom_text(aes(x=Dim1, y=Dim2), data=wtwin.diat.vec.sig, label=rownames(wtwin.diat.vec.sig), colour="orange", size=4) +
  geom_text(aes(x=Dim1, y=Dim2), data=wtwin.pig.vec.sig, label=rownames(wtwin.pig.vec.sig), colour="darkgreen", size=4) + 
  scale_x_continuous(limits=c(-2,1.5),breaks=seq(-2,1.5,0.5))+
  scale_y_continuous(limits=c(-1.5,2))+
  ggtitle("West Twin - Dimictic") +
  theme_bw(base_size=14)

cowplot::plot_grid(burnt.pcoa,finger.pcoa,wtwin.pcoa,ncol=3)
