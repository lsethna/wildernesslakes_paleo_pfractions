librarian::shelf(tidyverse,grid,vegan,cowplot,zoo)

rm(list=ls())

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/Paleo PCoA")

#read in master data
master_dat <- read.csv("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/raw_data/WL_paleo_masterdataset_14Oct2024.csv")

## ------------------------------------------------------------ ##
## ---- filter variables: pigments>0, diatoms>2% abundance ---- ##
## ------------------------------------------------------------ ##
vars <- data.frame(col.num=c(1:272),
                   var=colnames(master_dat))
#pigments 25:56
#diatoms 57:270
sel_d <- master_dat %>% select(colnames(master_dat[,57:270])) %>% 
  summarize(across(ach_microcephala:uln_ulna,max,na.rm=T)) #get maximum count for each species
sel_diat <- pivot_longer(sel_d,cols=ach_microcephala:uln_ulna) %>% filter(value>8) #filter to max value 8 (8 counts=2% rel. abund)
#16 total diatoms with >2% rel abund
sel_p <- master_dat %>% select(colnames(master_dat[,25:56])) %>% 
  summarize(across(chlide_a:car_z,max,na.rm=T))
sel_pig <- pivot_longer(sel_p,cols=chlide_a:car_z) %>% filter(value>0)
#20 different detectable pigments
#filter master data based on selected diatoms and pigments
master_v2 <- master_dat %>% select(lake,year_loess,depth, #id info
                                   10:12, #loi; perc inorg,org,calc
                                   17:24, #P and P fractions
                                   13,sel_diat$name, #BSi conc, select diatoms
                                   sel_pig$name)
glimpse(master_v2)

## ------------------------------- ##
## ---- interpolate variables ---- ##
## ------------------------------- ##
master_v3_interp <- master_v2 %>%
  group_by(lake) %>%
  mutate_at(vars(ex_p:b_car),funs(zoo::na.approx(.,method="constant",rule=2))) %>% #rule=2 means extend nearest values to leading/trailing NAs
  ungroup()
colnames(master_v3_interp)

write.csv(master_v3_interp,file="raw_data/interpolated_master_dat_15Oct24.csv")

## -------------- ##
## ---- PCoA ---- ##
## -------------- ##

#diatoms=16:31
#pigments=32:51

#standardize diatoms
diat.stand.hell <- decostand(master_v3_interp[,16:31], method="hellinger")
#diatom distance matrix
diat.dist <- dist(diat.stand.hell, method="euclidean")
diat.dist <- diat.dist/max(diat.dist)

#standardize pigments
pig.stand.hell <- decostand(master_v3_interp[,32:51], method="hellinger")
#pigment distance matrix
pig.dist <- dist(pig.stand.hell, method="euclidean")
pig.dist <- pig.dist/max(pig.dist)

#combine distance matrices
comb.dist <- diat.dist+pig.dist
#and scale
comb.mds <- cmdscale(d=comb.dist, k=2)
#add site IDs
comb.mds$lake <- master_v3_interp$lake
comb.mds$year <- master_v3_interp$year_loess
plot(comb.mds)
#get scores
site.comb <- data.frame(scores(comb.mds))
#add sample IDs
site.comb$lake <- master_v3_interp$lake
site.comb$year <- master_v3_interp$year_loess
site.comb$mix.regime <- master_v3_interp$mix.regime

#combined vectors for all vars
comb.vec <- envfit(comb.mds, env=master_v3_interp[,-3],na.rm=T) #include all paleo vars except depth
comb.vec <- data.frame(cbind(comb.vec$vectors$arrows, comb.vec$vectors$r, comb.vec$vectors$pvals))
comb.vec.sig <- subset(comb.vec, comb.vec$V4 <= 0.001 & comb.vec$V3 >= 0.5)
#just diatom vectors
diat.vec <- envfit(comb.mds,env=master_v3_interp[,16:31],na.rm=T)
diat.vec <- data.frame(cbind(diat.vec$vectors$arrows, diat.vec$vectors$r, diat.vec$vectors$pvals))
diat.vec.sig <- subset(diat.vec, diat.vec$V4 <= 0.001 & diat.vec$V3 >= 0.5)
#just pigment vectors
pig.vec <- envfit(comb.mds,env=master_v3_interp[,32:51],na.rm=T)
pig.vec <- data.frame(cbind(pig.vec$vectors$arrows, pig.vec$vectors$r, pig.vec$vectors$pvals))
pig.vec.sig <- subset(pig.vec, pig.vec$V4 <= 0.001 & pig.vec$V3 >= 0.5)

ggplot(aes(Dim1, Dim2,colour=lake), data=site.comb) +  
  #geom_text(aes(label=round(site.comb$year,0)), vjust=2, size=3, fontface="bold", show_guide=F) + 
  geom_path(aes(Dim1, Dim2),size=1, 
            arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), 
            lineend="round") + 
  theme_bw(base_size=14)

#plot PCoA paths for all lakes
ggplot(aes(Dim1, Dim2),data=site.comb) +  
  geom_text(aes(label=round(site.comb$year,0)), vjust=2, size=3, fontface="bold", show_guide=F) + 
  geom_path(aes(Dim1, Dim2,lty=lake), arrow=arrow(type="closed", ends="first",length=unit(0.05,"inches")), lineend="round") + 
  geom_segment(aes(x=rep(0, nrow(pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=pig.vec.sig) + 
  geom_text(aes(x=Dim1, y=Dim2), data=pig.vec.sig, label=rownames(pig.vec.sig), colour="darkgreen", size=4) +  
  geom_segment(aes(x=rep(0, nrow(diat.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(diat.vec.sig)),yend=Dim2), 
               colour="orange", data=diat.vec.sig) + 
  geom_text(aes(x=Dim1, y=Dim2), data=diat.vec.sig, label=rownames(diat.vec.sig), colour="orange", size=4) +
  theme_bw(base_size=14)

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
#create new dfs for each mixing regime
poly.master.dat <- subset(master_v3_interp,master_v3_interp$mix.regime=="Polymictic")
mixed.master.dat <- subset(master_v3_interp,master_v3_interp$mix.regime=="Well mixed")
di.master.dat <- subset(master_v3_interp,master_v3_interp$mix.regime=="Dimictic")

#standardize diatoms
diat.stand.hell.poly <- decostand(poly.master.dat[,8:23], method="hellinger")
diat.stand.hell.mixed <- decostand(mixed.master.dat[,8:23], method="hellinger")
diat.stand.hell.di <- decostand(di.master.dat[,8:23], method="hellinger")
#diatom distance matrix
diat.dist.poly <- dist(diat.stand.hell.poly, method="euclidean")
diat.dist.mixed <- dist(diat.stand.hell.mixed, method="euclidean")
diat.dist.di <- dist(diat.stand.hell.di, method="euclidean")

diat.dist.poly <- diat.dist.poly/max(diat.dist.poly)
diat.dist.mixed <- diat.dist.mixed/max(diat.dist.mixed)
diat.dist.di <- diat.dist.di/max(diat.dist.di)

#standardize pigments
pig.stand.hell.poly <- decostand(poly.master.dat[,24:43], method="hellinger")
pig.stand.hell.mixed <- decostand(mixed.master.dat[,24:43], method="hellinger")
pig.stand.hell.di <- decostand(di.master.dat[,24:43], method="hellinger")
#pigment distance matrix
pig.dist.poly <- dist(pig.stand.hell.poly, method="euclidean")
pig.dist.mixed <- dist(pig.stand.hell.mixed, method="euclidean")
pig.dist.di <- dist(pig.stand.hell.di, method="euclidean")

pig.dist.poly <- pig.dist.poly/max(pig.dist.poly)
pig.dist.mixed <- pig.dist.mixed/max(pig.dist.mixed)
pig.dist.di <- pig.dist.di/max(pig.dist.di)

#combine distance matrices
poly.comb.dist <- diat.dist.poly+pig.dist.poly
mixed.comb.dist <- diat.dist.mixed+pig.dist.mixed
di.comb.dist <- diat.dist.di+pig.dist.di

#and scale
poly.comb.mds <- cmdscale(d=poly.comb.dist, k=2)
mixed.comb.mds <- cmdscale(d=mixed.comb.dist, k=2)
di.comb.mds <- cmdscale(d=di.comb.dist, k=2)
#plot
plot(poly.comb.mds)
plot(mixed.comb.mds)
plot(di.comb.mds)
#get scores
poly.comb.scores <- data.frame(scores(poly.comb.mds))
mixed.comb.scores <- data.frame(scores(mixed.comb.mds))
di.comb.scores <- data.frame(scores(di.comb.mds))

#add site IDs
poly.comb.scores$lake <- poly.master.dat$lake
poly.comb.scores$year <- poly.master.dat$year_loess
mixed.comb.scores$lake <- mixed.master.dat$lake
mixed.comb.scores$year <- mixed.master.dat$year_loess
di.comb.scores$lake <- di.master.dat$lake
di.comb.scores$year <- di.master.dat$year_loess

#diatom and pigment vectors for each mixing regime
#polymictic diatom vectors
poly.diat.vec <- envfit(poly.comb.mds,env=poly.master.dat[,8:23],na.rm=T)
poly.diat.vec <- data.frame(cbind(poly.diat.vec$vectors$arrows, 
                                  poly.diat.vec$vectors$r, 
                                  poly.diat.vec$vectors$pvals))
poly.diat.vec.sig <- subset(poly.diat.vec, poly.diat.vec$V4 <= 0.001 & poly.diat.vec$V3 >= 0.5)
#well mixed diatom vectors
mixed.diat.vec <- envfit(mixed.comb.mds,env=mixed.master.dat[,8:23],na.rm=T)
mixed.diat.vec <- data.frame(cbind(mixed.diat.vec$vectors$arrows, 
                                  mixed.diat.vec$vectors$r, 
                                  mixed.diat.vec$vectors$pvals))
mixed.diat.vec.sig <- subset(mixed.diat.vec, mixed.diat.vec$V4 <= 0.001 & mixed.diat.vec$V3 >= 0.5)
#dimictic diatom vectors
di.diat.vec <- envfit(di.comb.mds,env=di.master.dat[,8:23],na.rm=T)
di.diat.vec <- data.frame(cbind(di.diat.vec$vectors$arrows, 
                                  di.diat.vec$vectors$r, 
                                  di.diat.vec$vectors$pvals))
di.diat.vec.sig <- subset(di.diat.vec, di.diat.vec$V4 <= 0.001 & di.diat.vec$V3 >= 0.5)
#polymictic pigment vectors
poly.pig.vec <- envfit(poly.comb.mds,env=poly.master.dat[,24:43],na.rm=T)
poly.pig.vec <- data.frame(cbind(poly.pig.vec$vectors$arrows, 
                                  poly.pig.vec$vectors$r, 
                                  poly.pig.vec$vectors$pvals))
poly.pig.vec.sig <- subset(poly.pig.vec, poly.pig.vec$V4 <= 0.001 & poly.pig.vec$V3 >= 0.5)
#well mixed pigment vectors
mixed.pig.vec <- envfit(mixed.comb.mds,env=mixed.master.dat[,24:43],na.rm=T)
mixed.pig.vec <- data.frame(cbind(mixed.pig.vec$vectors$arrows, 
                                   mixed.pig.vec$vectors$r, 
                                   mixed.pig.vec$vectors$pvals))
mixed.pig.vec.sig <- subset(mixed.pig.vec, mixed.pig.vec$V4 <= 0.001 & mixed.pig.vec$V3 >= 0.5)
#dimictic pigment vectors
di.pig.vec <- envfit(di.comb.mds,env=di.master.dat[,24:43],na.rm=T)
di.pig.vec <- data.frame(cbind(di.pig.vec$vectors$arrows, 
                                di.pig.vec$vectors$r, 
                                di.pig.vec$vectors$pvals))
di.pig.vec.sig <- subset(di.pig.vec, di.pig.vec$V4 <= 0.001 & di.pig.vec$V3 >= 0.5)

#create PCoA plots
#poly.pcoa <- 
ggplot(aes(Dim1, Dim2, colour=lake), data=poly.comb.scores) +  
  #geom_text(aes(label=round(poly.comb.scores$year,0)),size=3, fontface="bold",position=position_jitter(width=0.1,height=0.1)) + 
  geom_path(aes(Dim1, Dim2), arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), lineend="round") + 
  geom_segment(aes(x=rep(0, nrow(poly.diat.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(poly.diat.vec.sig)),yend=Dim2),
               colour="orange", data=poly.diat.vec.sig) +
  geom_text(aes(x=Dim1, y=Dim2), data=poly.diat.vec.sig, label=rownames(poly.diat.vec.sig), colour="orange", size=4) +
  geom_segment(aes(x=rep(0, nrow(poly.pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(poly.pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=poly.pig.vec.sig) + 
  geom_text(aes(x=Dim1, y=Dim2), data=poly.pig.vec.sig, label=rownames(poly.pig.vec.sig), colour="darkgreen", size=4) + 
  theme_bw(base_size=14)
#mixed.pcoa <- 
ggplot(aes(Dim1, Dim2, colour=lake), data=mixed.comb.scores) +  
  #geom_text(aes(label=round(mixed.comb.scores$year,0)),size=3, fontface="bold",position=position_jitter(width=0.1,height=0.1)) + 
  geom_path(aes(Dim1, Dim2), arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), lineend="round") + 
  geom_segment(aes(x=rep(0, nrow(mixed.diat.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(mixed.diat.vec.sig)),yend=Dim2),
               colour="orange", data=mixed.diat.vec.sig) +
  geom_text(aes(x=Dim1, y=Dim2), data=mixed.diat.vec.sig, label=rownames(mixed.diat.vec.sig), colour="orange", size=4) +
  geom_segment(aes(x=rep(0, nrow(mixed.pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(mixed.pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=mixed.pig.vec.sig) + 
  geom_text(aes(x=Dim1, y=Dim2), data=mixed.pig.vec.sig, label=rownames(mixed.pig.vec.sig), colour="darkgreen", size=4) + 
  theme_bw(base_size=14)
#di.pcoa <- 
ggplot(aes(Dim1, Dim2, colour=lake), data=di.comb.scores) +  
  #geom_text(aes(label=round(di.comb.scores$year,0)),size=3, fontface="bold",position=position_jitter(width=0.1,height=0.1)) + 
  geom_path(aes(Dim1, Dim2), arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), lineend="round") + 
  geom_segment(aes(x=rep(0, nrow(di.diat.vec.sig)),
                   xend=Dim1, y=rep(0, nrow(di.diat.vec.sig)),yend=Dim2),
               colour="orange", data=di.diat.vec.sig) +
  geom_text(aes(x=Dim1, y=Dim2), data=di.diat.vec.sig, label=rownames(di.diat.vec.sig), colour="orange", size=4) +
  geom_segment(aes(x=rep(0, nrow(di.pig.vec.sig)), 
                   xend=Dim1, y=rep(0, nrow(di.pig.vec.sig)),yend=Dim2), 
               colour="darkgreen", data=di.pig.vec.sig) + 
  geom_text(aes(x=Dim1, y=Dim2), data=di.pig.vec.sig, label=rownames(di.pig.vec.sig), colour="darkgreen", size=4) + 
  theme_bw(base_size=14)

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
