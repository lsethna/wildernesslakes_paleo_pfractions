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
#pigments 24:55
#diatoms 56:269
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
                                   #9:11, #loi; perc inorg,org,calc
                                   #16:23, #P and P fractions
                                   #14, #BSi flux
                                   sel_diat$name, #select diatoms
                                   sel_pig$name, #select pigments
                                   #274:276 #TOC, TN, TOC:TN

)  
glimpse(master_v2)

## ------------------------------- ##
## ---- interpolate variables ---- ##
## ------------------------------- ##
master_v3_interp <- master_v2 %>%
  group_by(lake) %>%
  #will replace NA values with interpolated values - done at every dated interval 
  mutate_at(vars(ast_formosa:b_car),funs(zoo::na.approx(.,method="constant",rule=2))) %>% #rule=2 means extend nearest values to leading/trailing NAs
  ungroup() %>%
  drop_na() #remove any NA values - beyond dated intervals

#quick check to make sure all dates and dates are included
ggplot(master_v3_interp,aes(x=year_loess,y=depth))+
  geom_line()+facet_wrap(~lake)

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions")
write.csv(master_v3_interp,file="raw_data/interpolated_diatom_pigment_dat_1Dec25.csv")

colnames(master_v3_interp)
#diatoms=4:29
#pigments=30:41

diat_pig_vars <- colnames(master_v3_interp)[4:41]

## ------------------------------------------------------- ##
## ---- create new variable names for pretty plotting ---- ##
## ------------------------------------------------------- ##
pretty_names <- function(x) {
  parts <- str_split_fixed(x, "_", 3)  # split each name into up to 3 parts
  
  pre  <- str_sub(parts[,1], 1, 3)                  # first 3 chars of first part
  mid  <- ifelse(parts[,2] != "", str_sub(parts[,2], 1, 2), "")  # first 2 chars of second part (if any)
  post <- ifelse(parts[,3] != "", str_sub(parts[,3], 1, 1), "")  # first char of third part (if any)
  
  str_to_upper(paste0(pre, mid, post))
}

pretty_name_code <- tibble(var = colnames(master_v3_interp)) %>% 
  mutate(code = case_when(var %in% diat_pig_vars ~ pretty_names(var), #change only diatoms and pigments
                          T ~ var)) %>%
  #manually change BCA to BCAR
  mutate(code = case_when(code=="BCA" ~ "BCAR",
                          T ~ code))

#check to make sure all codes are unique
length(unique(pretty_name_code$code))
#export
write.csv(pretty_name_code,file="pcoa_variable_abbr.csv")

###change master data to have pretty column names
name_map <- setNames(pretty_name_code$code, pretty_name_code$var)
# rename columns in the dataset
master_v3_interp <- master_v3_interp %>%
  rename_with(~ name_map[.x], .cols = names(name_map))


## --------------------------------------------------------------------------- ##
## -------------------------- PCoA all lakes together ------------------------ ##
## --------------------------------------------------------------------------- ##

#check column names and numbers
colnames(master_v3_interp)
#diatoms 4:29
#pigments 30:41

## Calculate distance matrices

#standardize diatoms
diat.stand.hell <- vegan::decostand(master_v3_interp[,4:29], method="hellinger")
#diatom distance matrix
diat.dist <- dist(diat.stand.hell, method="euclidean")
#normalize the diatom vectors to allow for comparison across proxies
diat.dist <- diat.dist/max(diat.dist)

#standardize pigments
pig.stand.hell <- decostand(master_v3_interp[,30:41], method="hellinger")
#pigment distance matrix
pig.dist <- dist(pig.stand.hell, method="euclidean")
pig.dist <- pig.dist/max(pig.dist)

#combine distance matrices
comb.dist <- diat.dist+pig.dist

## Principal Coordinates Analysis
pcoa_res <- cmdscale(comb.dist, k = 2, eig = TRUE)   # eig = TRUE for variance explained
site_scores <- as.data.frame(pcoa_res$points)
plot(site_scores)

colnames(site_scores) <- c("PCoA1", "PCoA2")

# add metadata
site_scores$lake <- master_v3_interp$lake
site_scores$year <- master_v3_interp$year_loess

## Extract variable loadings

# Diatoms
diat.vec <- envfit(pcoa_res$points, master_v3_interp[,4:29], na.rm = TRUE)
diat.df <- data.frame(
  diat.vec$vectors$arrows,
  r2 = diat.vec$vectors$r,
  p = diat.vec$vectors$pvals,
  varname = colnames(master_v3_interp)[4:29],
  type = "diatom"
)

# Pigments
pig.vec <- envfit(pcoa_res$points, master_v3_interp[,30:41], na.rm = TRUE)
pig.df <- data.frame(
  pig.vec$vectors$arrows,
  r2 = pig.vec$vectors$r,
  p = pig.vec$vectors$pvals,
  varname = colnames(master_v3_interp)[30:41],
  type = "pigment"
)

# Combine & apply significance filter
vecs <- bind_rows(diat.df, pig.df) %>%
  filter(p <= 0.05, r2 >= 0.5)

# Scale arrows by r² to indicate correlation strength
vecs_scaled <- vecs %>%
  mutate(
    Dim1 = Dim1 * r2,
    Dim2 = Dim2 * r2
  )
#these are the loadings for each of the variables

#plot PCoA paths for all lakes
ggplot(aes(PCoA1, PCoA2),data=site_scores) +  
  #lake trajectories 
  geom_path(aes(color=lake), size=0.75, 
            arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), 
            lineend="round") + 
  #vectors
  geom_segment(data = vecs_scaled,
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               color = "gray40") +
  geom_point(data = vecs_scaled,
             aes(x = Dim1, y = Dim2, fill = type),
             size = 3, pch = 21) +
  geom_text(data = vecs_scaled,
            aes(x = Dim1, y = Dim2, label = varname),
            size = 4,
            nudge_x = 0.05, nudge_y = 0.05) +
  #general aesthetics
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Lake",fill="Variable type") +
  scale_x_continuous(limits=c(-0.75,0.75))+
  scale_y_continuous(limits=c(-0.75,0.75))+
  scale_color_manual(values=c("#C1D9B7", #burnt
                              "#89CFF0", #dunnigan
                              "#007CAA", #elbow
                              "#665191", #etwin
                              "#d45087", #finger
                              "#7A871E", #flame
                              "#ebdcff", #smoke
                              "#ffd8c9" #wtwin
  )) +
  scale_fill_manual(values=c("gold","darkgreen")) +
  theme_classic(base_size=14) 


## ---------------------------------------------------------------------------- ##
## --------------------- PCoA for each mixing regime -------------------------- ##
## ---------------------------------------------------------------------------- ##

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
pcoa.plots <- list() #list to save PCoA plots

for (i in 1:length(mixing_regime)) {
  #create df for each mixing regime
  mix.master.dat <- subset(master_v3_interp,master_v3_interp$mix.regime==mixing_regime[i])
  
  #standardize diatoms and pigments... 
  diat.stand.hell.mix <- decostand(mix.master.dat[,4:29], method="hellinger")
  pig.stand.hell.mix <- decostand(mix.master.dat[,30:41], method="hellinger")
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
  diat.vec.mix <- envfit(comb.mds.mix,env=mix.master.dat[,4:29],na.rm=T)
  diat.vec.mix <- data.frame(cbind(diat.vec.mix$vectors$arrows, 
                                   diat.vec.mix$vectors$r, 
                                   diat.vec.mix$vectors$pvals),
                             var_type="diatoms")
  #pigments
  pig.vec.mix <- envfit(comb.mds.mix,env=mix.master.dat[,30:41],na.rm=T)
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
                                                                  lake=="finger" ~ "#d45087",
                                                                  lake=="flame" ~ "#7A871E",
                                                                  lake=="smoke" ~ "#ebdcff",
                                                                  lake=="wtwin" ~ "#ffd8c9")) 

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

cowplot::plot_grid(pcoa.plots[[2]],pcoa.plots[[1]],pcoa.plots[[3]],ncol=3)

## --------------------------------------------------------------------------- ##
## -------------------------- PCoA for each lake ----------------------------- ##
## --------------------------------------------------------------------------- ##

#set up variables to run PcoA for each lake in a loop
unique(master_v3_interp$lake)
#put in order for plotting
lakes = c("dunnigan","finger","burnt","smoke","elbow","etwin","flame","wtwin")
#make list of plot titles
lake_plot_titles <- c("Dunnigan","Finger","Burnt","Smoke","Elbow","East Twin","Flame","West Twin")

pcoa.scores_lake <- list() #list to save scores output
pcoa.variable.vectors.lake <- list() #list to save variable vectors
pcoa.plots.lake <- list() #list to save PCoA plots

for (i in 1:length(lakes)) {
  #create df for each mixing regime
  lake.master.dat <- subset(master_v3_interp,master_v3_interp$lake==lakes[i])
  
  #standardize diatoms and pigments... 
  diat.stand.hell.lake <- decostand(lake.master.dat[,4:29], method="hellinger")
  pig.stand.hell.lake <- decostand(lake.master.dat[,30:41], method="hellinger")
  #...to create distance matrices
  diat.dist.lake <- dist(diat.stand.hell.lake, method="euclidean")
  pig.dist.lake <- dist(pig.stand.hell.lake, method="euclidean")
  
  #combine distance matrices...
  comb.dist.lake <- diat.dist.lake+pig.dist.lake
  #...and scale
  comb.mds.lake <- cmdscale(d=comb.dist.lake, k=2)
  
  #get scores
  comb.scores.lake <- data.frame(scores(comb.mds.lake))
  #add sample IDs
  comb.scores.lake$lake <- lake.master.dat$lake
  comb.scores.lake$year <- lake.master.dat$year_loess
  
  #save to list
  pcoa.scores_lake[[i]] <- comb.scores.lake
  
  #calculate vectors for diatoms and pigments
  #diatom vectors
  diat.vec.lake <- envfit(comb.mds.lake,env=lake.master.dat[,4:29],na.rm=T)
  diat.vec.lake <- data.frame(cbind(diat.vec.lake$vectors$arrows, 
                                   diat.vec.lake$vectors$r, 
                                   diat.vec.lake$vectors$pvals),
                             var_type="diatoms")
  #pigments
  pig.vec.lake <- envfit(comb.mds.lake,env=lake.master.dat[,30:41],na.rm=T)
  pig.vec.lake <- data.frame(cbind(pig.vec.lake$vectors$arrows, 
                                  pig.vec.lake$vectors$r, 
                                  pig.vec.lake$vectors$pvals),
                            var_type="pigment")
  
  #keep only significant vectors for plotting
  diat.vec.lake.sig <- subset(diat.vec.lake, diat.vec.lake$V4 <= 0.05 & diat.vec.lake$V3 >= 0.5)
  pig.vec.lake.sig <- subset(pig.vec.lake, pig.vec.lake$V4 <= 0.05 & pig.vec.lake$V3 >= 0.5)
  
  #save these to list
  vec.lake <- rbind(diat.vec.lake.sig,pig.vec.lake.sig) %>% mutate(lake=lakes[i])
  pcoa.variable.vectors.lake[[i]] <- vec.lake
  
  # Scale arrows by r2 to show correlation strength
  diat.vec.lake.sig.scaled <- diat.vec.lake.sig %>%
    mutate(Dim1 = (Dim1*V3)/2,
           Dim2 = (Dim2*V3)/2)
  pig.vec.lake.sig.scaled <- pig.vec.lake.sig %>%
    mutate(Dim1 = (Dim1*V3)/2,
           Dim2 = (Dim2*V3)/2)
  
  #create PCoA biplot

  p <- 
    ggplot(aes(Dim1, Dim2),data=comb.scores.lake) +  
    #lake trajectories 
    geom_path(color="black", size=0.75, 
              arrow=arrow(type="closed", ends="first",length=unit(0.05,"inches")), lineend="round") + 
    #pigment vectors
    geom_segment(data=pig.vec.lake.sig.scaled,
                 aes(x=0,y=0,xend=Dim1,yend=Dim2), 
                 color="darkgreen") + 
    geom_text(data=pig.vec.lake.sig.scaled,
              aes(x=Dim1, y=Dim2),label=rownames(pig.vec.lake.sig.scaled), 
              color="darkgreen", size=3,vjust=-0.5,
              position=position_jitter(width=0.1,height=0.1)) +  
    #diatom vectors
    geom_segment(data=diat.vec.lake.sig.scaled,
                 aes(x=0,y=0,xend=Dim1,yend=Dim2), 
                 color="orange") + 
    geom_text(data=diat.vec.lake.sig.scaled,
              aes(x=Dim1, y=Dim2),label=rownames(diat.vec.lake.sig.scaled), 
              color="orange", size=3,vjust=-0.5,
              position=position_jitter(width=0.1,height=0.1)) +
    #general aesthetics
    ggtitle(lake_plot_titles[i]) +
    theme_classic()+
    theme(axis.title=element_blank(),
          title=element_text(size=11),
          axis.text=element_text(size=10))
  
  pcoa.plots[[i]] <- p
}

cowplot::plot_grid(plotlist=pcoa.plots,ncol=2,
                   labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),
                   label_fontface="plain")

#move rownames (important variables) to a column
pcoa.variable.vectors.lake <- lapply(pcoa.variable.vectors.lake,tibble::rownames_to_column, var="var")
pcoa.sig.variables.lake <- do.call(rbind,lapply(pcoa.variable.vectors.lake,as.data.frame))
unique(pcoa.sig.variables.lake$lake)
write.csv(pcoa.sig.variables.lake,file="sig_variables_pcoa_by_lake.csv")

## plot up PCoA loadings over time for each lake
pcoa.scores_lake_df <- dplyr::bind_rows(pcoa.scores_lake)
glimpse(pcoa.scores_lake_df)

#round year values to nearest 5, take average?
round_to_nearest_10 <- function(x) {
  10 * round(x / 10)
}
#test
round_to_nearest_10(1986)

pcoa.scores_lake_df %>%
  mutate(year_round = round_to_nearest_10(year)) %>% 
  group_by(lake,year_round) %>%
  summarize(Dim1 = median(Dim1)) %>%
ggplot(aes(x=year_round,y=Dim1,color=lake))+
  geom_vline(xintercept = 1900,lty="dashed")+
  geom_vline(xintercept = 1970,lty="dashed")+
  #geom_point()+
  geom_smooth(se=F)+
  #geom_line(size=1)+
  theme_classic()
