librarian::shelf(tidyverse,grid,vegan,cowplot,zoo)

rm(list=ls())

setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/Paleo PCoA")

#read in master data
master_dat <- read.csv("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions/raw_data/WL_paleo_masterdataset_20Nov2024.csv") %>%
  select(!X)
#Amelia's version
master_dat <- read_csv("~/Desktop/GitHub/wildernesslakes_paleo_pfractions/Master dataset/WL_paleo_masterdataset_20Nov2024.csv") 
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
sel_diat <- master_dat %>% select(colnames(master_dat[,c(56:270)])) %>% 
  pivot_longer(cols=ach_microcephala:uln_ulna) %>% filter(value>5) %>% #filter to 5% relative abundance (might want to only filter to 2% and at least in two samples)
  select(name) %>% distinct()
  
##
sel_p <- master_dat %>% select(colnames(master_dat[,24:56])) %>% 
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
  #will replace NA values with interpolated values - done at every dated interval 
  mutate_at(vars(percent_organic:b_car),funs(zoo::na.approx(.,method="constant",rule=2))) %>% #rule=2 means extend nearest values to leading/trailing NAs
  ungroup() %>%
  drop_na() #remove any NA values - beyond dated intervals

#quick check to make sure all dates and dates are included
ggplot(master_v3_interp,aes(x=year_loess,y=depth))+
  geom_line()+facet_wrap(~lake)

#to pull p-frac data
setwd("C:/Users/lsethna_smm/Documents/GitHub/wildernesslakes_paleo_pfractions")
write.csv(master_v3_interp,file="raw_data/interpolated_master_dat_4Sep25.csv")

## --------------------------------------------------------------------------- ##
## -------------------------- PCoA all lakes together ------------------------ ##
## --------------------------------------------------------------------------- ##
###### To make names simple--> A skipped######
colnames(master_v3_interp)
#diatoms=8:33
#pigments=34:45
#geochem=4:7

##create new variable names for pretty plotting
pretty_names <- function(x) {
  parts <- str_split_fixed(x, "_", 3)  # split each name into up to 3 parts
  
  pre  <- str_sub(parts[,1], 1, 3)                  # first 3 chars of first part
  mid  <- ifelse(parts[,2] != "", str_sub(parts[,2], 1, 2), "")  # first 2 chars of second part (if any)
  post <- ifelse(parts[,3] != "", str_sub(parts[,3], 1, 1), "")  # first char of third part (if any)
  
  str_to_upper(paste0(pre, mid, post))
}

pretty_name_code <- tibble(var = colnames(master_v3_interp)) %>% 
  mutate(code = case_when(row_number() %in% 8:45 ~ pretty_names(var), #change only diatoms and pigments
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

####### End of pretty names 


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
comb.mds <- cmdscale(d = comb.dist, k = 2, eig = TRUE)

#extract site scores
comb.scores <- as.data.frame(comb.mds$points)


#to calc varence
eigs <- comb.mds$eig
var1 <- round((eigs[1] / sum(eigs[eigs > 0])) * 100, 1)
var2 <- round((eigs[2] / sum(eigs[eigs > 0])) * 100, 1)

cat("Axis 1 =", var1, "%, Axis 2 =", var2, "% of total variance\n")



comb.mds <- cmdscale(d = comb.dist, k = 2)
plot(comb.mds)

#get scores
site.comb <- data.frame(scores(comb.mds))

#add sample IDs
site.comb$lake <- master_v3_interp$lake
site.comb$year <- master_v3_interp$year_loess

glimpse(site.comb)
#diatom vectors
diat.vec <- envfit(comb.mds,env=master_v3_interp[,8:33],na.rm=T)
diat.vec <- data.frame(cbind(diat.vec$vectors$arrows, diat.vec$vectors$r, diat.vec$vectors$pvals))
diat.vec.sig <- subset(diat.vec, diat.vec$V4 <= 0.05 & diat.vec$V3 >= 0.5) #only saves vectors that are significant (p<0.05) and are >50% correlated with ordination
#pigment vectors
pig.vec <- envfit(comb.mds,env=master_v3_interp[,34:45],na.rm=T)
pig.vec <- data.frame(cbind(pig.vec$vectors$arrows, pig.vec$vectors$r, pig.vec$vectors$pvals))
pig.vec.sig <- subset(pig.vec, pig.vec$V4 <= 0.05 & pig.vec$V3 >= 0.5)
#just geochem vectors
# geochem.vec <- envfit(comb.mds,env=master_v3_interp[,c(4:7)],na.rm=T)
# geochem.vec <- data.frame(cbind(geochem.vec$vectors$arrows, geochem.vec$vectors$r, geochem.vec$vectors$pvals))
# geochem.vec.sig <- subset(geochem.vec, geochem.vec$V4 <= 0.05 & geochem.vec$V3 >= 0.6)

#combine pigment and diatom vectors
diat.vec.sig$var <- "diatom"
pig.vec.sig$var <- "pigment"
comb.vec.sig <- rbind(diat.vec.sig,pig.vec.sig)

print(comb.vec.sig)


# Scale arrows by r2 to show correlation strength
comb.vec.sig.scaled <- comb.vec.sig %>%
  mutate(X1 = (Dim1*V3)/1.5,
         X2 = (Dim2*V3)/1.5)
# geochem.vec.sig.scaled <- geochem.vec.sig %>%
#   mutate(X1 = (X1*V3)/2,
#          X2 = (X2*V3)/2)


# Prepare the envfit output for diatoms and pigments
diat.vec <- envfit(comb.mds, env = master_v3_interp[,8:33], na.rm = TRUE)
pig.vec  <- envfit(comb.mds, env = master_v3_interp[,34:45], na.rm = TRUE)

# Combine vectors into a single data frame
diat.df <- data.frame(
  Dim1 = diat.vec$vectors$arrows[,1],
  Dim2 = diat.vec$vectors$arrows[,2],
  r    = diat.vec$vectors$r,
  pval = diat.vec$vectors$pvals,
  variable = rownames(diat.vec$vectors$arrows),
  var_type = "diatom"
)

pig.df <- data.frame(
  Dim1 = pig.vec$vectors$arrows[,1],
  Dim2 = pig.vec$vectors$arrows[,2],
  r    = pig.vec$vectors$r,
  pval = pig.vec$vectors$pvals,
  variable = rownames(pig.vec$vectors$arrows),
  var_type = "pigment"
)

vec.df <- rbind(diat.df, pig.df)

# Top loadings for PC1 and PC2 using absolute correlation (r) with the axis
top_loadings <- rbind(
  vec.df %>% arrange(desc(abs(Dim1))) %>% slice(1) %>% 
    mutate(axis = "PC1", total_variance = var1),
  vec.df %>% arrange(desc(abs(Dim2))) %>% slice(1) %>% 
    mutate(axis = "PC2", total_variance = var2)
)

top_loadings


# --- Envfit for diatoms and pigments ---
diat.vec <- envfit(comb.mds, env = master_v3_interp[,8:33], na.rm = TRUE)
pig.vec  <- envfit(comb.mds, env = master_v3_interp[,34:45], na.rm = TRUE)

# --- Combine vectors into one data frame ---
prepare_envfit_df <- function(envfit_obj, var_type) {
  data.frame(
    Dim1     = envfit_obj$vectors$arrows[,1],
    Dim2     = envfit_obj$vectors$arrows[,2],
    r        = envfit_obj$vectors$r,
    pval     = envfit_obj$vectors$pvals,
    variable = rownames(envfit_obj$vectors$arrows),
    var_type = var_type,
    stringsAsFactors = FALSE
  )
}

diat.df <- prepare_envfit_df(diat.vec, "diatom")
pig.df  <- prepare_envfit_df(pig.vec, "pigment")

vec.df <- rbind(diat.df, pig.df)

# --- Top loadings summary for each axis ---
# Store total variance for axes
total_var <- c(PC1 = var1, PC2 = var2)

# Function to get top loading for a given axis
get_top_loading <- function(df, axis_name, total_var) {
  axis_col <- ifelse(axis_name == "PC1", "Dim1", "Dim2")
  df %>%
    arrange(desc(abs(.data[[axis_col]]))) %>%
    slice(1) %>%
    mutate(axis = axis_name,
           total_variance = total_var[axis_name])
}

# Apply to all axes (PC1, PC2)
top_loadings <- bind_rows(
  get_top_loading(vec.df, "PC1", total_var),
  get_top_loading(vec.df, "PC2", total_var)
)

# View the summary table
top_loadings



write.csv(top_loadings, "/Users/a16512/Desktop/pcoa_top_loadings.csv", row.names = FALSE)

#plot PCoA paths for all lakes
ggplot(aes(X1, X2),data=site.comb) +  
  #lake trajectories 
  geom_path(aes(color=lake), size=0.75, 
            arrow=arrow(type="closed", ends="first",length=unit(0.1,"inches")), 
            lineend="round") + 
  #vectors
  geom_segment(data=comb.vec.sig.scaled,
               aes(x=0,y=0,xend=X1,yend=X2), 
               color="gray") + 
  geom_point(data=comb.vec.sig.scaled,
             aes(x=X1,y=X2,shape=var),color="gray",size=2)+
  geom_text(data=comb.vec.sig.scaled,
            aes(x=X1, y=X2),label=rownames(comb.vec.sig.scaled), size=4,
            position=position_jitter(width=0,height=0))+
  #general aesthetics
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Lake") +
  scale_x_continuous(limits=c(-0.75,0.75))+
  scale_y_continuous(limits=c(-0.75,0.75))+
  scale_color_manual(values=c("#d77932", #burnt
                              "#8b5a2b", #dunnigan
                              "#9a77b5", #elbow
                              "#6d91c1", #etwin
                              "#62497a", #finger
                              "#9b3d3d", #flame
                              "#e3b448", #smoke
                              "#3b5b7c" #wtwin
  )) +
  labs(
    x = paste0("PCoA Axis 1 (", var1, "%)"),
    y = paste0("PCoA Axis 2 (", var2, "%)"),
    color = "Lake"
  )+
  theme_classic(base_size=14) 

###calculating varience
eigs <- site.comb$eig
var1 <- round((eigs[1] / sum(eigs)) * 100, 1)
var2 <- round((eigs[2] / sum(eigs)) * 100, 1)

cat("Axis 1 =", var1, "%, Axis 2 =", var2, "% of total variance\n")


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
  diat.stand.hell.mix <- decostand(mix.master.dat[,8:33], method="hellinger")
  pig.stand.hell.mix <- decostand(mix.master.dat[,34:45], method="hellinger")
  #...to create distance matrices
  diat.dist.mix <- dist(diat.stand.hell.mix, method="euclidean")
  pig.dist.mix <- dist(pig.stand.hell.mix, method="euclidean")
  
  #combine distance matrices...
  comb.dist.mix <- diat.dist.mix+pig.dist.mix
  #...and scale
  comb.mds.mix <- cmdscale(d=comb.dist.mix, k=2, eig=TRUE)
  
  #extract coordinates
  comb.coords <- as.data.frame(comb.mds.mix$points)
  colnames(comb.coords) <- c("Dim1","Dim2")
  
  # Variance explained
  eig <- comb.mds.mix$eig
  var1 <- round((eig[1] / sum(eig)) * 100, 1)
  var2 <- round((eig[2] / sum(eig)) * 100, 1)
  
  comb.scores.mix <- comb.coords
  comb.scores.mix$lake <- mix.master.dat$lake
  comb.scores.mix$year <- mix.master.dat$year_loess
  comb.scores.mix$mix_regime <- mixing_regime[i]
  
  
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
  comb.scores.mix <- comb.scores.mix %>%
    mutate(color = case_when(
      lake == "flame"    ~ "#9b3d3d",
      lake == "burnt"    ~ "#d77932",
      lake == "smoke"    ~ "#e3b448",
      lake == "elbow"    ~ "#9a77b5",
      lake == "finger"   ~ "#62497a",
      lake == "etwin"    ~ "#6d91c1",
      lake == "wtwin"    ~ "#3b5b7c",
      lake == "dunnigan" ~ "#8b5a2b"
    ))
  
  
  
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
    labs(x = paste0("PCoA Axis 1 (", var1, "%)"),
         y = paste0("PCoA Axis 2 (", var2, "%)"), color = "Lake") +
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
  diat.stand.hell.lake <- decostand(lake.master.dat[,8:33], method="hellinger")
  pig.stand.hell.lake <- decostand(lake.master.dat[,34:45], method="hellinger")
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
  diat.vec.lake <- envfit(comb.mds.lake,env=lake.master.dat[,8:33],na.rm=T)
  diat.vec.lake <- data.frame(cbind(diat.vec.lake$vectors$arrows, 
                                   diat.vec.lake$vectors$r, 
                                   diat.vec.lake$vectors$pvals),
                             var_type="diatoms")
  #pigments
  pig.vec.lake <- envfit(comb.mds.lake,env=lake.master.dat[,34:45],na.rm=T)
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
  #colors
  comb.scores.mix <- comb.scores.mix %>%
    mutate(color = case_when(
      lake == "flame"    ~ "#9b3d3d",
      lake == "burnt"    ~ "#d77932",
      lake == "smoke"    ~ "#e3b448",
      lake == "elbow"    ~ "#9a77b5",
      lake == "finger"   ~ "#62497a",
      lake == "etwin"    ~ "#6d91c1",
      lake == "wtwin"    ~ "#3b5b7c",
      lake == "dunnigan" ~ "#8b5a2b"
    ))
  
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



pcoa.scores_lake_df %>%
  mutate(year_round = round_to_nearest_10(year)) %>% 
  group_by(lake, year_round) %>%
  summarize(Dim1 = median(Dim1), .groups = "drop") %>%
  ggplot(aes(x = year_round, y = Dim1, color = lake)) +
  geom_vline(xintercept = 1900, lty = "dashed") +
  geom_vline(xintercept = 1970, lty = "dashed") +
  geom_smooth(se = FALSE, size = 1.2) +
  theme_classic(base_size = 14) +
  scale_color_manual(
    values = c(
      "burnt"    = "#d77932",
      "dunnigan" = "#8b5a2b",
      "elbow"    = "#9a77b5",
      "etwin"    = "#6d91c1",
      "finger"   = "#62497a",
      "flame"    = "#9b3d3d",
      "smoke"    = "#e3b448",
      "wtwin"    = "#3b5b7c"
    )
  ) +
  labs(
    x = "Year",
    y = "PCoA Axis 1 (25.7%)",
    color = "Lake"
  )



#####Messing around with varience


# Define lake names and colors
lakes <- unique(master_v3_interp$lake)

lake_colors <- c(
  flame    = "#9b3d3d",
  burnt    = "#d77932",
  smoke    = "#e3b448",
  elbow    = "#9a77b5",
  finger   = "#62497a",
  etwin    = "#6d91c1",
  wtwin    = "#3b5b7c",
  dunnigan = "#8b5a2b"
)

# Empty lists
pcoa.plots <- list()
pcoa.scores_lake <- list()
pcoa.variable.vectors.lake <- list()

for (i in seq_along(lakes)) {
  
  # Subset data
  lake.master.dat <- subset(master_v3_interp, lake == lakes[i])
  
  # Standardize
  diat.stand.hell.lake <- decostand(lake.master.dat[,8:33], method="hellinger")
  pig.stand.hell.lake  <- decostand(lake.master.dat[,34:45], method="hellinger")
  
  # Distances
  diat.dist.lake <- dist(diat.stand.hell.lake)
  pig.dist.lake  <- dist(pig.stand.hell.lake)
  comb.dist.lake <- diat.dist.lake + pig.dist.lake
  
  # PCoA
  comb.mds.lake <- cmdscale(d=comb.dist.lake, k=2, eig=TRUE)
  
  # Extract variance
  eig <- comb.mds.lake$eig
  var1 <- round((eig[1] / sum(eig)) * 100, 1)
  var2 <- round((eig[2] / sum(eig)) * 100, 1)
  total_var <- round(((eig[1] + eig[2]) / sum(eig)) * 100, 1)
  
  # Combine with metadata
  comb.scores.lake <- as.data.frame(comb.mds.lake$points) %>%
    rename(Dim1 = V1, Dim2 = V2) %>%
    mutate(
      lake = lake.master.dat$lake,
      year = lake.master.dat$year_loess,
      decade = floor(year / 10) * 10,
      color = lake_colors[[lakes[i]]]
    )
  
  pcoa.scores_lake[[i]] <- comb.scores.lake
  
  # Fit vectors
  diat.vec.lake <- envfit(comb.mds.lake, lake.master.dat[,8:33], na.rm=TRUE)
  pig.vec.lake  <- envfit(comb.mds.lake, lake.master.dat[,34:45], na.rm=TRUE)
  
  # Extract & clean
  get_vec_df <- function(vec, var_type) {
    if (length(vec$vectors$arrows) == 0) return(data.frame())
    df <- data.frame(vec$vectors$arrows, r2 = vec$vectors$r, pval = vec$vectors$pvals)
    df <- subset(df, pval <= 0.05 & r2 >= 0.5)
    if (nrow(df) == 0) return(df)
    df <- df %>%
      mutate(Dim1 = (Dim1 * r2)/2, Dim2 = (Dim2 * r2)/2, var_type = var_type, lake = lakes[i],
             label = rownames(df))
    return(df)
  }
  
  diat.vec.lake.sig.scaled <- get_vec_df(diat.vec.lake, "diatom")
  pig.vec.lake.sig.scaled  <- get_vec_df(pig.vec.lake, "pigment")
  
  vec.lake <- rbind(diat.vec.lake.sig.scaled, pig.vec.lake.sig.scaled)
  pcoa.variable.vectors.lake[[i]] <- vec.lake
  
  # PCoA plot
  p <- ggplot(comb.scores.lake, aes(Dim1, Dim2)) +
    geom_path(color = lake_colors[[lakes[i]]], size = 0.8,
              arrow = arrow(type="closed", ends="first", length=unit(0.05,"inches")), lineend="round") +
    geom_point(color = lake_colors[[lakes[i]]], size = 2) +
    geom_text(aes(label = decade), color = "black", size = 3, vjust = -0.5, check_overlap = TRUE) +
    
    # Vectors
    geom_segment(data = pig.vec.lake.sig.scaled,
                 aes(x=0, y=0, xend=Dim1, yend=Dim2),
                 color="darkgreen", arrow=arrow(length=unit(0.1,"inches"))) +
    geom_text(data = pig.vec.lake.sig.scaled,
              aes(x=Dim1, y=Dim2, label=label),
              color="darkgreen", size=3, vjust=-0.5) +
    
    geom_segment(data = diat.vec.lake.sig.scaled,
                 aes(x=0, y=0, xend=Dim1, yend=Dim2),
                 color="orange", arrow=arrow(length=unit(0.1,"inches"))) +
    geom_text(data = diat.vec.lake.sig.scaled,
              aes(x=Dim1, y=Dim2, label=label),
              color="orange", size=3, vjust=-0.5) +
    
    labs(
      x = paste0("PCoA Axis 1 (", var1, "%)"),
      y = paste0("PCoA Axis 2 (", var2, "%)")
    ) +
    theme_classic(base_size = 13)
  
  pcoa.plots[[i]] <- p
}

p <- p + theme(
  plot.margin = margin(t = 20, r = 15, b = 15, l = 15)
)

# Combine plots
cowplot::plot_grid(
  plotlist = pcoa.plots,
  ncol = 2,
  label_fontface = "plain",
  rel_widths = c(1, 1),   # increases spacing horizontally
  rel_heights = c(1, 1)   # increases spacing vertically
)


#################Calculating the highest loading
####Across seperate lakes
# Create empty lists to store outputs
pcoa.plots <- list()
pcoa.scores_lake <- list()
pcoa.variable.vectors.lake <- list()
pcoa.top_loadings <- list()  # new list to hold top 3 loadings per lake

for (i in seq_along(lakes)) {
  
  # Subset data for the lake
  lake.master.dat <- subset(master_v3_interp, lake == lakes[i])
  
  # Standardize diatoms and pigments
  diat.stand.hell.lake <- decostand(lake.master.dat[,8:33], method="hellinger")
  pig.stand.hell.lake  <- decostand(lake.master.dat[,34:45], method="hellinger")
  
  # Create distance matrices and combine
  diat.dist.lake <- dist(diat.stand.hell.lake, method="euclidean")
  pig.dist.lake  <- dist(pig.stand.hell.lake,  method="euclidean")
  comb.dist.lake <- diat.dist.lake + pig.dist.lake
  
  # Perform PCoA with eigenvalues
  comb.mds.lake <- cmdscale(d=comb.dist.lake, k=2, eig=TRUE)
  
  # Extract coordinates and variance
  comb.coords <- as.data.frame(comb.mds.lake$points)
  colnames(comb.coords) <- c("Dim1","Dim2")
  eig <- comb.mds.lake$eig
  var1 <- round((eig[1] / sum(eig)) * 100, 1)
  var2 <- round((eig[2] / sum(eig)) * 100, 1)
  total_var <- round(((eig[1] + eig[2]) / sum(eig)) * 100, 1)
  
  # Add sample info and decade labels
  comb.scores.lake <- comb.coords %>%
    mutate(
      lake = lake.master.dat$lake,
      year = lake.master.dat$year_loess,
      decade = floor(year / 10) * 10,  # decade labels
      color = lake_colors[lake]
    )
  
  pcoa.scores_lake[[i]] <- comb.scores.lake
  
  # Fit environmental vectors
  diat.vec.lake <- envfit(comb.mds.lake, env=lake.master.dat[,8:33], na.rm=TRUE)
  pig.vec.lake  <- envfit(comb.mds.lake, env=lake.master.dat[,34:45], na.rm=TRUE)
  
  # Extract vector data
  diat.vec.lake <- data.frame(cbind(diat.vec.lake$vectors$arrows, diat.vec.lake$vectors$r, diat.vec.lake$vectors$pvals))
  pig.vec.lake  <- data.frame(cbind(pig.vec.lake$vectors$arrows,  pig.vec.lake$vectors$r,  pig.vec.lake$vectors$pvals))
  colnames(diat.vec.lake) <- colnames(pig.vec.lake) <- c("Dim1","Dim2","r2","pval")
  
  # Keep only significant vectors
  diat.vec.lake.sig <- subset(diat.vec.lake, pval <= 0.05 & r2 >= 0.5)
  pig.vec.lake.sig  <- subset(pig.vec.lake,  pval <= 0.05 & r2 >= 0.5)
  
  # Scale vectors
  diat.vec.lake.sig.scaled <- diat.vec.lake.sig %>%
    mutate(Dim1 = (Dim1 * r2)/2, Dim2 = (Dim2 * r2)/2)
  pig.vec.lake.sig.scaled <- pig.vec.lake.sig %>%
    mutate(Dim1 = (Dim1 * r2)/2, Dim2 = (Dim2 * r2)/2)
  
  # Combine vectors for saving
  vec.lake <- rbind(
    mutate(diat.vec.lake.sig.scaled, var_type="diatom", lake=lakes[i]),
    mutate(pig.vec.lake.sig.scaled,  var_type="pigment", lake=lakes[i])
  )
  pcoa.variable.vectors.lake[[i]] <- vec.lake
  
  # --- 🔍 NEW SECTION: Calculate highest loadings on Axis 1 and Axis 2 ---
  
  get_top_loadings <- function(df, var_group) {
    if (nrow(df) == 0) return(NULL)
    df <- df %>%
      mutate(
        Axis1_corr = abs(Dim1),
        Axis2_corr = abs(Dim2),
        variable = rownames(df),
        var_group = var_group
      )
    top1 <- df %>% arrange(desc(Axis1_corr)) %>% head(3)
    top2 <- df %>% arrange(desc(Axis2_corr)) %>% head(3)
    bind_rows(
      mutate(top1, axis = "PCoA1"),
      mutate(top2, axis = "PCoA2")
    )
  }
  
  top_diat <- get_top_loadings(diat.vec.lake, "diatom")
  top_pig  <- get_top_loadings(pig.vec.lake,  "pigment")
  
  top_all <- bind_rows(top_diat, top_pig) %>%
    mutate(lake = lakes[i], total_variance = total_var)
  
  pcoa.top_loadings[[i]] <- top_all
  
  # --- 🔚 Plot section ---
  
  p <- ggplot(comb.scores.lake, aes(Dim1, Dim2)) +
    geom_path(color=lake_colors[lakes[i]], size=0.8,
              arrow=arrow(type="closed", ends="first", length=unit(0.05,"inches")), lineend="round") +
    geom_point(color=lake_colors[lakes[i]], size=2) +
    geom_text(aes(label=decade), size=3, vjust=-0.5, color="black") +
    
    geom_segment(data=pig.vec.lake.sig.scaled, aes(x=0,y=0,xend=Dim1,yend=Dim2), 
                 color="darkgreen", arrow=arrow(length=unit(0.1,"inches"))) +
    geom_text(data=pig.vec.lake.sig.scaled, aes(x=Dim1, y=Dim2, label=rownames(pig.vec.lake.sig.scaled)), 
              color="darkgreen", size=3, vjust=-0.5) +
    
    geom_segment(data=diat.vec.lake.sig.scaled, aes(x=0,y=0,xend=Dim1,yend=Dim2), 
                 color="orange", arrow=arrow(length=unit(0.1,"inches"))) +
    geom_text(data=diat.vec.lake.sig.scaled, aes(x=Dim1, y=Dim2, label=rownames(diat.vec.lake.sig.scaled)), 
              color="orange", size=3, vjust=-0.5) +
    
    labs(
      x = paste0("PCoA Axis 1 (", var1, "%)"),
      y = paste0("PCoA Axis 2 (", var2, "%)"),
      title = paste0(lake_plot_titles[i], " — Total variance: ", total_var, "%")
    ) +
    theme_classic(base_size=13)
  
  pcoa.plots[[i]] <- p
}

# Combine loadings into one table for easy viewing
pcoa.top_loadings_df <- bind_rows(pcoa.top_loadings)

# View top variables for each lake
print(pcoa.top_loadings_df)

# Export if desired
write.csv(pcoa.top_loadings_df, "pcoa_top_loadings.csv", row.names = FALSE)
getwd()

write.csv(pcoa.top_loadings_df, "/Users/a16512/Desktop/pcoa_top_loadings.csv", row.names = FALSE)

########Across mixing regimes

# Add mixing regime
master_v3_interp <- master_v3_interp %>%
  mutate(mix.regime = case_when(
    lake %in% c("burnt","etwin","elbow","smoke") ~ "Polymictic",
    lake %in% c("dunnigan","finger") ~ "Well mixed",
    lake %in% c("wtwin","flame") ~ "Dimictic"
  ))

# Setup
mixing_regime <- unique(master_v3_interp$mix.regime)
pcoa.scores_mixing.regime <- list()
pcoa.variable.vectors <- list()
pcoa.plots <- list()
top_loadings <- data.frame()  # store top loading info for all regimes

# Loop over mixing regimes
for (i in seq_along(mixing_regime)) {
  
  # Subset data
  mix.master.dat <- subset(master_v3_interp, mix.regime == mixing_regime[i])
  
  # Standardize diatoms and pigments
  diat.stand <- decostand(mix.master.dat[,8:33], method = "hellinger")
  pig.stand  <- decostand(mix.master.dat[,34:45], method = "hellinger")
  
  # Distance matrices
  diat.dist <- dist(diat.stand, method = "euclidean")
  pig.dist  <- dist(pig.stand, method = "euclidean")
  
  # Combined distance
  comb.dist <- diat.dist + pig.dist
  
  # PCoA
  comb.mds <- cmdscale(d = comb.dist, k = 2, eig = TRUE)
  comb.coords <- as.data.frame(comb.mds$points)
  colnames(comb.coords) <- c("Dim1","Dim2")
  
  # Variance explained
  eig <- comb.mds$eig
  var1 <- round((eig[1] / sum(eig)) * 100, 1)
  var2 <- round((eig[2] / sum(eig)) * 100, 1)
  total_var <- round(((eig[1]+eig[2])/sum(eig))*100, 1)
  
  # Add sample info
  comb.coords <- comb.coords %>%
    mutate(lake = mix.master.dat$lake,
           year = mix.master.dat$year_loess,
           mix_regime = mixing_regime[i])
  
  pcoa.scores_mixing.regime[[mixing_regime[i]]] <- comb.coords
  
  # Envfit vectors
  diat.vec <- envfit(comb.mds, env = mix.master.dat[,8:33], na.rm = TRUE)
  pig.vec  <- envfit(comb.mds, env = mix.master.dat[,34:45], na.rm = TRUE)
  
  # Make data.frames
  diat.vec.df <- data.frame(diat.vec$vectors$arrows,
                            r = diat.vec$vectors$r,
                            pval = diat.vec$vectors$pvals,
                            var_type = "diatom")
  rownames(diat.vec.df) <- rownames(diat.vec$vectors$arrows)
  
  pig.vec.df <- data.frame(pig.vec$vectors$arrows,
                           r = pig.vec$vectors$r,
                           pval = pig.vec$vectors$pvals,
                           var_type = "pigment")
  rownames(pig.vec.df) <- rownames(pig.vec$vectors$arrows)
  
  # Combine
  vec.df <- rbind(diat.vec.df, pig.vec.df) %>%
    mutate(variable = rownames(.),
           mix_regime = mixing_regime[i])
  pcoa.variable.vectors[[mixing_regime[i]]] <- vec.df
  
  # Top loadings
  pc1_top <- vec.df %>% arrange(desc(abs(Dim1))) %>% slice(1) %>%
    mutate(axis = "PC1", var_group = var_type, total_variance = total_var)
  pc2_top <- vec.df %>% arrange(desc(abs(Dim2))) %>% slice(1) %>%
    mutate(axis = "PC2", var_group = var_type, total_variance = total_var)
  
  # Combine and add to summary
  top_loadings <- rbind(top_loadings, pc1_top, pc2_top)
  
  # Significant vectors for plotting
  sig.vec <- vec.df %>% filter(r >= 0.6, pval <= 0.05)
  
  # Scale vectors for plotting
  sig.vec <- sig.vec %>%
    mutate(Dim1 = Dim1 * r / 2,
           Dim2 = Dim2 * r / 2)
  
  # Colors
  comb.coords <- comb.coords %>%
    mutate(color = case_when(
      lake == "flame" ~ "#9b3d3d",
      lake == "burnt" ~ "#d77932",
      lake == "smoke" ~ "#e3b448",
      lake == "elbow" ~ "#9a77b5",
      lake == "finger" ~ "#62497a",
      lake == "etwin" ~ "#6d91c1",
      lake == "wtwin" ~ "#3b5b7c",
      lake == "dunnigan" ~ "#8b5a2b"
    ))
  
  # Plot
  p <- ggplot(comb.coords, aes(Dim1, Dim2)) +
    geom_path(aes(color = lake), size = 0.75,
              arrow = arrow(type = "closed", ends = "first", length = unit(0.05, "inches")),
              lineend = "round") +
    geom_point(aes(color = lake), size = 2) +
    geom_segment(data = sig.vec, aes(x = 0, y = 0, xend = Dim1, yend = Dim2, color = var_type),
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_text(data = sig.vec, aes(x = Dim1*1.2, y = Dim2*1.2, label = variable, color = var_type),
              size = 3, vjust = -0.5, position = position_jitter(width = 0.1, height = 0.1)) +
    labs(x = paste0("PCoA Axis 1 (", var1, "%)"),
         y = paste0("PCoA Axis 2 (", var2, "%)"),
         title = paste0(mixing_regime[i], " — Total variance: ", total_var, "%"),
         color = "Variable type / Lake") +
    scale_color_manual(values = c("diatom" = "orange", "pigment" = "darkgreen",
                                  unique(comb.coords$color))) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  
  pcoa.plots[[mixing_regime[i]]] <- p
}

# Inspect top loadings for all regimes
top_loadings


# Combine top loadings summary into a single data.frame
top_loadings_df <- do.call(rbind, top_loadings_summary)

# Optional: View the summary
print(top_loadings_df)

write.csv(top_loadings, "/Users/a16512/Desktop/pcoa_top_loadings.csv", row.names = FALSE)
