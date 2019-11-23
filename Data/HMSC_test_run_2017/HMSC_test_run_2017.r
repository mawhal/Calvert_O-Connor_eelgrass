## Test HMSC run at Seagrass Retreat
# started 23 November 2019

# load libraries
library(Hmsc)
library(tidyverse)

### read data
# response data
# setwd("~/Github/Calvert_O-Connor_eelgrass/Data/HMSC_test_run_2017")

Ygrazer <- read.csv("Hakai_2017_mesograzer_comm.csv")
meta <- read.csv("Hakai_2017_mesograzer_meta.csv")
# filter out sites we want
meta.comm <- cbind( meta, Ygrazer )
meta.comm <- meta.comm %>% 
  filter( site != "pruth_pocket" )
# separate again
meta <- meta.comm[,c(1:3)]
comm <- meta.comm[,-c(1:3)]
Ygrazer <- comm

# get rid of some taxa
comm <- comm[, -c(1,12)]

Ymicrobe <- read.csv("microbes_2017.csv")

Yfish

# explanatory data
X <- read.csv("Hakai_2017_Eelgrass_Biometrics_final.csv")
names(X)[3] <- "sample"

X1 <- X[,-c(4:50)]

X2 <- X1 %>%
  select(-detritus_wet_mass, -epi_algae_wet_mass, -seagrass_wet_mass, -avg_shoot_surface_area_m2, -shoots_m.2, -notes, -date)

Xab <- read.csv("biooracle_2017_small.csv")
# Xab$site <- gsub("_"," ",Xab$site)

X3 <- left_join(X2, Xab)

X4<- left_join(meta,X3)

Xgr <- X4[,-c(1:3)]


X5 <-left_join(Ymicrobe[,1:4], X3)

Xmi <- X5[,-c(1:4)]

# Random effects
studyDesign <- read.csv("pi.csv")
row.names(Xgr) <- studyDesign$Quadrat

#spatial data
spatial <- read.csv("spatial_2017.csv")



Xgr[is.na(Xgr)] <- 0
# Ygrazer[is.na(Y)] <- 0
Ygrazer <- as.matrix(Ygrazer)
Xgr <- as.data.frame(Xgr)
Xmi <- as.data.frame(Xmi)
#spatial <- as.data.frame(spatial)
#spatial$longitude <- spatial$longitude-min(spatial$longitude)

#spat <- data.frame(spat = sprintf('spatial_%.2d',1:78)) #spatial factor column for studyDesign
#studyDesign <- cbind(studyDesign, spat)
#studyDesign$Spatial <- factor(studyDesign$Quadrat)

rL1 = HmscRandomLevel(units = studyDesign$Quadrat)
rL2 = HmscRandomLevel(units = studyDesign$Site)
rL3 = HmscRandomLevel(units = studyDesign$Region)
#rL3= HmscRandomLevel(sData = spatial)

### building the model

## formula for fixed effects
XFormula = ~(detritus_dry_mass + epi_algae_dry_mass + seagrass_dry_mass + 
               lai + microepiphyte_wt_mg + nitraterange + tempmean + 
               temprange + curveloc + nitratemean + salinityrange + 
               salinitymean + dissoxmean + dissoxrange)

mgrazer <- Hmsc( Y = Ygrazer, 
           XData = Xgr, XFormula = XFormula, 
           distr = "poisson", studyDesign = studyDesign,
           ranLevels = list(Quadrat=rL1, 
                            Site=rL2,
                            Region=rL3))


### running the model
thin = 1
samples = 100
nChains = 1
set.seed(1)

# Run MCMC chains. took at least 12 hours on KS laptop

mod <- sampleMcmc(mgrazer, samples = samples , transient = 50, 
                  thin = 1, verbose = 10, nChains = nChains)




### model results

## variance partitioning
VP = computeVariancePartitioning(mod) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
plotVariancePartitioning(mod, VP = VP)

VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c("detritus","macroepiph","eel_biomass",
                           "LAI","microepiph",
                           "site","transect"), 
                         levels = rev(c("elevation","elev.square","year",
                                        "elev:year","elev2:year",
                                        "site","transect")), 
                         ordered = TRUE)) %>% 
  # mutate(effect = factor(c("elevation","elev.square","temp.anom.sum", "temp.anom.win",
  #                          "elev:year","elev2:year",
  #                          "site","transect","ty"), 
  #                        levels = rev(c("elevation","elev.square","temp.anom.sum", "temp.anom.win",
  #                                       "elev:year","elev2:year",
  #                                       "site","transect","ty")), 
  #                        ordered = TRUE)) %>% 
  # mutate(effect = factor(c("elevation","elev.square","temp.anom.sum", "temp.anom.win","transect","year","site"), 
  #                        levels = rev(c("elevation","elev.square","temp.anom.sum", "temp.anom.win","transect","year","site")), 
  #                        ordered = TRUE)) %>% 
  gather(key = species, value = variance, -effect) %>% 
  group_by(species) %>% 
  mutate(tempR2 = variance[effect == "year"])
mutate(tempR2 = variance[effect == "temp.anom.sum"])

hold <- VP.df %>% filter(effect == "temp.anomaly") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(m$Y))

windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon", viridis(7)), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")




