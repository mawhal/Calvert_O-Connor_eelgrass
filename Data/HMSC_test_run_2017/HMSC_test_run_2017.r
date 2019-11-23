## Test HMSC run at Seagrass Retreat
# started 23 November 2019

# load libraries
library(Hmsc)
library(tidyverse)

### read data
# response data
<<<<<<< HEAD
setwd("~/Github/Calvert_O-Connor_eelgrass/Data/HMSC_test_run_2017")

Ygrazer <- read.csv("Hakai_2017_mesograzer_comm.csv")





Ymicrobe <- read.csv("microbes_2017.csv")

Yfish

# explanatory data
X <- read.csv("Hakai_2017_Eelgrass_Biometrics_final.csv")
names(X)[3] <- "sample"

meta <- read.csv("Hakai_2017_mesograzer_meta.csv")
X1 <- X[,-c(4:50)]

X2 <- X1 %>%
  select(-detritus_wet_mass, -epi_algae_wet_mass, -seagrass_wet_mass, -avg_shoot_surface_area_m2, -shoots_m.2, -notes, -date)

Xab <- read.csv("biooracle_2017_small.csv")
Xab$site <- gsub("_"," ",Xab$site)

X3 <- left_join(X2, Xab)

X4<- left_join(meta,X3)

Xgr <- X4[,-c(1:3)]


X5 <-left_join(Ymicrobe[,1:4], X3)

Xmi <- X5[,-c(1:4)]

# Random effects
studyDesign <- read.csv("pi.csv")

#spatial data
spatial <- read.csv("spatial_2017.csv")



Xgr[is.na(Xgr)] <- 0
Y[is.na(Y)] <- 0
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
#rL3= HmscRandomLevel(sData = spatial)

### building the model

## formula for fixed effects
XFormula = ~()

mgrazer <- Hmsc( Y = Ygrazer, 
           XData = Xgr, XFormula = ~detritus_dry_mass + epi_algae_dry_mass + seagrass_dry_mass + lai + microepiphyte_wt_mg + nitraterange + tempmean + temprange + curveloc + nitratemean + salinityrange + salinitymean + dissoxmean + dissoxrange, distr = "poisson", studyDesign = studyDesign )


### running the model
thin = 1
samples = 100
nChains = 2
set.seed(1)

# Run MCMC chains. took at least 12 hours on KS laptop

mod <- sampleMcmc(mgrazer, samples = 100 , transient = 50, thin = 1, verbose = 10, nChains = 2)
