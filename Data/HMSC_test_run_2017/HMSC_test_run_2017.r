## Test HMSC run at Seagrass Retreat
# started 23 November 2019

# load libraries
library(Hmsc)
library(tidyverse)

### read data
# response data
setwd("~/Github/Calvert_O-Connor_eelgrass/Data/HMSC_test_run_2017")

Ygrazer <- read.csv("Hakai_2017_mesograzer_comm.csv")


Ymicrobe <- read.csv("metadata_microbes_2017 (2).csv")

Yfish

# explanatory data
X <- read.csv("Hakai_2017_Eelgrass_Biometrics_final.csv")
X1 <- X[,-c(4:50)]

X2 <- X1 %>%
  select(-detritus_wet_mass, -epi_algae_wet_mass, -seagrass_wet_mass, -avg_shoot_surface_area_m2, -shoots_m.2, -X, -X.1, -X.2, -X.3, -X.4, -X.5, -X.6,  -notes, -date)

Xab <- read.csv("biooracle_2017.csv")

#pi <- read.csv("pi.csv")

spatial <- read.csv("spatial_2017.csv")

# filter rows for each response

# random effects
# quad, site, space? Do we need quad and site if we have space?



### building the model

## formula for fixed effects
XFormula = ~()

mgrazer <- Hmsc( Y = Ygrazer, 
           XData = Xgrazer, XFormula = XFormula,
           distr = "poisson",
           studyDesign = studyDesign, 
           ranLevels = list( space=rLspace ) )


### running the model
thin = 1
samples = 100
nChains = 2
set.seed(1)


