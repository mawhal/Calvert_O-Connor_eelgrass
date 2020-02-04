### 2015 HMSC analysis on core bacteria ###
### Author: Bianca Trevizan Segovia, original HMSC code: Keila Stark ###
### Date created: January 07, 2020 ###

# and load it
library(Hmsc)
library(tidyverse)

### read all data
all_microbes_factors_2015 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/16S_retreat/hmsc_16S/2015/microbes_factors_2015.csv")
names(all_microbes_factors_2015 )

# response data
Ymicrobes_2015 <- all_microbes_factors_2015 %>% 
  dplyr::select(13:32)
names(Ymicrobes_2015)

# explanatory data = abiotic + seagrass metrics
Xenvironmental_2015 <- all_microbes_factors_2015 %>% 
  dplyr::select(7:12)
names(Xenvironmental_2015)

#factors representing random effects (quadrat, site, region level)
studyDesign <- all_microbes_factors_2015 %>% 
  dplyr::select(2:4)
names(studyDesign)

# spatial data (coordinates)
spatial <- all_microbes_factors_2015 %>% 
  dplyr::select(5,6)
spatial <- as.data.frame(spatial)
spatial$longitude <- spatial$longitude-min(spatial$longitude)
spatial <- unique(spatial)
rownames(spatial) <- unique(studyDesign$Site)

## Remove NA's / get into right format
Xenvironmental_2015[is.na(Xenvironmental_2015)] <- 0
Ymicrobes_2015[is.na(Ymicrobes_2015)] <- 0
Ymicrobes_2015 <- as.matrix(Ymicrobes_2015)
Xenvironmental_2015 <- as.data.frame(Xenvironmental_2015)

rL1 = HmscRandomLevel(units = studyDesign$Quadrat) #quadrat level
rL2 = HmscRandomLevel(units = studyDesign$Site)
rL3 = HmscRandomLevel(units = studyDesign$Region)
rL4= HmscRandomLevel(sData = spatial) # lat&long coordinates

### building the model

names(Xenvironmental_2015)
## formula for fixed effects
XFormula = ~(temperature + dissolved_oxygen_concentration + salinity + 
               lai_m2 + microepiphyte.y + biomass_dry_wt)

Ymicrobes_2015 <- as.matrix(Ymicrobes_2015)

mmicrobes <- Hmsc( Y = Ymicrobes_2015, 
           XData = Xenvironmental_2015, XFormula = XFormula, 
           distr = "poisson", studyDesign = studyDesign,
           ranLevels = list(Quadrat=rL1, 
                            #Site=rL2,
                            Region=rL3,
                            Site=rL4)) # Here we have to specify the structure of the data, which means that the rL4 spatial component has to contain the row names for it in the studyDesign (that is why in the other code, we had to run spat, to include the rownames for spatial there:  
#spat<- data.frame(spat = sprintf('spatial_%.2d',1:36)) #spatial factor column for studyDesign
#studyDesign <- cbind(studyDesign, spat)
# But because our lat long correspond to our sites, I could just use the Sites from the studyDesign instead, and don't really need to use it as random

### running the model
thin = 1
samples = 100
nChains = 2 # have to compute multiple chains to check if they all yield similar results
set.seed(1)

# Run MCMC chains. took at least 12 hours on KS laptop

mod <- sampleMcmc(mmicrobes, samples = samples , transient = 50, 
                  thin = thin, verbose = 10, nChains = nChains) # nParallel = 3 when running more than 1 Chain for the real model = running more than one chain allows us to check if all arer doing the same thing
# keila's parameters thin = 100 niter = 200000 nburn = 100000

### model results

## Assess model fit
preds = computePredictedValues(mod)
MF <- evaluateModelFit(hM = mod, predY = preds)

## variance partitioning
VP = computeVariancePartitioning(mod) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
plotVariancePartitioning(mod, VP = VP)
VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c("temperature","dissolved_oxygen_concentration",
                           "salinity", "lai_m2","microepiphyte.y",
                           "biomass_dry_wt","quadrat", "region", "spatial" ), 
                         levels = rev(c("temperature",
                                        "dissolved_oxygen_concentration",
                                        "salinity", "lai_m2","microepiphyte.y",
                                        "biomass_dry_wt","quadrat", 
                                        "region", "spatial")), 
                         ordered = TRUE)) %>% 
  gather(key = taxon, value = variance, -effect) %>% 
  group_by(taxon) 


VP.df$taxon <- factor(VP.df$taxon, 
                      levels = colnames(mod$Y)[order(colSums(mod$Y),decreasing = TRUE)], 
                      ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), taxon = colnames(mod$Y))

microbes <- ggplot(VP.df,aes(y = variance, x = taxon, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","magenta", "darkolivegreen","darkolivegreen1", "forestgreen", "khaki", "yellow","gold"))+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")
microbes

ggsave("~/PostDoc/projects/Hakai_Quadra_data_retreat/16S_retreat/hmsc_16S/2015/hmsc_microbes_2015_test.jpeg", plot = microbes,width=250, height=200, units="mm",dpi=300)

##### NOW EVALUATING #####
# To look at the parameter estimates, we extract the posterior distribution from the model object and convert it into a coda object.
mpost = convertToCodaObject(mod) 
summary(mpost$Beta)

# Checking MCMC diagnostics
## Let us now return to the parameters that guide posterior sampling in the MCMC algorithm. We first plot the trace plots of the β-parameters.
plot(mpost$Beta)
# par(mar=c(1,1,1,1)) if plot gives error, run this


