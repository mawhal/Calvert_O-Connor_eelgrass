####################################################################
###  Seagrass Microbial Data From Calvert Island 
###  Dataset begins in 2015
###  Samples from several sites (not always the same)
###
###  This code will run HMSC on the input microbial dataset
###  code by Bia Trevizan Segovia
###  started on   07 February 2020
###  
####################################################################

# load libraries
library(Hmsc)
library(tidyverse)
library(viridis)
library(corrplot)
library(vegan) 
  
### read data
# read data
microbes <- read.csv( "~/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/data_parfrey/16S/genus_16S_core_microbes.csv" )
names(microbes)

######################################################
######################## 2015 ########################
######################################################

## pick a year
microbes_2015 <- microbes %>% 
  filter( year==2015 )

microbes_2015 <- microbes_2015 %>%  
  mutate(quadrat_id = as.numeric(quadrat_id))

# explanatory data
environmental_coord <- read.csv("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R Code and Analysis/output data/merged_explanatory_20200214.csv")

## pick a year
keep_2015 <- c("2015")
environmental_2015 <- environmental_coord %>% 
  dplyr::filter(year %in% keep_2015)

### microbial sites for 2015 are: choked_inner, choked_sandspit, goose_south_east, goose_south_west, mcmullins_north, mcmullins_south, triquet_north, triquet_south
keep_sites <- c("choked_inner", "choked_sandspit", "goose_south_east", "goose_south_west", "mcmullins_north", "mcmullins_south", "triquet_north", "triquet_south")
environmental_2015_sites <- environmental_2015 %>% 
  dplyr::filter(site %in% keep_sites)

environmental_2015_sites <- environmental_2015_sites %>%  
  mutate(quadrat_id = as.numeric(quadrat_id))

### merge microbes and env_coord
microbes_metadata <- inner_join(environmental_2015_sites, microbes_2015, by = c("region","site", "quadrat_id"))
View(microbes_metadata)
# have to exclude mcmullins_north_5 because quadrat_biomass is NA
quadrat_with_NA <- c("mcmullins_north_5")
microbes_metadata <- microbes_metadata %>% 
  dplyr::filter(!quadrat %in% quadrat_with_NA)

# response data
names(microbes_metadata)
Ymicrobes <- microbes_metadata %>% 
  dplyr::select(27:61)
names(Ymicrobes)

## select columns for fixed and random effects
# explanatory data = abiotic + seagrass metrics
Xenvironmental <- microbes_metadata %>% 
  dplyr::select(temperature,salinity,depth,bed_area_m2, quadrat_biomass_g,quadrat_lai,quadrat_microepiphyte_mg,quadrat_macroalgae_g) #do,

# random effects (quadrat, site, region level)
studyDesign <- microbes_metadata %>% 
  dplyr::select(region, site, quadrat)
  
studyDesign <-  as.data.frame(studyDesign)
studyDesign$region <- factor(studyDesign$region)
studyDesign$site <- factor(studyDesign$site)
studyDesign$quadrat <- factor(studyDesign$quadrat)

# spatial data from lats and longs
# make longitude positive
spatial <- microbes_metadata %>% 
  dplyr::select(lat, long)
spatial$long <- spatial$long-min(spatial$long)
spatial <- as.data.frame(spatial)
spatial <- unique(spatial)
rownames(spatial) <- unique(studyDesign$site)

# formats
Ymicrobes <- as.matrix(Ymicrobes)
Xenvironmental <- as.data.frame(Xenvironmental)
#View(Xenvironmental)

# random effects structure
rL1 = HmscRandomLevel( units = unique(studyDesign$quadrat) ) #quadrat level
rL2 = HmscRandomLevel( units = unique(studyDesign$site) )
rL3 = HmscRandomLevel( units = unique(studyDesign$region) )
rL4 = HmscRandomLevel( sData = spatial )


## formula for fixed effects
XFormula = ~( temperature+salinity+depth+bed_area_m2+ quadrat_biomass_g+quadrat_lai+quadrat_microepiphyte_mg+quadrat_macroalgae_g)

# the model
mmicrobe_2015 <- Hmsc( Y = Ymicrobes, 
                   XData = Xenvironmental, XFormula = XFormula, 
                   distr = "poisson", studyDesign = studyDesign,
                   ranLevels = list(quadrat=rL1, 
                                    #site=rL2,
                                    region=rL3,
                                    site=rL4) )

### running the model
thin = 100
samples = 100
nChains = 2 # have to compute multiple chains to check if they all yield similar results
set.seed(1)

# Run MCMC chains. took at least 12 hours on KS laptop

mod <- sampleMcmc(mmicrobe_2015, samples = 2000 , transient = 100, 
                  thin = 100, verbose = 10, nChains = 2) # nParallel = 3 when running more than 1 Chain for the real model = running more than one chain allows us to check if all arer doing the same thing
# keila's = sampleMcmc(hM, samples = 2000, transient = 100, thin = 100, verbose = 10) (edited)
# saving
model <- "microbes_16S_2015"
ModelDir <- paste0( here::here(), "/R Code and Analysis/output data")
filename = file.path(ModelDir, paste("model_",as.character(model),
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples),".Rdata",sep = ""))
save(mod,file=filename)


load(filename)


## model results
m <- mod

# We evaluate MCMC convergence in terms of four kinds of parameters that we are especially interested in: the species niches `Beta`, the influence of traits on species niches `Gamma`, residual species associations `Omega`, and the strength of phylogenetic signal `rho`. As there are `r ns` species, the matrix `Omega` is a `r ns` x `r ns` matrix with `r ns*ns` elements. To avoid excessive computational times, we evaluate MCMC convergence for `Omega` only for a subsample of 100 randomly selected species pairs.
mpost = convertToCodaObject(m)



## Assess model fit
preds = computePredictedValues(m)
MF <- evaluateModelFit(hM = m, predY = preds)


## parameter estimates
postBeta = getPostEstimate(m, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature", "salinity", "quadrat_biomass_g", 
                                "quadrat_lai", "quadrat_microepiphyte_mg" ),
                            levels = c("intercept", "temperature", "do", "salinity", "quadrat_biomass_g", 
                                       "quadrat_lai", "quadrat_microepiphyte_mg" ), 
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 7), 
                          levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)],
                          ordered = TRUE)

windows(12,4)
ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")

# pull taxa with positive estimated response to temperature anomaly
taxa.pos.anom <- as.character(pos.neg[ pos.neg$parameter=='quadrat_lai' & pos.neg$value>0, 'species'])

# calculate mean and variance of parameter esimates
pbdf <- data.frame( t(postBeta$mean), taxon=colnames(postBeta$mean) )
names(pbdf) <- c("intercept","temp","do","sal",
                 "biomass","LAI","microepi","taxon")
# ## Add some basic trait information
# pbdf$alga <- "alga"
# pbdf$alga[c(3,15,20,56,57,62,68,82)] <- "invert"
# # More specific groups


coef_plot <- pbdf %>% 
  select( -taxon ) %>%
  gather( coefficient, value )

windows(6,4)
ggplot( coef_plot, aes(y=value, x=factor(1))) +
  facet_wrap(~coefficient, scales="free_y") +
  geom_hline(yintercept=0) +
  stat_summary( fun.data=mean_cl_boot, geom='errorbar',
                size=1, width=0.1, position=position_dodge(width=0.9) )  +
  stat_summary( fun.y=mean, geom='point',
                size=3, position=position_dodge(width=0.9) )  +
  ylab('coefficient') + xlab('') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

## variance partitioning
VP = computeVariancePartitioning(m) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP)

VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c("temp","do","sal",
                            "biomass","LAI","microepi",
                           "quadrat","region","site"), 
                         levels = rev(c("temp","do","sal",
                                        "biomass","LAI","microepi",
                                        "quadrat","region","site")), 
                         ordered = TRUE)) %>% 
  gather(key = species, value = variance, -effect) %>% 
  group_by(species) %>% 
  mutate(tempR2 = variance[effect == "temp"])


hold <- VP.df %>% filter(effect == "temp") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(m$Y))

windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkgreen","forestgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")

## associations
OmegaCor = computeAssociations(m)
supportLevel = 0.95
# choose the random variable to plot
rlevel = 3
toPlot = ((OmegaCor[[rlevel]]$support>supportLevel) 
          + (OmegaCor[[rlevel]]$support<(1-supportLevel))>0)*OmegaCor[[rlevel]]$mean
# reorder species matrix
plotorder <- order( postBeta$mean[5,], decreasing = TRUE )
toPlot <- toPlot[ plotorder, plotorder]
# reorder automatically
library(lessR)
mynewcor <- corReorder( toPlot, order="hclust", nclusters=4 )
# windows(12,12)
corrplot( mynewcor, method = "color", 
          col = colorRampPalette(c("blue","white","red"))(200),
          title = paste("random effect level:", m$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )



######################################################
######################## 2017 ########################
######################################################

## pick a year
microbes_2017 <- microbes %>% 
  filter( year==2017 )

microbes_2017 <- microbes_2017 %>%  
  mutate(quadrat_id = as.numeric(quadrat_id))

# explanatory data
environmental_coord <- read.csv("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R Code and Analysis/output data/merged_explanatory.csv")

## pick a year
keep_2017 <- c("2017")
environmental_2017 <- environmental_coord %>% 
  dplyr::filter(year %in% keep_2017)

### microbial sites for 2017 are: choked_inner, choked_sandspit, pruth_bay, pruth_pocket, triquet_north, triquet_south
keep_sites <- c("choked_inner", "choked_sandspit", "pruth_pocket", "pruth_bay", "triquet_north", "triquet_south")
environmental_2017_sites <- environmental_2017 %>% 
  dplyr::filter(site %in% keep_sites)

environmental_2017_sites <- environmental_2017_sites %>%  
  mutate(quadrat_id = as.numeric(quadrat_id))

### merge microbes and env_coord
microbes_metadata <- inner_join(environmental_2017_sites, microbes_2017, by = c("region","site", "quadrat_id"))
View(microbes_metadata)

# have to exclude rows with NA values
microbes_metadata <-  na.omit(microbes_metadata)

# response data
names(microbes_metadata)
Ymicrobes <- microbes_metadata %>% 
  dplyr::select(27:46)
names(Ymicrobes)

## select columns for fixed and random effects
# explanatory data = abiotic + seagrass metrics
Xenvironmental <- microbes_metadata %>% 
  dplyr::select(temperature,salinity,depth,bed_area_m2, quadrat_biomass_g,quadrat_lai,quadrat_microepiphyte_mg,quadrat_macroalgae_g) #do,

# random effects (quadrat, site, region level)
studyDesign <- microbes_metadata %>% 
  dplyr::select(region, site, quadrat)

studyDesign <-  as.data.frame(studyDesign)
studyDesign$region <- factor(studyDesign$region)
studyDesign$site <- factor(studyDesign$site)
studyDesign$quadrat <- factor(studyDesign$quadrat)

# spatial data from lats and longs
# make longitude positive
spatial <- microbes_metadata %>% 
  dplyr::select(lat, long)
spatial$long <- spatial$long-min(spatial$long)
spatial <- as.data.frame(spatial)
# formats
Ymicrobes <- as.matrix(Ymicrobes)
Xenvironmental <- as.data.frame(Xenvironmental)
#View(Xenvironmental)

# random effects structure
rL1 = HmscRandomLevel( units = unique(studyDesign$quadrat) ) #quadrat level
rL2 = HmscRandomLevel( units = unique(studyDesign$site) )
rL3 = HmscRandomLevel( units = unique(studyDesign$region) )
rL4 = HmscRandomLevel( sData = spatial )


## formula for fixed effects
XFormula = ~( temperature+salinity+depth+bed_area_m2+ quadrat_biomass_g+quadrat_lai+quadrat_microepiphyte_mg+quadrat_macroalgae_g)

# the model
mmicrobe_2015 <- Hmsc( Y = Ymicrobes, 
                       XData = Xenvironmental, XFormula = XFormula, 
                       distr = "poisson", studyDesign = studyDesign,
                       ranLevels = list(quadrat=rL1, 
                                        #site=rL2,
                                        region=rL3,
                                        site=rL4) )

### running the model
thin = 100
samples = 100
nChains = 2 # have to compute multiple chains to check if they all yield similar results
set.seed(1)

# Run MCMC chains. took at least 12 hours on KS laptop

mod <- sampleMcmc(mmicrobe_2015, samples = samples , transient = 50, 
                  thin = thin, verbose = 10, nChains = nChains) # nParallel = 3 when running more than 1 Chain for the real model = running more than one chain allows us to check if all arer doing the same thing
#(hM, samples = 2000, transient = 100, thin = 100, verbose = 10)

# saving
model <- "grazer2017"
ModelDir <- paste0( here::here(), "/R Code and Analysis/output data")
filename = file.path(ModelDir, paste("model_",as.character(model),
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples),".Rdata",sep = ""))
save(mod,file=filename)


load(filename)


## model results
m <- mod

# We evaluate MCMC convergence in terms of four kinds of parameters that we are especially interested in: the species niches `Beta`, the influence of traits on species niches `Gamma`, residual species associations `Omega`, and the strength of phylogenetic signal `rho`. As there are `r ns` species, the matrix `Omega` is a `r ns` x `r ns` matrix with `r ns*ns` elements. To avoid excessive computational times, we evaluate MCMC convergence for `Omega` only for a subsample of 100 randomly selected species pairs.
mpost = convertToCodaObject(m)



## Assess model fit
preds = computePredictedValues(m)
MF <- evaluateModelFit(hM = m, predY = preds)


## parameter estimates
postBeta = getPostEstimate(m, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature", "do", "salinity", "quadrat_biomass_g", 
                              "quadrat_lai", "quadrat_microepiphyte_mg" ),
                            levels = c("intercept", "temperature", "do", "salinity", "quadrat_biomass_g", 
                                       "quadrat_lai", "quadrat_microepiphyte_mg" ), 
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 7), 
                          levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)],
                          ordered = TRUE)

windows(12,4)
ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")

# pull taxa with positive estimated response to temperature anomaly
taxa.pos.anom <- as.character(pos.neg[ pos.neg$parameter=='quadrat_lai' & pos.neg$value>0, 'species'])

# calculate mean and variance of parameter esimates
pbdf <- data.frame( t(postBeta$mean), taxon=colnames(postBeta$mean) )
names(pbdf) <- c("intercept","temp","do","sal",
                 "biomass","LAI","microepi","taxon")
# ## Add some basic trait information
# pbdf$alga <- "alga"
# pbdf$alga[c(3,15,20,56,57,62,68,82)] <- "invert"
# # More specific groups


coef_plot <- pbdf %>% 
  select( -taxon ) %>%
  gather( coefficient, value )

windows(6,4)
ggplot( coef_plot, aes(y=value, x=factor(1))) +
  facet_wrap(~coefficient, scales="free_y") +
  geom_hline(yintercept=0) +
  stat_summary( fun.data=mean_cl_boot, geom='errorbar',
                size=1, width=0.1, position=position_dodge(width=0.9) )  +
  stat_summary( fun.y=mean, geom='point',
                size=3, position=position_dodge(width=0.9) )  +
  ylab('coefficient') + xlab('') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




## variance partitioning
VP = computeVariancePartitioning(m) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP)

VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c("temp","do","sal",
                           "biomass","LAI","microepi",
                           "quadrat","region","site"), 
                         levels = rev(c("temp","do","sal",
                                        "biomass","LAI","microepi",
                                        "quadrat","region","site")), 
                         ordered = TRUE)) %>% 
  gather(key = species, value = variance, -effect) %>% 
  group_by(species) %>% 
  mutate(tempR2 = variance[effect == "temp"])


hold <- VP.df %>% filter(effect == "temp") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(m$Y)[order(colSums(m$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(m$Y))

windows(8,5)
ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkgreen","forestgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")






## associations
OmegaCor = computeAssociations(m)
supportLevel = 0.95
# choose the random variable to plot
rlevel = 3
toPlot = ((OmegaCor[[rlevel]]$support>supportLevel) 
          + (OmegaCor[[rlevel]]$support<(1-supportLevel))>0)*OmegaCor[[rlevel]]$mean
# reorder species matrix
plotorder <- order( postBeta$mean[5,], decreasing = TRUE )
toPlot <- toPlot[ plotorder, plotorder]
# reorder automatically
library(lessR)
mynewcor <- corReorder( toPlot, order="hclust", nclusters=4 )
# windows(12,12)
corrplot( mynewcor, method = "color", 
          col = colorRampPalette(c("blue","white","red"))(200),
          title = paste("random effect level:", m$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )


