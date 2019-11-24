## Test HMSC run at Seagrass Retreat
# started 23 November 2019

# load libraries
library(Hmsc)
library(tidyverse)
library(viridis)
library(corrplot)

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
# get rid of some taxa
comm <- comm[, -c(1,12,17,19,25)] # get rid of Acidiacea, Gastropoda, Lottoidea, Neogastropoda, Rissoidae
# rename to Ygrazer
Ygrazer <- comm
# order columns based on total abundance
Ygrazer <- Ygrazer[ , rev(order(colSums(Ygrazer))) ]

Ymicrobe <- read.csv("microbes_2017.csv")

Yfish

# explanatory data
X <- read.csv("Hakai_2017_Eelgrass_Biometrics_final.csv")
names(X)[3] <- "sample"

X1 <- X[,-c(4:50)]

X2 <- X1 %>%
  select(-detritus_wet_mass, -epi_algae_wet_mass, -seagrass_wet_mass, -avg_shoot_surface_area_m2, -shoots_m.2, -notes, -date)

# Bio-ORACLE data
Xab <- read.csv("biooracle_2017_small.csv")
# run PCA on Bio_ORACLE data
pca <- prcomp( Xab[,-1] )

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
fixnames <- c("detritus","macroepiph","eel_biomass",
              "LAI","microepiph","nitrate_range",
              "tempmean","temprange","curveloc","nitratemean",
              "salinityrange","salinitymean","dissoxmean",
              "dissoxrange")
XFormula = ~(detritus_dry_mass + epi_algae_dry_mass + seagrass_dry_mass + 
               lai + microepiphyte_wt_mg)
fixnames <- c("detritus","macroepiph","eel_biomass",
              "LAI","microepiph")

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

## Assess model fit
preds = computePredictedValues(mod)
MF <- evaluateModelFit(hM = mod, predY = preds)




## parameter estimates
postBeta = getPostEstimate(mod, parName = "Beta")
# windows(5,8)
# plotBeta(m, post = postBeta, param = "Sign", supportLevel = 0.95, mar=c(7,11,0,6))
postBeta$mean[, c("Gammaridean","Isopoda")]

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor( c("intercept", fixnames ),
                            levels = c("intercept", fixnames ), 
                            ordered = TRUE)
pos.neg$taxon <- factor(rep(colnames(postBeta$mean), each = length(fixnames)+1), 
                          levels = colnames(mod$Y)[order(colSums(mod$Y),decreasing = TRUE)],
                          ordered = TRUE)

windows(12,4)
ggplot(pos.neg, aes(y = parameter, x = taxon, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 45,size=9, hjust = 1))+
  xlab(label = "")+
  ylab(label = "")
ggsave( "hmsc_grazer_betas.jpg" )


## variance partitioning
VP = computeVariancePartitioning(mod) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(mod, VP = VP)

VP.df <- as.data.frame(VP$vals) %>% 
  mutate(effect = factor(c(fixnames, "quadrat", "site", "region" ), 
                         levels = rev(c(fixnames, "quadrat", "site", "region" )), 
                         ordered = TRUE)) %>% 
  gather(key = taxon, value = variance, -effect) %>% 
  group_by(taxon) %>% 
  mutate(massR2 = variance[effect == "eel_biomass"])

# could use something life this to add variance explained for each fixed and random effect
hold <- VP.df %>% filter(effect == "eel_biomass") %>% arrange(desc(massR2))

VP.df$taxon <- factor(VP.df$taxon, 
                        levels = colnames(mod$Y)[order(colSums(mod$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), taxon = colnames(mod$Y))

windows(8,5)
ggplot(VP.df,aes(y = variance, x = taxon, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","magenta", viridis(length(fixnames))), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")
ggsave( "hmsc_grazer_varpart.jpg", width=250, height=175, units="mm",dpi=300 )



## Co-occurence

## associations
OmegaCor = computeAssociations(mod)
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
          title = paste("random effect level:", mod$rLNames[rlevel]), mar=c(0,0,0.5,0), tl.cex=0.6 )
