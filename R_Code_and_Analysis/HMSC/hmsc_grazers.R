####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will run HMSC on the input grazer dataset
###  code by Matt Whalen
###  started on   04 February 2020
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
m <- read_csv( "../Data/R Code for Data Prep/master data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan

### select sites
# filter taxa and sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon4))
# # Four core sites
# mfilt <- m %>%
#   filter( is.na(remove), !is.na(taxon4),
#           site %in% c("inner choked","sandspit",
#                       "triquet north","triquet south") )


## pick a year
muse <- mfilt %>% 
  filter( year==2017 )

# which data to use
muse <- muse


# summarize taxon counts per sample
m.sum <- muse %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon4 ) %>% 
  summarize( abundance=length(size) )



# make a community dataset
m.meta <- m.sum %>% 
  spread( taxon4, abundance, fill=0 )

meta <- data.frame(m.meta[,c(1:3)])
# meta$year <- factor( meta$year, ordered=T )
comm <- data.frame(m.meta[,-c(1:3)])
names(comm) <- make.cepnames( names(comm) )

# rename quadrats
meta <- meta %>% 
  mutate( quadrat_id = as.numeric(meta$sample), quadrat = paste( site,quadrat_id, sep="_" ) )


# get rid of some taxa
# rename to Ygrazer
Ygrazer <- comm
Ygrazer <- Ygrazer[ , rev(order(colSums(Ygrazer))) ]

# abundance cutoff?
Ygrazer <- Ygrazer[, colSums(Ygrazer) > 3]

###
# explanatory data
explanatory <- read_csv("output data/merged_explanatory_20200214.csv")

#  
ab <- explanatory
#filter out missing quads
ab <- ab %>% filter( !is.na(quadrat) )
# merge with metadata
X <- left_join( meta, ab, by=c("year","site","quadrat_id","quadrat") )
# merge lat long back in
ab %>% mutate(long=-long) %>% select(site,lat,long) %>% distinct()

# filter out NA values from grazer dataset
nas <- !is.na(X$quadrat_shoot_density) & !is.na(X$quadrat_biomass_g) & !is.na(X$quadrat_lai)
Ygrazer <- Ygrazer[nas,]
# filter out missing shoot densities
X <- X %>%  filter( !is.na(quadrat_shoot_density),!is.na(quadrat_biomass_g),!is.na(quadrat_lai)  )

## select columns for fixed and random effects
# explanatory data = abiotic + seagrass metrics
Xenvironmental <- X %>% 
  select( temperature, salinity, 
          depth=depth..m..chart.datum., area=bed_area_m2,
          biomass=quadrat_biomass_g, lai=quadrat_lai, 
          microepi=quadrat_microepiphyte_mg, 
          macro=quadrat_macroalgae_g, 
          ) 

# random effects (quadrat, site, region level)
studyDesign <- X %>% 
  mutate( region=unlist(lapply(strsplit(site,"_"),function(z)z[1])),
          long=-long ) %>% 
  select( region, lat, long, site, quadrat )
studyDesign <-  as.data.frame(studyDesign)
studyDesign$region <- factor(studyDesign$region)
studyDesign$site <- factor(studyDesign$site)
studyDesign$quadrat <- factor(studyDesign$quadrat)

# spatial data from lats and longs
# make longitude positive
spatial <- studyDesign  %>% select(site,lat,long) %>% distinct()
spatial <- as.data.frame(spatial)
rownames(spatial) <- spatial$site
spatial <- spatial[,2:3]


# formats
Ygrazer <- as.matrix(Ygrazer)
Xenvironmental <- as.data.frame(Xenvironmental)



# random effects structure
rL1 = HmscRandomLevel( units = unique(studyDesign$site_quad) ) #quadrat level
rL2 = HmscRandomLevel( units = unique(studyDesign$site) )
rL3 = HmscRandomLevel( units = unique(studyDesign$region) )
rL4 = HmscRandomLevel( sData = spatial )




## formula for fixed effects
XFormula = ~( temperature + salinity + depth + area + biomass +  
                lai + microepi + macro  ) # add bed area when available


# the model
mgrazer <- Hmsc( Y = Ygrazer, 
                   XData = Xenvironmental, XFormula = XFormula, 
                   distr = "poisson", studyDesign = studyDesign,
                   ranLevels = list(quadrat=rL1, 
                                    #Site=rL2,
                                    region=rL3,
                                    site=rL4) )



### running the model
thin = 100
samples = 1000
transient = 100
nChains = 4 
set.seed(1)

# Run MCMC chains. took at least 12 hours on KS laptop

mod <- sampleMcmc(mgrazer, samples = samples , transient = transient, 
                  thin = thin, verbose = 10, nChains = nChains) # nParallel = 3 when running more than 1 Chain for the real model = running more than one chain allows us to check if all arer doing the same thing


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

#
plot( mpost$Beta )


## Assess model fit
preds = computePredictedValues(m)
MF <- evaluateModelFit(hM = m, predY = preds)


## parameter estimates
postBeta = getPostEstimate(m, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth", "area", "biomass",
                                "lai", "microepi", "macro" ),
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
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
names(pbdf) <- c("intercept","temp","sal","depth","area",
                 "biomass","LAI","microepi","macro", "taxon")
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
  mutate(effect = factor(c("temp","sal","depth","area",
                            "biomass","LAI","microepi","macro",
                           "quadrat","region","site"),
                         levels = rev(c("temp","sal","depth","area",
                                        "biomass","LAI","microepi","macro",
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
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "")






## associations
OmegaCor = computeAssociations(m)
supportLevel = 0.95
# choose the random variable to plot
rlevel = 1
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


