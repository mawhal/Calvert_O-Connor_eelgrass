###################################################################
###  Hakai + O'Connor Seagrass Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code produces summaries of eelgrass variables
###    inlcuding those measured at shoot and quadrat levels
###
###  Question: do we have information about the shoots that were
###            swabbed for microbes?
###
### 
###  started by Matt Whalen  on 19 April 2020
###  
###################################################################


## libraries
library(tidyverse)


# what to do with 2015? Different number of shoots, not included in quadrat level dataset
#   can maybe use single shoot data here??

## read data
# quadrats
qall <- read_csv( "output_data/emily_hakai_quad_combined_20200205_bia_reviewed.csv" )
qsel <- qall %>% select( year, site, sample=quadrat_id, shoot.density=quadrat_shoot_density, lai=quadrat_lai, biomass=quadrat_biomass_g )
# shoots
sall <- read_csv( "../Data/R Code for Data Prep/master data/O'Connor_hakai_seagrass_MASTER_shoots.csv" )
ssel <- sall %>% select( year, site, sample=id, shoot.length, shoot.width )
# bed area and site depth
abiotic <- read_csv( "../metadata/00_Hakai_UBC_metadata_MASTER - geolocation_site.csv" )
abiotic$site_name <- tolower(gsub( " ","_", abiotic$site_name ))
abiotic$site_name[ abiotic$site_name =="inner_choked"] <- "choked_inner"
abiotic$site_name[ abiotic$site_name =="sandspit"] <- "choked_sandspit"
abiotic <- abiotic %>% 
  filter( project=="O'Connor" ) %>% 
  select( site=site_name, lat, long, depth_m = `depth (m, chart datum)`, area_m2=area_m2_2016 )

# merge
qs <- full_join( qsel, ssel )
qs[ apply( qs,1, function(z) any(is.na(z)) ), ] 
 # complete cases
qs <- qs[ complete.cases( qs ), ]
qsa <- left_join( qs, abiotic )


# 
psych::pairs.panels( select(qsa,shoot.density,biomass,shoot.length,shoot.width,lai), scale = F )
psych::pairs.panels( qsa %>% select(shoot.density,biomass,shoot.length,shoot.width,lai) %>% 
                       mutate(shoot.density=log(shoot.density), biomass=log(biomass), lai=log(lai)), scale = F )
pp <- select(qsa, depth_m, area_m2, shoot.density, biomass, shoot.length, shoot.width)
ppl <- pp %>% mutate( area_m2=log(area_m2), shoot.density=log(shoot.density), biomass=log(biomass) )
pp <- select(qsa, shoot.density, biomass, shoot.length, shoot.width)
ppl <- pp %>% mutate(  shoot.density=log(shoot.density), biomass=log(biomass) )
psych::pairs.panels( ppl, scale = T )


# PCA
pca <- prcomp( ppl )
summary(pca)
biplot(pca,choices = c(1,2), pc.biplot = F, scale=0.5 )
biplot(pca,choices = c(2,3))
ppp <- data.frame( qsa, pca$x )
ggplot( ppp, aes(x=PC1,y=PC2,col=site)) + geom_point()
library( ggfortify )
autoplot( prcomp(ppl), data=ppp, col='site', 
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )
ggsave( "figs/eelgrass_PCA_biotic12.png", width=5, height=3 )
autoplot( prcomp(ppl), data=ppp, col='site', x=2, y=3,
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )
ggsave( "figs/eelgrass_PCA_biotic23.png", width=5, height=3 )
autoplot( prcomp(ppl), data=ppp, col='site', x=1, y=4,
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )
ggsave( "figs/eelgrass_PCA_14.png", width=5, height=3 )
autoplot( prcomp(ppl), data=ppp, col='site', x=4, y=5,
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )

# repeat but average across samples within a site
ppm <- qsa %>% 
  group_by(year,site) %>% 
  summarize_all( mean ) %>% 
  ungroup()
pp <- select(ppm, depth_m, area_m2, shoot.density, biomass, shoot.length, shoot.width)
ppl <- pp %>% mutate( area_m2=log(area_m2), shoot.density=log(shoot.density), biomass=log(biomass) )
pp <- select(ppm, shoot.density, biomass, shoot.length, shoot.width)
ppl <- pp %>% mutate( shoot.density=log(shoot.density), biomass=log(biomass) )

psych::pairs.panels( ppl, scale = T )


# PCA
pca <- prcomp( ppl )
summary(pca)
autoplot( prcomp(ppl), data=ppm, col='site', 
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )
ggsave( "figs/eelgrass_PCA_bioticmean12.png", width=5, height=3 )
autoplot( prcomp(ppl), data=ppm, x=2, y=3,
          col='site', 
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )
ggsave( "figs/eelgrass_PCA_bioticmean23.png", width=5, height=3 )

qsa$site <- fct_reorder(qsa$site,qsa$shoot.length)
ggplot( qsa, aes(x=site,y=shoot.length)) + geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ppm$site <- fct_reorder(ppm$site,ppm$shoot.length)
ggplot( ppm, aes(x=site,y=shoot.length)) + geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(axis.text=element_text(size=6))
vegan::make.cepnames(unique(ppm$site))
ggsave( "figs/eelgrass_shoot_length.png", width=4,height=3 )



# replace LAI with other variables
pp <- select(qsa, depth_m,  biomass, lai ) #area_m2,
ppl <- pp %>% mutate( biomass=log(biomass), lai=log(lai) ) #area_m2=log(area_m2), 
# PCA
pca <- prcomp( ppl )
summary(pca)
autoplot( prcomp(ppl), data=qsa, col='site', 
          loadings = T, loadings.colour="blue",
          loadings.label = T, loadings.label.colour="black", 
          loadings.label.size = 2,  
          scale = 1  )
qsa$site <- fct_reorder(qsa$site,qsa$biomass)
ggplot( qsa, aes(x=site,y=biomass)) + geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_y_continuous(trans="log10") +
  theme(axis.text=element_text(size=6))
ggsave( "figs/eelgrass_biomass.png", width=4,height=3 )
