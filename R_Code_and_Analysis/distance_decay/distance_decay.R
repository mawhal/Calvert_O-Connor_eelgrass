################Preliminary Analyses to explore the impact of distance between sites on grazer dissimilarity
##Started by Coreen April 2020
##This script  makes distance matrices between sites and plots Bray-Curtis grazer dissimilarity between each pair of sites against dstance between each pair of sites for 2014-207
##For some reason i did each year totally separately. No idea why
## updated by Whalen on 15 May 2020. Keeping analysis separated by date for now
## updated by Bia on 07 July 2020 to add ASV level microbes, update for corrected microbial tables and add title

library(vegan)
library(tidyverse)
# library(distances)
library(fields)
library(cowplot)


# Geographic distance between site pairs ---------------------------------------
Hakaispatial <- read.csv("metadata/00_Hakai_UBC_metadata_MASTER - geolocation_site.csv")
Hakaispatial1 <- Hakaispatial %>%
  select(site_name, lat, long)

##Change site names that don't match master grazer data
Hakaispatial1$site_name <- recode(Hakaispatial1$site_name, 
                                  "inner_choked" = "choked_inner", "sandspit" = "choked_sandspit")
coords <- Hakaispatial %>%  select( long, lat )
spdf <- SpatialPointsDataFrame( coords, Hakaispatial1 )

## Make distance matrix
Hakai.distance <- rdist.earth( coords[,c('long','lat')] )
Hakai.geog <- as.data.frame(Hakai.distance)
# Hakai.distance<- distances(Hakaispatial1, id_variable = "site_name", dist_variables = NULL)
# Hakai.geog <-as.data.frame(as.matrix(Hakai.distance))

###Now make the distance matricies long
# Hakai.geog$Sites1 <- rownames(Hakai.geog)
colnames(Hakai.geog) <- Hakaispatial1$site
Hakai.geog$Sites1 <- Hakaispatial1$site


##Data frame of distances between all site pairs

Hakai.geographic.distance <- Hakai.geog %>% 
  gather(Sites2, Geog_Distance, - Sites1)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = T)




#### MACROEUKARYOTES
# pick a taxonomic level
level <- "finest"
# folder location
path <- "R_Code_and_Analysis/betadiversity/Bray-Curtis/"

##### 2016 grazers ------------------------------------------------------------
dist16 <- read_csv( paste0(path,"2016_macroeuk_braycurtis_",level,".csv") )
# just take the upper portion of the distance matrix so we don't repeat the numbers
dist16 <- as.data.frame(as.matrix(dist16))
dist16[lower.tri(dist16,diag = T)] <- NA
meta16 <- read_csv( paste0(path,"2016_macroeuk_metadata.csv") )
meta16$site <- unlist( lapply( strsplit( meta16$sample, split = "_"), function(z) paste(z[1:(length(z)-1)],collapse = "_") ) )
# script used to calculate Bray-Curtis have shorten column names to make the distance matrix more compact.
meta16$samp.short <- vegan::make.cepnames(meta16$sample)

###Now make the distance matrix long
dist16$Sites1 <- colnames(dist16)
dist16.collapse <- dist16 %>% 
  gather(Sites2, Community_Distance, -Sites1)

# add sites from metadata
dist16.sites <- left_join( dist16.collapse, select(meta16, site, sample, samp.short), by=c("Sites1" = "samp.short") )
dist16.sites <- left_join( dist16.sites, select(meta16, site, sample, samp.short), by=c("Sites2" = "samp.short") )
dist16.sites <- dist16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16.distance <- dist16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance <- left_join(dist16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

### plots 
# Graph1 <- Hakai.2016.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 Macroeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 Macroeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)

BC16 <- Hakai.2016.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 macroeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)



##### 2017 grazers ------------------------------------------------------------
dist17 <- read_csv( paste0(path,"2017_macroeuk_braycurtis_",level,".csv") )
# just take the upper portion of the distance matrix so we don't repeat the numbers
dist17 <- as.data.frame(as.matrix(dist17))
dist17[lower.tri(dist17,diag = T)] <- NA
meta17 <- read_csv( paste0(path,"2017_macroeuk_metadata.csv") )
meta17$site <- unlist( lapply( strsplit( meta17$sample, split = "_"), function(z) paste(z[1:(length(z)-1)],collapse = "_") ) )
# script used to calculate Bray-Curtis have shorten column names to make the distance matrix more compact.
meta17$samp.short <- vegan::make.cepnames(meta17$sample)


###Now make the distance matricies long

dist17$Sites1 <- colnames(dist17)
dist17.collapse <- dist17 %>% 
  gather(Sites2, Community_Distance, - Sites1)

# add sites from metadata
dist17.sites <- left_join( dist17.collapse, select(meta17, site, sample, samp.short), by=c("Sites1" = "samp.short") )
dist17.sites <- left_join( dist17.sites, select(meta17, site, sample, samp.short), by=c("Sites2" = "samp.short") )
dist17.sites <- dist17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist17.distance <- dist17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance <- left_join(dist17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

### plots
# Graph1 <- Hakai.2017.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 Macroeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 Macroeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 macroeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


##### 2015 grazers ------------------------------------------------------------
dist15 <- read_csv( paste0(path,"2015_macroeuk_braycurtis_",level,".csv") )
# just take the upper portion of the distance matrix so we don't repeat the numbers
dist15 <- as.data.frame(as.matrix(dist15))
dist15[lower.tri(dist15,diag = T)] <- NA
meta15 <- read_csv( paste0(path,"2015_macroeuk_metadata.csv") )
meta15$site <- unlist( lapply( strsplit( meta15$sample, split = "_"), function(z) paste(z[1:(length(z)-1)],collapse = "_") ) )
# script used to calculate Bray-Curtis have shorten column names to make the distance matrix more compact.
meta15$samp.short <- vegan::make.cepnames(meta15$sample)

###Now make the distance matrix long
dist15$Sites1 <- colnames(dist15)
dist15.collapse <- dist15 %>% 
  gather(Sites2, Community_Distance, - Sites1)

# add sites from metadata
dist15.sites <- left_join( dist15.collapse, select(meta15, site, sample, samp.short), by=c("Sites1" = "samp.short") )
dist15.sites <- left_join( dist15.sites, select(meta15, site, sample, samp.short), by=c("Sites2" = "samp.short") )
dist15.sites <- dist15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist15.distance <- dist15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance <- left_join(dist15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

### plots
# Graph1 <- Hakai.2015.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 Macroeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 Macroeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 macroeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)



##### 2014 grazers ------------------------------------------------------------
dist14 <- read_csv( paste0(path,"2014_macroeuk_braycurtis_",level,".csv") )
# just take the upper portion of the distance matrix so we don't repeat the numbers
dist14 <- as.data.frame(as.matrix(dist14))
dist14[lower.tri(dist14,diag = T)] <- NA
meta14 <- read_csv( paste0(path,"2014_macroeuk_metadata.csv") )
meta14$site <- unlist( lapply( strsplit( meta14$sample, split = "_"), function(z) paste(z[1:(length(z)-1)],collapse = "_") ) )
# script used to calculate Bray-Curtis have shorten column names to make the distance matrix more compact.
meta14$samp.short <- vegan::make.cepnames(meta14$sample)

###Now make the distance matrix long
dist14$Sites1 <- colnames(dist14)
dist14.collapse <- dist14 %>% 
  gather(Sites2, Community_Distance, - Sites1)

# add sites from metadata
dist14.sites <- left_join( dist14.collapse, select(meta14, site, sample, samp.short), by=c("Sites1" = "samp.short") )
dist14.sites <- left_join( dist14.sites, select(meta14, site, sample, samp.short), by=c("Sites2" = "samp.short") )
dist14.sites <- dist14.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist14.distance <- dist14.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2014.distance <- left_join(dist14.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

### plots
# Graph1 <- Hakai.2014.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2014 Macroeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2014.distance %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2014 Macroeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC14 <- Hakai.2014.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2014 macroeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)

windows(12,3)
macro <- cowplot::plot_grid( BC14, BC15, BC16, BC17, ncol=4)
ggsave( paste0("R_Code_and_Analysis/distance_decay/BCdecay_macroeuk_",level,".png"), width = 12, height = 3  )





#### MICROBES

#### Prokaryotes - 16S

## ASV LEVEL
# pick a year
year <- 2015

#load 16S microbial distance matrix ASV
bc_16S15_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_braycurtis.csv") )
# bc_16S15_ASV <- bc_16S15_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_16S15_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_16S15_ASV$Sites1 <- colnames(bc_16S15_ASV)
dist16S15.collapse <- bc_16S15_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S15.sites <- left_join( dist16S15.collapse, select(bc_16S15_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S15.sites <- left_join( dist16S15.sites, select(bc_16S15_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S15.sites <- dist16S15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S15.distance <- dist16S15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance.16S <- left_join(dist16S15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2015.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2016

#load 16S microbial distance matrix ASV
bc_16S16_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_braycurtis.csv") )
# bc_16S16_ASV <- bc_16S16_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_16S16_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_16S16_ASV$Sites1 <- colnames(bc_16S16_ASV)
dist16S16.collapse <- bc_16S16_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S16.sites <- left_join( dist16S16.collapse, select(bc_16S16_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S16.sites <- left_join( dist16S16.sites, select(bc_16S16_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S16.sites <- dist16S16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S16.distance <- dist16S16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance.16S <- left_join(dist16S16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2016.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC16 <- Hakai.2016.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2017

#load 16S microbial distance matrix ASV
bc_16S17_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_braycurtis.csv") )
# bc_16S17_ASV <- bc_16S17_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_16S17_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_16S17_ASV$Sites1 <- colnames(bc_16S17_ASV)
dist16S17.collapse <- bc_16S17_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S17.sites <- left_join( dist16S17.collapse, select(bc_16S17_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S17.sites <- left_join( dist16S17.sites, select(bc_16S17_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S17.sites <- dist16S17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S17.distance <- dist16S17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance.16S <- left_join(dist16S17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2017.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2018

#load 16S microbial distance matrix ASV
bc_16S18_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_braycurtis.csv") )
# bc_16S18_ASV <- bc_16S18_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_16S18_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_16S_",year,"_metadata.csv") )




###Now make the distance matrices long
bc_16S18_ASV$Sites1 <- colnames(bc_16S18_ASV)
dist16S18.collapse <- bc_16S18_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S18.sites <- left_join( dist16S18.collapse, select(bc_16S18_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S18.sites <- left_join( dist16S18.sites, select(bc_16S18_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S18.sites <- dist16S18.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S18.distance <- dist16S18.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2018.distance.16S <- left_join(dist16S18.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2018.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2018.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC18 <- Hakai.2018.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2018 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


windows(12,3)
title <-  ggdraw() + draw_label("ASV level",fontface = 'bold', size = 14, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
ASV_16S_raw <- cowplot::plot_grid( BC15, BC16, BC17, BC18, ncol=4)
ASV_16S <- plot_grid(title, plots,ncol = 1,rel_heights = c(0.05, 1)) # rel_heights values control vertical title margins
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_prokaryote_ASV.png"), width = 12, height = 5  )


## GENUS LEVEL
# pick a year
year <- 2015

#load 16S microbial distance matrix GENUS
bc_16S15_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_braycurtis.csv") )
# bc_16S15_genus <- bc_16S15_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_16S15_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_16S15_genus$Sites1 <- colnames(bc_16S15_genus)
dist16S15.collapse <- bc_16S15_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S15.sites <- left_join( dist16S15.collapse, select(bc_16S15_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S15.sites <- left_join( dist16S15.sites, select(bc_16S15_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S15.sites <- dist16S15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S15.distance <- dist16S15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance.16S <- left_join(dist16S15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2015.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2016

#load 16S microbial distance matrix GENUS
bc_16S16_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_braycurtis.csv") )
# bc_16S16_genus <- bc_16S16_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_16S16_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_16S16_genus$Sites1 <- colnames(bc_16S16_genus)
dist16S16.collapse <- bc_16S16_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S16.sites <- left_join( dist16S16.collapse, select(bc_16S16_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S16.sites <- left_join( dist16S16.sites, select(bc_16S16_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S16.sites <- dist16S16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S16.distance <- dist16S16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance.16S <- left_join(dist16S16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2016.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC16 <- Hakai.2016.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2017

#load 16S microbial distance matrix GENUS
bc_16S17_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_braycurtis.csv") )
# bc_16S17_genus <- bc_16S17_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_16S17_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_16S17_genus$Sites1 <- colnames(bc_16S17_genus)
dist16S17.collapse <- bc_16S17_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S17.sites <- left_join( dist16S17.collapse, select(bc_16S17_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S17.sites <- left_join( dist16S17.sites, select(bc_16S17_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S17.sites <- dist16S17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S17.distance <- dist16S17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance.16S <- left_join(dist16S17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2017.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2018

#load 16S microbial distance matrix GENUS
bc_16S18_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_braycurtis.csv") )
# bc_16S18_genus <- bc_16S18_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_16S18_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_metadata.csv") )




###Now make the distance matrices long
bc_16S18_genus$Sites1 <- colnames(bc_16S18_genus)
dist16S18.collapse <- bc_16S18_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S18.sites <- left_join( dist16S18.collapse, select(bc_16S18_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S18.sites <- left_join( dist16S18.sites, select(bc_16S18_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S18.sites <- dist16S18.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S18.distance <- dist16S18.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2018.distance.16S <- left_join(dist16S18.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2018.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2018.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC18 <- Hakai.2018.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2018 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


title <-  ggdraw() + draw_label("genus level",fontface = 'bold', size = 14, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
genus_16S_raw  <- cowplot::plot_grid( BC15, BC16, BC17, BC18, ncol=4)
genus_16S <- plot_grid(title, genus_16S_raw,ncol = 1,rel_heights = c(0.05, 1)) # rel_heights values control vertical title margins
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_prokaryote_genus.png"), width = 12, height = 5  )



### FAMILY LEVEL
# pick a year
year <- 2015

#load 16S microbial distance matrix family
bc_16S15_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_braycurtis.csv") )
# bc_16S15_family <- bc_16S15_family %>% 
#   dplyr::rename("sample" = "X1")
bc_16S15_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_metadata.csv") )




###Now make the distance matrices long
bc_16S15_family$Sites1 <- colnames(bc_16S15_family)
dist16S15.collapse <- bc_16S15_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S15.sites <- left_join( dist16S15.collapse, select(bc_16S15_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S15.sites <- left_join( dist16S15.sites, select(bc_16S15_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S15.sites <- dist16S15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S15.distance <- dist16S15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance.16S <- left_join(dist16S15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2015.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2016

#load 16S microbial distance matrix family
bc_16S16_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_braycurtis.csv") )
# bc_16S16_family <- bc_16S16_family %>% 
#   dplyr::rename("sample" = "X1")
bc_16S16_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_metadata.csv") )




###Now make the distance matrices long
bc_16S16_family$Sites1 <- colnames(bc_16S16_family)
dist16S16.collapse <- bc_16S16_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S16.sites <- left_join( dist16S16.collapse, select(bc_16S16_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S16.sites <- left_join( dist16S16.sites, select(bc_16S16_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S16.sites <- dist16S16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S16.distance <- dist16S16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance.16S <- left_join(dist16S16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2016.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC16 <- Hakai.2016.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2017

#load 16S microbial distance matrix family
bc_16S17_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_braycurtis.csv") )
# bc_16S17_family <- bc_16S17_family %>% 
#   dplyr::rename("sample" = "X1")
bc_16S17_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_16S17_family$Sites1 <- colnames(bc_16S17_family)
dist16S17.collapse <- bc_16S17_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S17.sites <- left_join( dist16S17.collapse, select(bc_16S17_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S17.sites <- left_join( dist16S17.sites, select(bc_16S17_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S17.sites <- dist16S17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S17.distance <- dist16S17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance.16S <- left_join(dist16S17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2017.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2018

#load 16S microbial distance matrix family
bc_16S18_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_braycurtis.csv") )
# bc_16S18_family <- bc_16S18_family %>% 
#   dplyr::rename("sample" = "X1")
bc_16S18_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_16S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_16S18_family$Sites1 <- colnames(bc_16S18_family)
dist16S18.collapse <- bc_16S18_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist16S18.sites <- left_join( dist16S18.collapse, select(bc_16S18_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist16S18.sites <- left_join( dist16S18.sites, select(bc_16S18_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist16S18.sites <- dist16S18.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist16S18.distance <- dist16S18.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2018.distance.16S <- left_join(dist16S18.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2018.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 prokaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2018.distance.16S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 prokaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC18 <- Hakai.2018.distance.16S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2018 prokaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


title <-  ggdraw() + draw_label("family level",fontface = 'bold', size = 14, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
family_16S_raw <- cowplot::plot_grid( BC15, BC16, BC17, BC18, ncol=4)
family_16S <- plot_grid(title, family_16S_raw,ncol = 1,rel_heights = c(0.05, 1)) # rel_heights values control vertical title margins
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_prokaryote_family.png"), width = 12, height = 5  )



### microeukaryotes - 18S

## ASV LEVEL
# pick a year
year <- 2015

#load 18S microbial distance matrix ASV
bc_18S15_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_braycurtis.csv") )
# bc_18S15_ASV <- bc_18S15_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_18S15_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S15_ASV$Sites1 <- colnames(bc_18S15_ASV)
dist18S15.collapse <- bc_18S15_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S15.sites <- left_join( dist18S15.collapse, select(bc_18S15_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S15.sites <- left_join( dist18S15.sites, select(bc_18S15_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S15.sites <- dist18S15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S15.distance <- dist18S15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance.18S <- left_join(dist18S15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2015.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2016

#load 18S microbial distance matrix ASV
bc_18S16_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_braycurtis.csv") )
# bc_18S16_ASV <- bc_18S16_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_18S16_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S16_ASV$Sites1 <- colnames(bc_18S16_ASV)
dist18S16.collapse <- bc_18S16_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S16.sites <- left_join( dist18S16.collapse, select(bc_18S16_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S16.sites <- left_join( dist18S16.sites, select(bc_18S16_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S16.sites <- dist18S16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S16.distance <- dist18S16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance.18S <- left_join(dist18S16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2016.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC16 <- Hakai.2016.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2017

#load 18S microbial distance matrix ASV
bc_18S17_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_braycurtis.csv") )
# bc_18S17_ASV <- bc_18S17_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_18S17_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_18S17_ASV$Sites1 <- colnames(bc_18S17_ASV)
dist18S17.collapse <- bc_18S17_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S17.sites <- left_join( dist18S17.collapse, select(bc_18S17_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S17.sites <- left_join( dist18S17.sites, select(bc_18S17_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S17.sites <- dist18S17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S17.distance <- dist18S17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance.18S <- left_join(dist18S17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2017.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2018

#load 18S microbial distance matrix ASV
bc_18S18_ASV <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_braycurtis.csv") )
# bc_18S18_ASV <- bc_18S18_ASV %>% 
#   dplyr::rename("sample" = "X1")
bc_18S18_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/ASV_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S18_ASV$Sites1 <- colnames(bc_18S18_ASV)
dist18S18.collapse <- bc_18S18_ASV %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S18.sites <- left_join( dist18S18.collapse, select(bc_18S18_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S18.sites <- left_join( dist18S18.sites, select(bc_18S18_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S18.sites <- dist18S18.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S18.distance <- dist18S18.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2018.distance.18S <- left_join(dist18S18.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2018.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2018.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC18 <- Hakai.2018.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2018 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


windows(12,3)
title <-  ggdraw() + draw_label("ASV level",fontface = 'bold', size = 14, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
ASV_18S_raw <- cowplot::plot_grid( BC15, BC16, BC17, BC18, ncol=4)
ASV_18S <- plot_grid(title, plots,ncol = 1,rel_heights = c(0.05, 1)) # rel_heights values control vertical title margins
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_microeuk_ASV.png"), width = 12, height = 5  )





## GENUS LEVEL
# pick a year
year <- 2015

#load 18S microbial distance matrix GENUS
bc_18S15_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_braycurtis.csv") )
# bc_18S15_genus <- bc_18S15_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_18S15_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S15_genus$Sites1 <- colnames(bc_18S15_genus)
dist18S15.collapse <- bc_18S15_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S15.sites <- left_join( dist18S15.collapse, select(bc_18S15_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S15.sites <- left_join( dist18S15.sites, select(bc_18S15_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S15.sites <- dist18S15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S15.distance <- dist18S15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance.18S <- left_join(dist18S15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2015.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2016

#load 18S microbial distance matrix GENUS
bc_18S16_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_braycurtis.csv") )
# bc_18S16_genus <- bc_18S16_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_18S16_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S16_genus$Sites1 <- colnames(bc_18S16_genus)
dist18S16.collapse <- bc_18S16_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S16.sites <- left_join( dist18S16.collapse, select(bc_18S16_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S16.sites <- left_join( dist18S16.sites, select(bc_18S16_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S16.sites <- dist18S16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S16.distance <- dist18S16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance.18S <- left_join(dist18S16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2016.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC16 <- Hakai.2016.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2017

#load 18S microbial distance matrix GENUS
bc_18S17_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_braycurtis.csv") )
# bc_18S17_genus <- bc_18S17_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_18S17_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_18S17_genus$Sites1 <- colnames(bc_18S17_genus)
dist18S17.collapse <- bc_18S17_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S17.sites <- left_join( dist18S17.collapse, select(bc_18S17_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S17.sites <- left_join( dist18S17.sites, select(bc_18S17_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S17.sites <- dist18S17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S17.distance <- dist18S17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance.18S <- left_join(dist18S17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2017.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2018

#load 18S microbial distance matrix GENUS
bc_18S18_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_braycurtis.csv") )
# bc_18S18_genus <- bc_18S18_genus %>% 
#   dplyr::rename("sample" = "X1")
bc_18S18_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S18_genus$Sites1 <- colnames(bc_18S18_genus)
dist18S18.collapse <- bc_18S18_genus %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S18.sites <- left_join( dist18S18.collapse, select(bc_18S18_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S18.sites <- left_join( dist18S18.sites, select(bc_18S18_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S18.sites <- dist18S18.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S18.distance <- dist18S18.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2018.distance.18S <- left_join(dist18S18.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2018.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2018.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC18 <- Hakai.2018.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2018 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


windows(12,3)
title <-  ggdraw() + draw_label("genus level",fontface = 'bold', size = 14, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
genus_18S_raw  <- cowplot::plot_grid( BC15, BC16, BC17, BC18, ncol=4)
genus_18S <- plot_grid(title, genus_18S_raw,ncol = 1,rel_heights = c(0.05, 1)) # rel_heights values control vertical title margins
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_microeuk_genus.png"), width = 12, height = 5  )





### FAMILY LEVEL
# pick a year
year <- 2015

#load 18S microbial distance matrix family
bc_18S15_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_braycurtis.csv") )
# bc_18S15_family <- bc_18S15_family %>% 
#   dplyr::rename("sample" = "X1")
bc_18S15_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_metadata.csv") )



###Now make the distance matrices long
bc_18S15_family$Sites1 <- colnames(bc_18S15_family)
dist18S15.collapse <- bc_18S15_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S15.sites <- left_join( dist18S15.collapse, select(bc_18S15_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S15.sites <- left_join( dist18S15.sites, select(bc_18S15_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S15.sites <- dist18S15.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S15.distance <- dist18S15.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2015.distance.18S <- left_join(dist18S15.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2015.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2015.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2015 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC15 <- Hakai.2015.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2015 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2016

#load 18S microbial distance matrix family
bc_18S16_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_braycurtis.csv") )
# bc_18S16_family <- bc_18S16_family %>% 
#   dplyr::rename("sample" = "X1")
bc_18S16_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_18S16_family$Sites1 <- colnames(bc_18S16_family)
dist18S16.collapse <- bc_18S16_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S16.sites <- left_join( dist18S16.collapse, select(bc_18S16_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S16.sites <- left_join( dist18S16.sites, select(bc_18S16_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S16.sites <- dist18S16.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S16.distance <- dist18S16.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2016.distance.18S <- left_join(dist18S16.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2016.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2016.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2016 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC16 <- Hakai.2016.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2016 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2017

#load 18S microbial distance matrix family
bc_18S17_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_braycurtis.csv") )
# bc_18S17_family <- bc_18S17_family %>% 
#   dplyr::rename("sample" = "X1")
bc_18S17_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_18S17_family$Sites1 <- colnames(bc_18S17_family)
dist18S17.collapse <- bc_18S17_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S17.sites <- left_join( dist18S17.collapse, select(bc_18S17_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S17.sites <- left_join( dist18S17.sites, select(bc_18S17_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S17.sites <- dist18S17.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S17.distance <- dist18S17.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2017.distance.18S <- left_join(dist18S17.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2017.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2017.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2017 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC17 <- Hakai.2017.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2017 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


# pick a year
year <- 2018

#load 18S microbial distance matrix family
bc_18S18_family <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_braycurtis.csv") )
# bc_18S18_family <- bc_18S18_family %>% 
#   dplyr::rename("sample" = "X1")
bc_18S18_meta <- read_csv(paste0("R_Code_and_Analysis/mantel/family_18S_",year,"_metadata.csv") )


###Now make the distance matrices long
bc_18S18_family$Sites1 <- colnames(bc_18S18_family)
dist18S18.collapse <- bc_18S18_family %>% 
  gather(Sites2, Community_Distance, -Sites1)


# add sites from metadata
dist18S18.sites <- left_join( dist18S18.collapse, select(bc_18S18_meta, site, site_quadrat_id, labels), by=c("Sites1" = "labels") )
dist18S18.sites <- left_join( dist18S18.sites, select(bc_18S18_meta, site, site_quadrat_id, labels), by=c("Sites2" = "labels") )
dist18S18.sites <- dist18S18.sites %>% select( Sites1=site.x, Sites2=site.y, Community_Distance )

dist18S18.distance <- dist18S18.sites %>%
  # separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  # separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)


### Unite into one data frame
Hakai.2018.distance.18S <- left_join(dist18S18.distance,Hakai.geographic.distance, by = "Site.Pair") %>% 
  filter( !is.na(Community_Distance) )

# ### plots
# Graph1 <- Hakai.2018.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites1),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 microeukaryotes")+
#   geom_smooth(method = lm)
# Graph2 <- Hakai.2018.distance.18S %>%
#   drop_na(Geog_Distance) %>%
#   drop_na(Community_Distance)%>%
#   ggplot(aes(x = Geog_Distance, y = Community_Distance))+
#   theme_classic()+
#   geom_point(aes(colour = Sites2),alpha=0.25)+
#   xlab("Geographic distance (km)")+
#   ylab("B-C dissimilarity 2018 microeukaryotes")+
#   geom_smooth(method = lm)
# plot_grid(Graph1, Graph2, nrow = 2)


BC18 <- Hakai.2018.distance.18S %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  theme_classic()+
  geom_point(alpha=0.25)+
  xlab("Geographic distance (km)")+
  ylab("B-C dissimilarity\n2018 microeukaryotes")+
  ylim( c(0,1) ) + xlim( c(0,41) ) +
  geom_smooth(method = lm)


windows(12,3)
title <-  ggdraw() + draw_label("family level",fontface = 'bold', size = 14, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
family_18S_raw <- cowplot::plot_grid( BC15, BC16, BC17, BC18, ncol=4)
family_18S <- plot_grid(title, family_18S_raw,ncol = 1,rel_heights = c(0.05, 1)) # rel_heights values control vertical title margins
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_microeuk_family.png"), width = 12, height = 5  )


######################################
### saving all at the finest level ###
######################################
title_ASV <-  ggdraw() + draw_label("Finest taxonomic level (ASV for microbes)",fontface = 'bold', size = 18, x = 0.35, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
finest_ASV <- cowplot::plot_grid (ASV_16S_raw, ASV_18S_raw,macro, ncol=1)
finest_ASV_title <- plot_grid(title_ASV, finest_ASV, ncol = 1,rel_heights = c(0.05, 1))
finest_ASV_title
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_all_finest_ASV.png"), width = 12, height = 10  )

title_genus <-  ggdraw() + draw_label("Finest taxonomic level (genus for microbes)",fontface = 'bold', size = 18, x = 0.35, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
finest_genus <- cowplot::plot_grid (ASV_16S_raw, ASV_18S_raw,macro, ncol=1)
finest_genus_title <- plot_grid(title_genus, finest_genus, ncol = 1,rel_heights = c(0.05, 1))
finest_genus_title
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_all_finest_genus.png"), width = 12, height = 10  )

########################################
### saving all at the coarsest level ###
########################################
title_family <-  ggdraw() + draw_label("Family level",fontface = 'bold', size = 18, x = 0.5, hjust = 0) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
family <- cowplot::plot_grid (family_16S_raw, family_18S_raw,macro, ncol=1)
family_title <- plot_grid(title_family, family, ncol = 1,rel_heights = c(0.05, 1))
family_title
ggsave(paste0("R_Code_and_Analysis/distance_decay/BCdecay_all_family.png"), width = 12, height = 10  )
