################Preliminary Analyses to explore the impact of distance between sites on grazer dissimilarity
##Started by Coreen April 2020
##This script  makes distance matrices between sites and plots Bray-Curtis grazer dissimilarity between each pair of sites against dstance between each pair of sites for 2014-207
##For some reason i did each year totally separately. No idea why

library(vegan)
library(tidyverse)
library(distances)
library(cowplot)


# Geographic distance between site pairs ---------------------------------------

Hakaispatial <- read.csv("metadata/00_Hakai_UBC_metadata_MASTER - geolocation_site.csv")

Hakaispatial1 <- Hakaispatial %>%
  select(site_name,lat, long)

##Change site names that don't match master grazer data

Hakaispatial1$site_name <- recode(Hakaispatial1$site_name, "inner_choked" = "choked_inner", "sandspit" = "choked_sandspit")

## Make distance matrix

Hakai.distance<- distances(Hakaispatial1, id_variable = "site_name", dist_variables = NULL)

###Now make the distance matricies long

Hakai.geog <-as.data.frame(as.matrix(Hakai.distance))

Hakai.geog$Sites1 <- rownames(Hakai.geog)

##Data frame of distances between all site pairs

Hakai.geographic.distance <- Hakai.geog %>% 
  gather(Sites2, Geog_Distance, - Sites1)%>%
  unite("Site.Pair", Sites1, Sites2, sep = "-", remove = FALSE)



# 2016 grazers ------------------------------------------------------------

grazers_all <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")

###Spread the data for an nmds

grazers2016 <- grazers_all %>%
  filter(year == "2016")%>%
  select(-X, -date, -size)%>%
  group_by(site, sample)%>%
  add_count(taxon, name = "taxon.count")%>%
  distinct()

grazers2016.for.nmds <- grazers2016%>%
  group_by(site, sample)%>%
  spread(taxon, taxon.count, fill = 0)

###Now do the nmds

#Remove identifier columns
sppfor.mds <- grazers2016.for.nmds %>%
  ungroup(site, sample)%>%
  select(-site, -sample, -year)

spp.collapse <- grazers2016.for.nmds
###Some columns and rows have nothing in them- remove
###now do the mds

spp.mds<-metaMDS(sppfor.mds, trymax = 100)

stressplot(spp.mds)

## make distance variables ordered


##make the site column groupings
data.scores <- as.data.frame(scores(spp.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.fra

grpsite <-spp.collapse$site  #  add the grp variable created earlier
grp.sample <-spp.collapse$sample
data.scores$site <- grpsite
data.scores$sample <- grp.sample

##make the species column groupings
species.scores <- as.data.frame(scores(spp.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

##### Now make community distance matrix between points based on nmds

###Plot the nmds

data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  #facet_wrap(vars(res_area))+
  #geom_point(aes(size = 2))+
  theme_classic()+
  title("2016 Grazers")+
  geom_point(aes(colour = site, size = .5))+
  xlab("nmds1")+
  ylab("nmds2")

data.scores1 <- data.scores %>%
  unite("Site.Sample", site, sample, sep = "-", remove = TRUE)

Hakai.community.distance<- distances(data.scores1, id_variable = "Site.Sample", dist_variables = NULL)


###Now make the distance matricies long

Hakai.community <-as.data.frame(as.matrix(Hakai.community.distance))

Hakai.community$Sites1 <- rownames(Hakai.community)


Hakai.community.collapse <- Hakai.community %>% 
  gather(Sites2, Community_Distance, - Sites1)

Hakai.species.distance <- Hakai.community.collapse%>%
  separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Site_1, Site_2, sep = "-", remove = FALSE)

### Unite into one data frame

Hakai.2016.distance <- left_join(Hakai.geographic.distance, Hakai.species.distance, by = "Site.Pair")

###plot it- I don't know of a good way to show the sites pairs because each pair yields way too many colours so I have a plot of each of site 1 and 2 of each pair

Graph1 <- Hakai.2016.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  geom_jitter()+
  theme_classic()+
  geom_jitter(aes(colour = Site_1))+
  xlab("Relative Distance")+
  ylab("B-C Dissimilarity 2016 Grazers")+
  geom_smooth(method = lm)

Graph2 <- Hakai.2016.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  #filter(age.class == "geriatric")%>%
  #filter(conception.yday > 100) %>%
  #group_by(year)%>%
  #drop_na(age.class)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  #facet_wrap(vars(res_area))+
  geom_jitter()+
  theme_classic()+
  geom_jitter(aes(colour = Site_2))+
  xlab("Relative Distance")+
  ylab("B-C Dissimilarity 2016 Grazers")+
  geom_smooth(method = lm)

plot_grid(Graph1, Graph2, nrow = 2)



# 2017 Grazers --------------------------------------------------------------------
grazers_all <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")


###Spread the data for an nmds
grazers2017 <- grazers_all %>%
  filter(year == "2017")%>%
  select(-X, -date, -size)%>%
  group_by(site, sample)%>%
  add_count(taxon, name = "taxon.count")%>%
  distinct()

grazers2017.for.nmds <- grazers2017%>%
  group_by(site, sample)%>%
  spread(taxon, taxon.count, fill = 0)

###Now do the nmds

#Remove identifier columns
sppfor.mds <- grazers2017.for.nmds %>%
  ungroup(site, sample)%>%
  select(-site, -sample, -year)

spp.collapse <- grazers2017.for.nmds
###Some columns and rows have nothing in them- remove
###now do the mds

spp.mds<-metaMDS(sppfor.mds, trymax = 100)

stressplot(spp.mds)


############## make distance variables ordered


##make the site column groupings
data.scores <- as.data.frame(scores(spp.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.fra

grpsite <-spp.collapse$site  #  add the grp variable created earlier
grp.sample <-spp.collapse$sample
data.scores$site <- grpsite
data.scores$sample <- grp.sample
#head(data.scores)
##make the species column groupings
species.scores <- as.data.frame(scores(spp.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
#head(species.scores)  #look at the data

# 2015 Grazers --------------------------------------------------------------------


grazers_all <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")


###Spread the data for an nmds

grazers2015 <- grazers_all %>%
  filter(year == "2015")%>%
  select(-X, -date, -size)%>%
  group_by(site, sample)%>%
  add_count(taxon, name = "taxon.count")%>%
  distinct()


grazers2015.for.nmds <- grazers2015%>%
  group_by(site, sample)%>%
  spread(taxon, taxon.count, fill = 0)

###Now do the nmds

#Remove identifier columns
sppfor.mds <- grazers2015.for.nmds %>%
  ungroup(site, sample)%>%
  select(-site, -sample, -year)

spp.collapse <- grazers2015.for.nmds
###Some columns and rows have nothing in them- remove
###now do the mds

spp.mds<-metaMDS(sppfor.mds, trymax = 200)

stressplot(spp.mds)


############## make distance variables ordered


##make the site column groupings
data.scores <- as.data.frame(scores(spp.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.fra

grpsite <-spp.collapse$site  #  add the grp variable created earlier
grp.sample <-spp.collapse$sample
data.scores$site <- grpsite
data.scores$sample <- grp.sample

##make the species column groupings
species.scores <- as.data.frame(scores(spp.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

##### Now make community distance matrix between points based on nmds

###Plot the nmds


data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  #facet_wrap(vars(res_area))+
  #geom_point(aes(size = 2))+
  theme_classic()+
  title("2015 Grazers")+
  geom_point(aes(colour = site, size = .5))+
  xlab("nmds1")+
  ylab("nmds2")

data.scores1 <- data.scores %>%
  unite("Site.Sample", site, sample, sep = "-", remove = TRUE)

Hakai.community.distance<- distances(data.scores1, id_variable = "Site.Sample", dist_variables = NULL)

print(Hakai.community.distance)

###Now make the distance matricies long

Hakai.community <-as.data.frame(as.matrix(Hakai.community.distance))

Hakai.community$Sites1 <- rownames(Hakai.community)


Hakai.community.collapse <- Hakai.community %>% 
  gather(Sites2, Community_Distance, - Sites1)

Hakai.species.distance <- Hakai.community.collapse%>%
  separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Site_1, Site_2, sep = "-", remove = FALSE)

### Unite into one data frame

Hakai.2015.distance <- left_join(Hakai.geographic.distance, Hakai.species.distance, by = "Site.Pair")

######plot it- I don't know of a good way to show the sites pairs because each pair yields way too many colours so I have a plot of each of site 1 and 2 of each pair

library(cowplot)


Graph1 <- Hakai.2015.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  geom_jitter()+
  theme_classic()+
  geom_jitter(aes(colour = Site_1))+
  xlab("Relative Distance")+
  ylab("B-C Dissimilarity 2015 Grazers")+
  geom_smooth(method = lm)

Graph2 <- Hakai.2015.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  geom_jitter()+
  theme_classic()+
  geom_jitter(aes(colour = Site_2))+
  xlab("Relative Distance")+
  ylab("B-C Dissimilarity 2015 Grazers")+
  geom_smooth(method = lm)

plot_grid(Graph1, Graph2, nrow = 2)


# 2014 Grazers --------------------------------------------------------------------

grazers_all <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")


###Spread the data for an nmds

grazers2014 <- grazers_all %>%
  filter(year == "2014")%>%
  select(-X, -date, -size)%>%
  group_by(site, sample)%>%
  add_count(taxon, name = "taxon.count")%>%
  distinct()


grazers2014.for.nmds <- grazers2014%>%
  group_by(site, sample)%>%
  spread(taxon, taxon.count, fill = 0)

###Now do the nmds

#Remove identifier columns
sppfor.mds <- grazers2014.for.nmds %>%
  ungroup(site, sample)%>%
  select(-site, -sample, -year)

spp.collapse <- grazers2014.for.nmds
###Some columns and rows have nothing in them- remove
###now do the mds

spp.mds<-metaMDS(sppfor.mds, trymax = 200)

stressplot(spp.mds)


############## make distance variables ordered


##make the site column groupings
data.scores <- as.data.frame(scores(spp.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.fra
grpsite <-spp.collapse$site  #  add the grp variable created earlier
grp.sample <-spp.collapse$sample
data.scores$site <- grpsite
data.scores$sample <- grp.sample
#head(data.scores)
##make the species column groupings
species.scores <- as.data.frame(scores(spp.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

##### Now make community distance matrix between points based on nmds

###Plot the nmds


data.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  theme_classic()+
  title("2014 Grazers")+
  geom_point(aes(colour = site, size = .5))+
  xlab("nmds1")+
  ylab("nmds2")

data.scores1 <- data.scores %>%
  unite("Site.Sample", site, sample, sep = "-", remove = TRUE)

Hakai.community.distance<- distances(data.scores1, id_variable = "Site.Sample", dist_variables = NULL)

print(Hakai.community.distance)

###Now make the distance matricies long

Hakai.community <-as.data.frame(as.matrix(Hakai.community.distance))

Hakai.community$Sites1 <- rownames(Hakai.community)


Hakai.community.collapse <- Hakai.community %>% 
  gather(Sites2, Community_Distance, - Sites1)

Hakai.species.distance <- Hakai.community.collapse%>%
  separate(Sites1, c("Site_1", "Sample_1"), sep = "-", remove = TRUE)%>%
  separate(Sites2, c("Site_2", "Sample_2"), sep = "-", remove = TRUE)%>%
  unite("Site.Pair", Site_1, Site_2, sep = "-", remove = FALSE)

### Unite into one data frame

Hakai.2014.distance <- left_join(Hakai.geographic.distance, Hakai.species.distance, by = "Site.Pair")

###plot it- I don't know of a good way to show the sites pairs because each pair yields way too many colours so I have a plot of each of site 1 and 2 of each pair

library(cowplot)


Graph1 <- Hakai.2014.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  geom_jitter()+
  theme_classic()+
  geom_jitter(aes(colour = Site_1))+
  xlab("Relative Distance")+
  ylab("B-C Dissimilarity 2014 Grazers")+
  geom_smooth(method = lm)

Graph2 <- Hakai.2014.distance %>%
  drop_na(Geog_Distance) %>%
  drop_na(Community_Distance)%>%
  ggplot(aes(x = Geog_Distance, y = Community_Distance))+
  geom_jitter()+
  theme_classic()+
  geom_jitter(aes(colour = Site_2))+
  xlab("Relative Distance")+
  ylab("B-C Dissimilarity 2014 Grazers")+
  geom_smooth(method = lm)

plot_grid(Graph1, Graph2, nrow = 2)