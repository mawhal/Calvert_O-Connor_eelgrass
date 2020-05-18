### 2017 Microbial abundance data (18S) to use in MANTEL analysis ###
### Author: Bianca Trevizan Segovia ###
### Date created: January 14, 2020 ###
### Date modified: February 14, 2020 ###
### modified by Whalen on 18 May 2020: save metadata along with Bray-Curtis distance matrix

library(dplyr)
library(vegan)
library(usedist)


allsamples_18S <- read_csv("R_Code_and_Analysis/mantel/genus_level_18S_MASTER.csv")

allsamples_18S <- allsamples_18S[order(allsamples_18S$year, allsamples_18S$site),]

### filter only meso quadrat surveys
meso_quadrat <- c("meso_quadrat")
zostera_meso_quadrat <- allsamples_18S %>% 
  dplyr::filter(survey_type %in% meso_quadrat)

### filter only zostera leaf_old samples
leaf_old <- c("leaf_old")
zostera_old <- zostera_meso_quadrat %>% 
  dplyr::filter(sample_type %in% leaf_old)

zostera_old  <- zostera_old  %>% 
  dplyr::mutate(site=recode(site,
                            "choked_south_pigu" = "choked_inner",
                            "goose_southeast" = "goose_south_east",
                            "goose_southwest" = "goose_south_west",
                            "mcmullin_north" = "mcmullins_north",
                            "mcmullin_south" = "mcmullins_south",
                            "pruth_bay_south" = "pruth_bay"))

###################################
############### 2015 ##############
###################################

### filter only 2015 samples
only_2015 <- c("2015")
zostera_old_2015 <- zostera_old %>% 
  dplyr::filter(year %in% only_2015)
zostera_old_2015$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_18S_2015  <- zostera_old_2015 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,meso_quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_18S_2015 <- microbes_18S_2015 %>%
  select(site_quadrat_id, everything())

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_18S_2015$labels <- vegan::make.cepnames(microbes_18S_2015$site_quadrat_id)

### move labels from last column to first
microbes_18S_2015 <- microbes_18S_2015 %>%
  select(labels, everything())

names(microbes_18S_2015)
#abundance only
microbes_18S_2015_abund <- microbes_18S_2015 %>%
  select(c(11:ncol(microbes_18S_2015)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_18S_2015 <- vegdist(microbes_18S_2015_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_18S_2015 <- as.matrix(dist_microbes_18S_2015)
# mat_microbes_18S_2015[upper.tri(mat_microbes_18S_2015, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_18S_2015$labels
rownames(mat_microbes_18S_2015) <- site_id
colnames(mat_microbes_18S_2015) <- site_id
write_csv( data.frame(mat_microbes_18S_2015), "R_Code_and_Analysis/mantel/genus_18S_2015_braycurtis.csv" )

# save metadata
meta_18S_2015 <- microbes_18S_2015 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_18S_2015, "R_Code_and_Analysis/mantel/genus_18S_2015_metadata.csv" )

###################################
############### 2016 ##############
###################################

### filter only 2016 samples
only_2016 <- c("2016")
zostera_old_2016 <- zostera_old %>% 
  dplyr::filter(year %in% only_2016)
zostera_old_2016$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_18S_2016  <- zostera_old_2016 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,meso_quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_18S_2016 <- microbes_18S_2016 %>%
  select(site_quadrat_id, everything())

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_18S_2016$labels <- vegan::make.cepnames(microbes_18S_2016$site_quadrat_id)

### move labels from last column to first
microbes_18S_2016 <- microbes_18S_2016 %>%
  select(labels, everything())

names(microbes_18S_2016)
#abundance only
microbes_18S_2016_abund <- microbes_18S_2016 %>%
  select(c(11:ncol(microbes_18S_2016)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_18S_2016 <- vegdist(microbes_18S_2016_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_18S_2016 <- as.matrix(dist_microbes_18S_2016)
# mat_microbes_18S_2016[upper.tri(mat_microbes_18S_2016, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_18S_2016$labels
rownames(mat_microbes_18S_2016) <- site_id
colnames(mat_microbes_18S_2016) <- site_id
write_csv( data.frame(mat_microbes_18S_2016), "R_Code_and_Analysis/mantel/genus_18S_2016_braycurtis.csv" )

# save metadata
meta_18S_2016 <- microbes_18S_2016 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_18S_2016, "R_Code_and_Analysis/mantel/genus_18S_2016_metadata.csv" )


###################################
############### 2017 ##############
###################################

### filter only 2017 samples
only_2017 <- c("2017")
zostera_old_2017 <- zostera_old %>% 
  dplyr::filter(year %in% only_2017)
zostera_old_2017$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_18S_2017  <- zostera_old_2017 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,meso_quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_18S_2017 <- microbes_18S_2017 %>%
  select(site_quadrat_id, everything())

## drop NA rows *** samples that did not work
microbes_18S_2017 <- microbes_18S_2017 %>% drop_na()

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_18S_2017$labels <- vegan::make.cepnames(microbes_18S_2017$site_quadrat_id)

### move labels from last column to first
microbes_18S_2017 <- microbes_18S_2017 %>%
  select(labels, everything())

names(microbes_18S_2017)
#abundance only
microbes_18S_2017_abund <- microbes_18S_2017 %>%
  select(c(11:ncol(microbes_18S_2017)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_18S_2017 <- vegdist(microbes_18S_2017_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_18S_2017 <- as.matrix(dist_microbes_18S_2017)
# mat_microbes_18S_2017[upper.tri(mat_microbes_18S_2017, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_18S_2017$labels
rownames(mat_microbes_18S_2017) <- site_id
colnames(mat_microbes_18S_2017) <- site_id
write_csv( data.frame(mat_microbes_18S_2017), "R_Code_and_Analysis/mantel/genus_18S_2017_braycurtis.csv" )

# save metadata
meta_18S_2017 <- microbes_18S_2017 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_18S_2017, "R_Code_and_Analysis/mantel/genus_18S_2017_metadata.csv" )


###################################
############### 2018 ##############
###################################

### filter only 2018 samples
only_2018 <- c("2018")
zostera_old_2018 <- zostera_old %>% 
  dplyr::filter(year %in% only_2018)
zostera_old_2018$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_18S_2018  <- zostera_old_2018 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,meso_quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_18S_2018 <- microbes_18S_2018 %>%
  select(site_quadrat_id, everything())

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_18S_2018$labels <- vegan::make.cepnames(microbes_18S_2018$site_quadrat_id)

### move labels from last column to first
microbes_18S_2018 <- microbes_18S_2018 %>%
  select(labels, everything())

names(microbes_18S_2018)
#abundance only
microbes_18S_2018_abund <- microbes_18S_2018 %>%
  select(c(11:ncol(microbes_18S_2018)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_18S_2018 <- vegdist(microbes_18S_2018_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_18S_2018 <- as.matrix(dist_microbes_18S_2018)
# mat_microbes_18S_2018[upper.tri(mat_microbes_18S_2018, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_18S_2018$labels
rownames(mat_microbes_18S_2018) <- site_id
colnames(mat_microbes_18S_2018) <- site_id
write_csv( data.frame(mat_microbes_18S_2018), "R_Code_and_Analysis/mantel/genus_18S_2018_braycurtis.csv" )

# save metadata
meta_18S_2018 <- microbes_18S_2018 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_18S_2018, "R_Code_and_Analysis/mantel/genus_18S_2018_metadata.csv" )

