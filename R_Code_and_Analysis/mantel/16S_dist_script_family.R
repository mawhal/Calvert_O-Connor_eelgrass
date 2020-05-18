### 2017 Microbial abundance data (16S) to use in MANTEL analysis ###
### Author: Bianca Trevizan Segovia ###
### Date created: January 13, 2020 ###
### Date modified: February 24, 2020 ### new file 16S 2016 not rarafied
### modified by Whalen on 18 May 2020: save metadata along with Bray-Curtis distance matrix

# library(dplyr)
library(tidyverse)
library(vegan)


allsamples <- read_csv("R_Code_and_Analysis/mantel/family_level_16S_MASTER_REMOVED_CONT_NOT_RAREFIED.csv")

allsamples <- allsamples[order(allsamples$year, allsamples$site),]

### filter only meso quadrat surveys
meso_quadrat <- c("meso_quadrat")
zostera_meso_quadrat <- allsamples %>% 
  dplyr::filter(survey_type %in% meso_quadrat)

### filter only zostera leaf_old samples
leaf_old <- c("leaf_old")
zostera_old <- zostera_meso_quadrat %>% 
  dplyr::filter(sample_type %in% leaf_old)
zostera_old$sample_type

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
microbes_16S_2015  <- zostera_old_2015 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_16S_2015 <- microbes_16S_2015 %>%
  dplyr::select(site_quadrat_id, everything())

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_16S_2015$labels <- vegan::make.cepnames(microbes_16S_2015$site_quadrat_id)

### move labels from last column to first
microbes_16S_2015 <- microbes_16S_2015 %>%
  dplyr::select(labels, everything())

names(microbes_16S_2015)
#abundance only
microbes_16S_2015_abund <- microbes_16S_2015 %>%
  dplyr::select(c(16:ncol(microbes_16S_2015)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_16S_2015 <- vegdist(microbes_16S_2015_abund, method = "bray")


# set new matrix from dist matrix discarding upper triangle
mat_microbes_16S_2015 <- as.matrix(dist_microbes_16S_2015)
# mat_microbes_16S_2015[upper.tri(mat_microbes_16S_2015, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_16S_2015$labels
rownames(mat_microbes_16S_2015) <- site_id
colnames(mat_microbes_16S_2015) <- site_id
write_csv( data.frame(mat_microbes_16S_2015), "R_Code_and_Analysis/mantel/family_16S_2015_braycurtis.csv" )

# save metadata
meta_16S_2015 <- microbes_16S_2015 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_16S_2015, "R_Code_and_Analysis/mantel/family_16S_2015_metadata.csv" )



###################################
############### 2016 ##############
###################################

### ********** contamination, so a lot of samples did not work ********** ###

### filter only 2016 samples
only_2016 <- c("2016")
zostera_old_2016 <- zostera_old %>% 
  dplyr::filter(year %in% only_2016)
zostera_old_2016$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_16S_2016  <- zostera_old_2016 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_16S_2016 <- microbes_16S_2016 %>%
  dplyr::select(site_quadrat_id, everything())

## drop NA rows *** samples that did not work
microbes_16S_2016 <- microbes_16S_2016 %>% drop_na()

### remove mcmullins 2016
remove_mcmullins <- c("mcmullins")
microbes_16S_2016 <- microbes_16S_2016 %>% 
  dplyr::filter(!region %in% remove_mcmullins )

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_16S_2016$labels <- vegan::make.cepnames(microbes_16S_2016$site_quadrat_id)

### move labels from last column to first
microbes_16S_2016 <- microbes_16S_2016 %>%
  dplyr::select(labels, everything())

names(microbes_16S_2016)
#abundance only
microbes_16S_2016_abund <- microbes_16S_2016 %>%
  dplyr::select(c(16:ncol(microbes_16S_2016)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_16S_2016 <- vegdist(microbes_16S_2016_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_16S_2016 <- as.matrix(dist_microbes_16S_2016)
# mat_microbes_16S_2016[upper.tri(mat_microbes_16S_2016, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_16S_2016$labels
rownames(mat_microbes_16S_2016) <- site_id
colnames(mat_microbes_16S_2016) <- site_id
write_csv( data.frame(mat_microbes_16S_2016), "R_Code_and_Analysis/mantel/family_16S_2016_braycurtis.csv" )

# save metadata
meta_16S_2016 <- microbes_16S_2016 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_16S_2016, "R_Code_and_Analysis/mantel/family_16S_2016_metadata.csv" )

###################################
############### 2017 ##############
###################################

### filter only 2017 samples
only_2017 <- c("2017")
zostera_old_2017 <- zostera_old %>% 
  dplyr::filter(year %in% only_2017)
zostera_old_2017$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_16S_2017  <- zostera_old_2017 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_16S_2017 <- microbes_16S_2017 %>%
  dplyr::select(site_quadrat_id, everything())

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_16S_2017$labels <- vegan::make.cepnames(microbes_16S_2017$site_quadrat_id)

### move labels from last column to first
microbes_16S_2017 <- microbes_16S_2017 %>%
  dplyr::select(labels, everything())

names(microbes_16S_2017)
#abundance only
microbes_16S_2017_abund <- microbes_16S_2017 %>%
  dplyr::select(c(16:ncol(microbes_16S_2017)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_16S_2017 <- vegdist(microbes_16S_2017_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_16S_2017 <- as.matrix(dist_microbes_16S_2017)
# mat_microbes_16S_2017[upper.tri(mat_microbes_16S_2017, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_16S_2017$labels
rownames(mat_microbes_16S_2017) <- site_id
colnames(mat_microbes_16S_2017) <- site_id
write_csv( data.frame(mat_microbes_16S_2017), "R_Code_and_Analysis/mantel/family_16S_2017_braycurtis.csv" )

# save metadata
meta_16S_2017 <- microbes_16S_2017 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_16S_2017, "R_Code_and_Analysis/mantel/family_16S_2017_metadata.csv" )

###################################
############### 2018 ##############
###################################

### filter only 2018 samples
only_2018 <- c("2018")
zostera_old_2018 <- zostera_old %>% 
  dplyr::filter(year %in% only_2018)
zostera_old_2018$year

### concatenate site and quadrat_id to create spatial quadrat column
microbes_16S_2018  <- zostera_old_2018 %>% 
  dplyr::mutate(site_quadrat_id = paste(site,quadrat_id, sep="_")) 

### move site_quadrat_id from last column to first
microbes_16S_2018 <- microbes_16S_2018 %>%
  dplyr::select(site_quadrat_id, everything())

## drop NA rows *** samples that did not work
microbes_16S_2018 <- microbes_16S_2018 %>% drop_na()

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
microbes_16S_2018$labels <- vegan::make.cepnames(microbes_16S_2018$site_quadrat_id)

### move labels from last column to first
microbes_16S_2018 <- microbes_16S_2018 %>%
  dplyr::select(labels, everything())

names(microbes_16S_2018)
#abundance only
microbes_16S_2018_abund <- microbes_16S_2018 %>%
  dplyr::select(c(16:ncol(microbes_16S_2018)))

#abundance data frame - bray curtis dissimilarity
dist_microbes_16S_2018 <- vegdist(microbes_16S_2018_abund, method = "bray")

# set new matrix from dist matrix discarding upper triangle
mat_microbes_16S_2018 <- as.matrix(dist_microbes_16S_2018)
# mat_microbes_16S_2018[upper.tri(mat_microbes_16S_2018, diag = TRUE)] <- NA
# add labels 
site_id <- microbes_16S_2018$labels
rownames(mat_microbes_16S_2018) <- site_id
colnames(mat_microbes_16S_2018) <- site_id
write_csv( data.frame(mat_microbes_16S_2018), "R_Code_and_Analysis/mantel/family_16S_2018_braycurtis.csv" )

# save metadata
meta_16S_2018 <- microbes_16S_2018 %>% select(year,site,site_quadrat_id,labels)
write_csv( meta_16S_2018, "R_Code_and_Analysis/mantel/family_16S_2018_metadata.csv" )
