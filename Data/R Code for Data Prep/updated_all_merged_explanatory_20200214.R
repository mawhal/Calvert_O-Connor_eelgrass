### Merging all seagrass metrics, abiotic data and coordinates to use in HMSC analysis ###
### Author: Bianca Trevizan Segovia ###
### Date created: February 04, 2020 ###

library(dplyr)
library(vegan)

########################                
##### ABIOTIC DATA #####
########################    

# load data
abiotic <- read.csv("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R Code and Analysis/output data/Bia_reviewed_O'Connor_hakai_seagrass_MASTER_abiotic.csv")

# abiotic variables of interest = temperature, dissolved_oxygen_concentration, salinity, pH
abiotic_variables <- abiotic %>% 
  dplyr::select(date, site, year, depth, salinity, temperature, ph)

#############################                
##### SPATIAL LAT LONG ######
#############################
coordinates <- read.csv("~/PostDoc/projects/Calvert_O-Connor_eelgrass/metadata/coordinates.csv", header=T)

#################################                
##### SEAGRASS METRICS DATA #####
#################################    

##### QUADRAT LEVEL ####
# Emily's summarised seagrass metrics for all years with LAI
eelgrass_metrics <- read.csv("~/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/data_oconnor/seagrass_emily_files_final_code/hakai_quad_combined_20200214.csv", header=T)

############################################################
##### JOIN ABIOTICS, COORDINATES AND  SEAGRASS METRICS DATA ##### 
############################################################

# inner_join = retain only rows that matches
# because abiotics is already reflecting only species data, I'll use that
environmental<- full_join(abiotic_variables, eelgrass_metrics, by = c("site", "year"))

environmental_spatial <- full_join(coordinates, environmental, by = c("site"))

### concatenate site and quadrat_id to create quadrat for random effect
environmental_spatial  <- environmental_spatial  %>%
  dplyr::mutate(quadrat = paste(site,quadrat_id, sep="_"))

### concatenate quadrat and year to create quadrat_year for when running the model with all years
environmental_spatial  <- environmental_spatial  %>%
  dplyr::mutate(quadrat_year = paste(quadrat,year, sep="_"))

environmental_spatial  <- environmental_spatial[order(environmental_spatial$year),] 

write.csv(environmental_spatial, file = "~/PostDoc/projects/Calvert_O-Connor_eelgrass/R Code and Analysis/output data/merged_explanatory.csv", quote = F, row.names = F)

