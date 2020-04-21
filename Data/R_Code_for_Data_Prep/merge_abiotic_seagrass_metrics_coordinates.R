### Merging all seagrass metrics, abiotic data and coordinates to use in HMSC analysis ###
### Author: Bianca Trevizan Segovia ###
### Date created: February 04, 2020 ###

library(dplyr)
library(vegan)

########################                
##### ABIOTIC DATA #####
########################    

# load data
abiotic <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_abiotic_20200214.csv")

# abiotic variables of interest = temperature, dissolved_oxygen_concentration, salinity, pH
abiotic_variables <- abiotic %>% 
  dplyr::select(date, site, year, depth, salinity, temperature, ph)

#############################                
##### SPATIAL LAT LONG ######
#############################
coordinates <- read.csv("metadata/coordinates.csv", header=T)

#################################                
##### SEAGRASS METRICS DATA #####
#################################    

##### QUADRAT LEVEL ####
# Emily's summarised seagrass metrics for all years with LAI
eelgrass_metrics <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_seagrass_metrics_20200214.csv", header=T)

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

write.csv(environmental_spatial, file = "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", quote = F, row.names = F)

