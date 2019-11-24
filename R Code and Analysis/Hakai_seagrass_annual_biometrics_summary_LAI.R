###################################################################
###  Hakai + O'Connor Seagrass Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will summarize biometric data on the quadrat level
###       to calculate leaf surface area, average across all shoots,  
###       and then calculate an estimate of  Leaf Area Index 
###  Leaf Area Index (LAI) = total eelgrass surface area per square meter
###  This is a complementary aspect of living space to eelgrass biomass
### 
###  code by Matt Whalen + Keila Stark
###  started on   23 November 2019
###  
###################################################################


## libraries
library(tidyverse)


# what to do with 2015? Different number of shoots, not included in quadrat level dataset
#   can maybe use single shoot data here??

## read data
qall <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_quadrats.csv" )


## calculate surface area per shoot
q.sa.mean <- qall %>% 
  mutate( SA = shoot.length*shoot.width ) %>% 
  # average across replicate shoots
  group_by(site, id) %>% 
  summarize_all( mean, na.rm=T )




