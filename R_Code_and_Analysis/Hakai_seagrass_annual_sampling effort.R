####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code show the sampling effort for each year
###     included: number of sites, quadrats, shoots measured, grazer
###               samples sorted
###  code by Matt Whalen
###  started on   21 Novemeber 2019
###  
####################################################################



# load libraries
library( tidyverse )


# read data
# mesograzers
m <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# eelgrass biometrics
q <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_quadrats.csv" )
s <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_shoots.csv" )
# abiotic
a <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_abiotic.csv" )



# summaries
# number of grazer samples sorted
msum <- m %>% 
  group_by( year, site ) %>% 
  summarize( quads=length(unique(sample)) )

msum %>% 
  spread( site, quads, fill=0 )

# number of quads
q %>% 
  group_by( year, site ) %>% 
  summarize( quads=length(unique(id)))

# number of shoots measured per quad
q %>% 
  group_by( year, site, id ) %>% 
  summarize( shoots=length(unique(sheath.length)) )

# number of shoos for epiphytes
s %>% 
  group_by( year, site, id ) %>% 
  summarize(  shoots=length(unique(shoot.length)) )
