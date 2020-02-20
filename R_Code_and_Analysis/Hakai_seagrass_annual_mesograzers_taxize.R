####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will collect taxonomic information on all epifaunal 
###      taxa in the dataset
###  code by Matt Whalen
###  started on   12 December 2019
###  
####################################################################


## load libraries
library(tidyverse)
library(taxize)

# API key for taxize package
ENTREZ_KEY = "a06da45de96c4b7d680015d4c5b46694d908" 


# read data
taxa.raw <- read_csv("output data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv")
# remove rows that we can't use
taxa.rem <- taxa.raw$taxon2[ is.na(taxa.raw$remove) ]
# replace periods with spaces
taxa.space <- gsub( "[.]"," ", taxa.rem )
# delete sp designations
taxa.sp <- gsub( "sp","", taxa.space )
# remove trailing white space
taxa.white <- trimws(taxa.sp)
# remove NAs
taxa.final <- taxa.white[ !is.na(taxa.white) ]
# unique entries only
taxa.clean <- sort(unique(taxa.final))

# get taxonomy for all taxa
taxa.classify <- classification( taxa.clean, db="worms")#, ENTREZ_KEY = "a06da45de96c4b7d680015d4c5b46694d908"  )
# rbind all entries
taxa.bind <- rbind(taxa.classify)

# spread out the taxonomic levels as columns
taxa.spread <- taxa.bind %>% 
  select( -id ) %>% 
  spread( rank,name )

# write the taxonmy to disk
write_csv( taxa.spread, "output data/O'Connor_hakai_seagrass_epifauna_taxonomy.csv")


