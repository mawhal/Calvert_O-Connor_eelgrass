####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will prepare O'Connor data for HMSC test
###  code by Matt Whalen
###  started on   23 Novemeber 2019
###  
####################################################################



# load libraries
library( tidyverse )


# read data
m <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan

### select sites
# filter taxa and sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon4))
# # Four core sites
# mfilt <- m %>%
#   filter( is.na(remove), !is.na(taxon4),
#           site %in% c("inner choked","sandspit",
#                       "triquet north","triquet south") )


## 2017 data to start
m17 <- mfilt %>% 
  filter( year==2017 )

# which data to use
muse <- m17


# summarize taxon counts per sample
m.sum <- muse %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon4 ) %>% 
  summarize( abundance=length(size) )



# make a community dataset
m.meta <- m.sum %>% 
  spread( taxon4, abundance, fill=0 )

meta <- data.frame(m.meta[,c(1:3)])
# meta$year <- factor( meta$year, ordered=T )
comm <- data.frame(m.meta[,-c(1:3)])
names(comm) <- make.cepnames( names(comm) )



# write metadata and community data to disk
write_csv( comm, "output data/Hakai_2017_mesograzer_comm.csv" )
write_csv( meta, "output data/Hakai_2017_mesograzer_meta.csv" )