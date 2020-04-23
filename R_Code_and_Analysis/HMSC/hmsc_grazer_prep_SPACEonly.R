####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will prepare data for running HMSC models 
###      on grazer community data
###  code by Matt Whalen
###  started on   18 February 2020
###  
####################################################################



# load libraries
library(Hmsc)
library(tidyverse)
library(viridis)
library(corrplot)
library(vegan)


### read data
# read data
m <- read_csv( "../Data/R Code for Data Prep/master data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output_data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan

### select sites
# filter taxa and sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon4))


### more data to read
# explanatory data
explanatory <- read_csv("output_data/merged_explanatory_20200214.csv")

# prepare for merging 
ab <- explanatory
# filter out missing quads
ab <- ab %>% filter( !is.na(quadrat) )


########   Data preparation and model specification


### function to prepare model data and model specification for any given year
model.prep <- function( community=mfilt, year=2015 ){ 
  
  pick = year
  muse <- community %>% 
    filter( year==pick )
  
  
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
  
  # rename quadrats
  meta <- meta %>% 
    mutate( quadrat_id = as.numeric(as.factor(meta$sample)), quadrat = paste( site,quadrat_id, sep="_" ) )
  
  
  # get rid of some taxa
  # rename to Ygrazer
  Ygrazer <- comm
  Ygrazer <- Ygrazer[ , rev(order(colSums(Ygrazer))) ]
  
  # abundance cutoff?
  Ygrazer <- Ygrazer[, colSums(Ygrazer) > 3]
  
  
  # merge with metadata
  ab2 <- ab %>% 
    select( region, site, lat, long) %>% 
    distinct()
  ab2[ab2$site=="goose_north",c("lat","long")] <- c(52.001869,	-128.425817)
  X <- left_join( meta, ab2, by=c("site") )
  
  
  # fixed and random effects (quadrat, site, region level)
  X <- X %>% 
    mutate( region=unlist(lapply(strsplit(site,"_"),function(z)z[1])),
            long=-long ) %>% 
    select( region, lat, long, site, quadrat )
  X <-  as.data.frame(X)
  X$region <- factor(X$region)
  X$site <- factor(X$site)
  X$quadrat <- factor(X$quadrat)
  
  # spatial data from lats and longs
  # make longitude positive
  spatial <- X  %>% select(site,lat,long) %>% distinct()
  spatial <- as.data.frame(spatial)
  rownames(spatial) <- spatial$site
  spatial <- spatial[,2:3]
  
  studyDesign <- X
  
  
  # formats
  Ygrazer <- as.matrix(Ygrazer)
  Xenvironmental <- as.data.frame(X)
  
  
  
  # random effects structure
  rL1 = HmscRandomLevel( units = unique(studyDesign$site_quad) ) #quadrat level
  # rL2 = HmscRandomLevel( units = unique(studyDesign$site) )
  # rL3 = HmscRandomLevel( units = unique(studyDesign$region) )
  rL4 = HmscRandomLevel( sData = spatial )
  
  
  
  ## formula for fixed effects
  XFormula = ~( region + site  ) # add bed area when available
  
  
  # the model
  mgrazer <- Hmsc( Y = Ygrazer, 
                     XData = Xenvironmental, XFormula = XFormula, 
                     distr = "poisson", studyDesign = studyDesign,
                     ranLevels = list(quadrat=rL1, 
                                      #Site=rL2,
                                      # region=rL3,
                                      site=rL4) )
  
  return(mgrazer)

}




# which years are available?
years <- sort(unique(mfilt$year))
# get rid of 2014
years <- years[years != 2014]
# # for now remove 2015 because quadrats not yet resolved
# years <- years[years != 2015]

mgrazers <- lapply(years, function(z) model.prep(year=z) )
mgrazer17 <- model.prep( year=2016 )


# saving
model <- paste0("grazer",years)
ModelDir <- paste0( here::here(), "/R_Code_and_Analysis/output_data")
filename = file.path(ModelDir, paste("model_",as.character(model),
                                     "_specification.Rdata",sep = ""))

mapply( function(x,y) save(x,file=y), x=mgrazers, y=filename )



