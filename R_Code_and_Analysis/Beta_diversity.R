####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will investigate beta diversity of grazer (and microbe)
###        community data
###
### started 3 March 2020 by Matt WHalen
###  
####################################################################


library(tidyverse)
library(adespatial)


## Goal calculate local contributions to beta diversity in different years (do we get the same answer?)

# steps - split data into years, make a list of community matrices, perform calculation
#

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

# filter taxa and sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon4))

# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon4 ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean <- m.sum %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, taxon4 ) %>% 
  summarize( abundance=mean(abundance))


# make a community dataset
m.meta <- m.mean %>% 
  spread( taxon4, abundance, fill=0 )

# split by year
comm.years <- split( m.meta, m.meta$year )
lapply(comm.zero, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:3)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:3)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:2])
d <- lapply(comm.zero, function(z) z[,-c(1:2)])

# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBD <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBD <- data.frame( do.call( rbind, LCBD ) )
names(LCBD) <- c("LCBD","p.adj")
LCBD <- data.frame( m.meta[,1:2], LCBD )

windows(3.5,4)
ggplot( data=LCBD, aes( x=year, y=factor(site), size=LCBD, color=p.adj )) + geom_point(alpha=1) +
  # scale_x_discrete( limits=unique(LCBD$method)[c(6,1,3,5)], 
                    # labels=c("visual\nsessile","meta\nsessile","meta\nmobile","visual\nmobile") ) + 
  ylab( "Site" ) + theme_bw( ) 

