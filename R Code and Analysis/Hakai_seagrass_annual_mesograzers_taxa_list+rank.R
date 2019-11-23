####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will be used to refine taxa lists 
###     in part by looking at rank abundance patterns over time
###  code by Matt Whalen
###  started on   14 Novemeber 2019
###  
####################################################################



# load libraries
library( tidyverse )
library( vegan )
library( BiodiversityR )


# read data
m <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan
# m$taxon <- gsub( "[.]", " ", m$taxon )
m %>% filter(is.na(taxon4)) %>% select(taxon) %>% unique %>% unlist

# filter taxa and sites
mfilt <- m %>% 
  filter( is.na(remove), !is.na(taxon4), site %in% c("inner chocked","sandspit","triquet north","triquet south") )

# summarize taxon counts per sample
m.sum <- mfilt %>% 
  unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( ID, year, taxon4 ) %>% 
  summarize( abundance=length(size) )

# make a community dataset
m.meta <- m.sum %>% 
  spread( taxon4, abundance, fill=0 )

meta <- data.frame(m.meta[,c(1,2)])
meta$year <- factor( meta$year, ordered=T )
comm <- data.frame(m.meta[,-c(1,2)])
names(comm) <- make.cepnames( names(comm) )



# make rank abundance curves
ra1 <- rankabundance(comm)

par(mar=c(5,4,2,2)+0.1,xpd=T)
rankabunplot(ra1, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=c(1,12000), xlim=c(0,30))

# different curves by year
rankabuncomp(comm, y=meta, factor='year', 
             scale='proportion', legend=FALSE)

# top taxa in each year
ra14 <- rankabundance( comm[meta$year == 2014,] )
ra15 <- rankabundance( comm[meta$year == 2015,] )
ra16 <- rankabundance( comm[meta$year == 2016,] )
ra17 <- rankabundance( comm[meta$year == 2017,] )

windows(12,7)
yrange=c(1,100000)
xrange=c(0,30)
par(mfrow=c(2,2),mar=c(5,4,1,1)+0.1)
rankabunplot(ra14, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=yrange, xlim=xrange, main='2014')
rankabunplot(ra15, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=yrange, xlim=xrange, main='2015')
rankabunplot(ra16, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=yrange, xlim=xrange, main='2016')
rankabunplot(ra17, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=yrange, xlim=xrange, main='2017')
