####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will compare epifaunal communities in basic ways
###      by comparing summary patterns of richness, diversity, composition
###  code by Matt Whalen
###  started on   21 Novemeber 2019
###  
####################################################################



# load libraries
library( tidyverse )
library( vegan )
library( cowplot )

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


# summarize taxon counts per sample
m.sum <- mfilt %>% 
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

# richness and diversity
meta$richness <- specnumber(comm)
meta$shannon <- diversity(comm)

# add geolocation
geo <- read_csv( "../metadata/coreen_hakai_sites_location_2015.csv" )
geo <- geo %>% 
  select( site, lat, long )

meta <- left_join( meta, geo )
# rename and reorder sites
meta$site <- gsub("inner choked","choked_inner",meta$site)
meta$site <- gsub("sandspit","choked_sandspit",meta$site)
meta$site <- gsub(" ","_",meta$site)
meta$site <- fct_reorder( meta$site, -meta$lat )


# nmds of comparable data
# - first save the bray-curtis distances
meta17 <- meta %>% 
  filter(year==2016) %>% 
  unite(sample, site, sample, sep="_")
sample.names <- make.cepnames(meta17$sample)
write_csv( meta17, "output data/2017_grazer_metadata.csv")
commdist17 <- vegdist( comm[meta$year==2017,], method = "bray" )
commdist17 <- as.matrix(commdist17)
rownames(commdist17) <- sample.names
colnames(commdist17) <- sample.names

write_csv( data.frame(commdist17), "output data/2016_grazer_braycurtis.csv")

# NMDS
mds <- metaMDS( comm, distance="bray", k=7 )
mds # very high stress with 2 axes 

# extract points for first two axes
meta <- cbind(meta, mds$points[,c(1,2)])

# convex hull
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

meta$site2 <- factor( meta$site, ordered=T )
man.col <- viridis::viridis(5)
cols <- c( 'black','black',man.col[2],man.col[2],man.col[2],
           man.col[3],man.col[3],man.col[4],man.col[4],man.col[5] )
ggplot( meta, aes(MDS1,MDS2,col=site)) + 
  # stat_chull( fill=NA ) +
  stat_ellipse( aes(lty=site), lwd=1 ) +
  geom_point( size=3 ) +
  scale_linetype_manual( values=c(1,2,1,2,3,1,2,1,2,1)) +
  scale_color_manual( values=cols  )
