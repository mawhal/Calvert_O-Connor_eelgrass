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
m <- read_csv( "output_data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output_data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan


# filter taxa and sites
# mfilt <- m %>%
#   filter( is.na(remove), !is.na(taxon4))
# Four core sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon4) )#,
          site %in% c("inner choked","sandspit",
                      "triquet north","triquet south") )


# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon4 ) %>% 
  summarize( abundance=length(size) )


# read in quadrat level data
q <- read_csv( "output_data/O'Connor_hakai_seagrass_MASTER_quadrats.csv" )
q <- q %>% 
  select( year, site, sample=id, shoot.count, biomass ) %>% 
  group_by( year, site, sample ) %>% 
  summarize( shoot.count=mean(shoot.count, na.rm=T), 
             biomass=mean(biomass, na.rm=T) )
sort(unique(q$site))
q$site[ q$site == "triquet bay"] <- "triquet south"
q$site[ q$site == "sand spit"] <- "sandspit"
q$site[ q$site %in% c("inner choked i5", "choked south pigu","flat island")] <- "inner choked"
q$sample <- toupper( q$sample )
q$sample[ q$year==2015 & q$site=="inner choked" ] <- LETTERS[1:6]


# merge mesograzers and quadrats 
mq <- left_join( m.sum, q  )



# make a community dataset
m.meta <- mq %>% 
  spread( taxon4, abundance, fill=0 )

meta <- data.frame(m.meta[,c(1:5)])
# meta$year <- factor( meta$year, ordered=T )
comm <- data.frame(m.meta[,-c(1:5)])
names(comm) <- make.cepnames( names(comm) )

# richness and diversity
meta$richness <- specnumber(comm)
meta$shannon <- diversity(comm)

# add geolocation
geo <- read_csv( "../metadata/coreen_hakai_sites_location_2015.csv" )
geo <- geo %>% 
  select( site, lat, long )

meta <- left_join( meta, geo )
meta$site <- fct_reorder( meta$site, meta$lat )


# total counts per samples
m.sample <- m.sum %>% 
  group_by( year, site, sample ) %>%
  summarize( totalA=sum(abundance) )

m.sample %>% 
  group_by( year,site ) %>% 
  summarize( length(sample) )
  
# total counts per sites
m.totals <- m.sum %>% 
  group_by( year, site ) %>% 
  summarize( totalA=sum(abundance), samples=length(unique(sample))) %>% 
  mutate( count.per.sample=totalA/samples )

# density per eelgrass biomass sample across all taxa
meta <- left_join( meta, m.sample )
meta <- meta %>% 
  mutate( density.per.biomass=totalA/biomass )

# 



pa <- ggplot( m.sample, aes(x=site,y=totalA) ) + facet_wrap(~year, ncol=4) +
  geom_boxplot() +
  geom_point( data=m.totals, aes(x=site,y=count.per.sample),
              size=3, col='slateblue') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10() + ylab( "total epifauna abundance" )
pb <- ggplot( meta, aes(x=site,y=density.per.biomass) ) + facet_wrap(~year, ncol=4) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab( "epifauna density per g eelgrass" )
pc <- ggplot( meta, aes(x=site,y=richness) ) + facet_wrap(~year, ncol=4) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab( "epifauna taxon richness" )
pd <- ggplot( meta, aes(x=site,y=shannon) ) + facet_wrap(~year, ncol=4) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab( "epifauna taxon Shannon diversity" )
pe <- ggplot( meta, aes(x=site,y=biomass) ) + facet_wrap(~year, ncol=4) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab( "eelgrass biomass (grams)" )
pf <- ggplot( meta, aes(x=site,y=shoot.count) ) + facet_wrap(~year, ncol=4) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab( "eelgrass shoot density" )

plot_grid( pe, pf, pa, pb,  ncol = 1 )
plot_grid( pc, pd,  ncol = 1 )

# nmds of comparable data

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


ggplot( meta, aes(MDS1,MDS2,col=site)) + 
  stat_chull( fill=NA ) +
  geom_point( size=3 ) +
  scale_color_manual( values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"))
