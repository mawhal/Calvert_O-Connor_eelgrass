####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will compare epifaunal communities using package mobr
###  code by Matt Whalen
###  started on   12 December 2019
###  
####################################################################



# load libraries
library( tidyverse )
library( mobr )

# colors
man.col <- viridis::viridis(5)
cols <- c( man.col[1],man.col[2],man.col[3],man.col[4],man.col[5] )

# read data
m <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan


# filter taxa and sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon2))
# # Four core sites
# mfilt <- m %>%
#   filter( is.na(remove), !is.na(taxon2),
#           site %in% c("inner choked","sandspit",
#                       "triquet north","triquet south") )


# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon2 ) %>% 
  summarize( abundance=length(size) )




# read in quadrat level data
q <- read_csv( "output data/O'Connor_hakai_seagrass_MASTER_quadrats.csv" )
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


# start with 2017
m17 <- mq %>% 
  filter( year==2016 )



# make a community dataset
m.meta <- m17 %>% 
  spread( taxon2, abundance, fill=0 )

# identify metadata and data columns
meta <- data.frame(m.meta[,c(1:5)])
# meta$year <- factor( meta$year, ordered=T )
comm <- data.frame(m.meta[,-c(1:5)])
names(comm) <- vegan::make.cepnames( names(comm) )


# add geolocation
geo <- read_csv( "../metadata/00_Hakai_UBC_metadata_MASTER - geolocation_site.csv" )
geo <- geo %>% 
  select( site=site_name, lat, long, depth=`depth (m, chart datum)`, 
          bed_area = area_m2_2016 )
geo$site <- gsub("_"," ", geo$site)

meta <- left_join( meta, geo )
# meta$site <- fct_reorder( meta$site, meta$lat )
meta$sample <- as.numeric( meta$sample )

# calculate lat and long for individual samples (quadrats)
# length of one degree latitude in meters at 51.7 degrees
degree_meter_lat <- 111261.61
# length of one degree longitude in meters at 51.7 degrees
degree_meter_long <- 69136.13


# a 3x3 grid each separated by two meters, with center set to zero
mat_lat <- rep(c(-2,0,2),each=3)
mat_long <- rep(c(-2,0,2),3)

grid_lat  <- (mat_lat * 1/degree_meter_lat)
grid_long <- (mat_long * 1/degree_meter_long)

lat_long_convert <- data.frame( sample=1:9, grid_lat, grid_long )

meta <- left_join(meta, lat_long_convert)

plot_attr <- meta %>% 
  mutate( y = lat + grid_lat, x = long + grid_long ) %>% 
  select( site, shoot.count, biomass, depth, bed_area, x, y)

# rename sites
plot_attr$site <- factor( plot_attr$site, 
                          levels=c("pruth pocket","sandspit","inner choked",
                                   "triquet south", "triquet north","goose south west"),
                          labels=c("pruth","sand","inner","triq_n","triq_s","goose_sw"))
plot_attr <- plot_attr %>% 
  arrange(site)
  

# MOBR analysis
# prep data
inv_mob_in <- make_mob_in( comm, plot_attr )

# simple plots
plot_rarefaction( inv_mob_in, "site", "spat" )
par( mfrow=c(1,2) )
plot_abu( inv_mob_in, "site", 'rad', col=c("black",cols) )
plot_abu( inv_mob_in, "site", 'sad', col=c("black",cols) )

# two-scale analysis
inv_stats <- get_mob_stats(inv_mob_in, group_var = "site",
                           index = c("N", "S", "S_n", "S_PIE","pct_rare"),
                           n_perm = 200)
windows(11,4)
plot(inv_stats, index="N", col=cols )
plot(inv_stats, index="S", col=cols)
plot(inv_stats, index="S_n", col=cols)
plot(inv_stats, index="S_PIE", col=cols)
plot(inv_stats, index="pct_rare", col=cols)

# Multi-scale analysis
inv_deltaS = get_delta_stats(inv_mob_in, 'site', ref_group='inner',
                             type='discrete', log_scale=TRUE, n_perm = 200)

plot( inv_deltaS, "sand","inner", display="ddelta S" )

par(mfrow=c(1,5), mar=c(5,4,2,1)+0.1)
overlap_effects( inv_deltaS, "goose_sw")

env_raw=inv_mob_in$env
inv_deltaS = get_delta_stats(inv_mob_in, group_var="site", env_var='bed_area', 
                             type='continuous', log_scale=TRUE, n_perm = 200)
