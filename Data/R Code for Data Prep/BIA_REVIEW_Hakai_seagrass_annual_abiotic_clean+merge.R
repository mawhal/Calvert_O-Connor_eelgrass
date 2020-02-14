###################################################################
###  Hakai + O'Connor Seagrass Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will clean and merge raw data from different years 
###  code by Matt Whalen
###  started on   22 January 2018
###  last updated 27 September 2018
###  
###################################################################

## Script goals and notes
# read, clean, and merge abiotic data across years

# libraries
library(tidyverse)
library(lubridate)
library(reshape2)

#### Read data -------------------------------------------------------------------------------
### 
# person of record: Nicole Knight
d14 <- NULL
##
##
## 2015
##
##
##
d15 <- read_csv( "~/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/data_oconnor/abiotic/hakai_abiotic_2015.csv" )
a15 <- d15 %>% 
  dplyr::mutate( date=mdy(date) ) %>%
  dplyr::select( date, site, lat=latitude, long=longitude, depth=mean_depth, do=dissolved_oxygen_concentration, salinity, 
          temperature, pressure, par, fluorometry_chlorophyll, turbidity, chl_a, phaeo,
          nitrate_nitrite_ugL, phosphate_ugL, silicate_ugL, bed_area_m2 ) %>%
  distinct()
a15$site <- gsub( x=a15$site, pattern= "_", replacement=" "  )
# rename sites so they are constent across time
a15$site[ a15$site=="choked south pigu" ] <- "choked_inner"
a15$site[ a15$site=="choked sandspit" ] <- "choked_sandspit"
a15$site[ a15$site=="triquet bay south" ] <- "triquet_south"
a15$site[ a15$site=="triquet north" ] <- "triquet_north"
a15$site[ a15$site=="goose southeast" ] <- "goose_south_east"
a15$site[ a15$site=="goose southwest" ] <- "goose_south_west"
a15$site[ a15$site=="mcmullin south" ] <- "mcmullins_south"
a15$site[ a15$site=="mcmullin north" ] <-  "mcmullins_north"
a15$site[ a15$site=="pruth pocket" ] <-  "pruth_pocket"

# NO INFORMATION FOR PRUTH BAY FROM 2015
#
##
##
## 2016
##
##
##
d16 <- read_csv( "~/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/data_oconnor/abiotic/hakai_abiotic_2016.csv" )
geo <- as.data.frame( do.call(rbind, strsplit( d16$coordinates,split="N, ", fixed=TRUE )) )
names(geo) <- c("lat","long")
lat <- gsub( "[(]", "", geo$lat)
latsplit <- strsplit( lat, split="[.]" )
geo$lat <- unlist(lapply( latsplit, function(z) as.numeric(paste(z[2],z[3],sep="."))/60 + as.numeric(z[1]) ))
long <- gsub( "W)", "", geo$long)
longsplit <- strsplit( long, split="[.]" )
geo$long <- unlist(lapply( longsplit, function(z) -( as.numeric(paste(z[2],z[3],sep="."))/60 + as.numeric(z[1]) ) ))
d16.geo <- cbind( d16, geo )
a16 <- d16.geo %>% 
  dplyr::mutate( date=dmy(date) ) %>%
  dplyr::select( date, site, lat, long, depth=`depth (m)`, do=`do (mg/l)`,
          salinity=`salinity (ppt)`, temperature = `temperature (c)`,
          site.depth = `max depth at site (m)`)
# rename sites so they are constent across time
a16$site[ a16$site=="flat island" ] <- "choked_inner"
a16$site[ a16$site=="sandspit" ] <- "choked_sandspit"
a16$site[ a16$site=="triquet bay" ] <- "triquet_south"
a16$site[ a16$site=="triquet north" ] <- "triquet_north"
a16$site[ a16$site=="goose southeast" ] <- "goose_south_east"
a16$site[ a16$site=="goose southwest" ] <- "goose_south_west"
a16$site[ a16$site=="pruth pocket" ] <-  "pruth_pocket"
a16$site[ a16$site=="pruth bay" ] <-  "pruth_bay"

# pick the deepest sample for each site
u16 <- list()
for(i in 1:length(unique(a16$site)) ){
  tmp <- a16[ a16$site == unique(a16$site)[i], ]
  u16[[i]] <- if( all(!is.na(tmp$depth)) ) {
    tmp[ tmp$depth == max(tmp$depth,na.rm=T), ]
  } else {
  tmp
  }}
# combine list
a16 <- do.call( rbind, u16 )
# looks like two YSI over different depths in a bed (maybe at different parts of the tide?) - average these
a16 <- a16 %>% 
  group_by(date,site) %>% 
  summarize( lat=mean(lat,na.rm=T), long=mean(long,na.rm=T),
             depth=mean(depth,na.rm=T), site.depth = mean(site.depth, na.rm=T), 
             do=mean(do,na.rm=T), salinity=mean(salinity,na.rm=T),
             temperature=mean(temperature,na.rm=T)  )
##
##
## 2017
##
##
d17 <- read_csv( "~/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/data_oconnor/abiotic/corrected_hakai_abiotic_2017.csv" )
a17 <- d17 %>%
  dplyr::mutate( date=dmy(date) ) %>%
  dplyr::select( date, site, depth=`depth(m)`, temperature=`temp(c)`, salinity=`salinity(ppt)`, ph )
# rename sites so they are constent across time
a17$site <- gsub( x=a17$site, pattern= "_", replacement=" "  )
a17$site[ a17$site=="choked interior i5" ] <- "choked_inner"
a17$site[ a17$site=="sand spit" ] <- "choked_sandspit"
a17$site[ a17$site=="triquet bay" ] <- "triquet_south"
a17$site[ a17$site=="triquet north" ] <- "triquet_north"
a17$site[ a17$site=="pruth pocket" ] <-  "pruth_pocket"
a17$site[ a17$site=="pruth bay" ] <-  "pruth_bay"

a17$site
# pick the deepest sample for each site
u17 <- list()
for(i in 1:length(unique(a17$site)) ){
  tmp <- a17[ a17$site == unique(a17$site)[i], ]
  u17[[i]] <- if( all(!is.na(tmp$depth)) ) {
    tmp[ tmp$depth == max(tmp$depth,na.rm=T), ]
  } else {
    tmp
  }}
# combine list
a17 <- do.call( rbind, u17 )
##
##
## 2018
##
##
d18 <- read_csv( "~/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/data_oconnor/abiotic/corrected_hakai_abiotic_2018.csv" )
# remove blank rows
d18 <- d18 %>% filter( !is.na(`temp (c)`) )
a18 <- d18 %>%
  dplyr::mutate( date=dmy(date) ) %>%
  dplyr::select( date, site, depth=`sampling_depth (m)`, temperature=`temp (c)`,
          salinity=`salinity (ppt)`,site.depth=`bottom depth (m)`, ph=pH )
# pick the deepest sample for each site
u18 <- list()
for(i in 1:length(unique(a18$site)) ){
  tmp <- a18[ a18$site == unique(a18$site)[i], ]
  u18[[i]] <- if( all(!is.na(tmp$depth)) ) {
    tmp[ tmp$depth == max(tmp$depth,na.rm=T), ]
  } else {
    tmp
  }}
# combine list
a18 <- do.call( rbind, u18 )
# rename sites so they are constent across time
a18$site[ a18$site=="choked I5" ] <- "choked_inner"
a18$site[ a18$site=="sandspit" ] <- "choked_sandspit"
a18$site[ a18$site=="triquet bay" ] <- "triquet_south"
a18$site[ a18$site=="triquet north" ] <- "triquet_north"
a18$site[ a18$site=="pruth pocket" ] <-  "pruth_pocket"
a18$site[ a18$site=="pruth bay" ] <-  "pruth_bay"


##
##
# compare site names and columns names across years
a  <- list( a15, a16, a17, a18 )
lapply( a, function(z) unique(z$site) )
lapply( a, names )
##
##
#### MERGE THE DATA ----------------------------------------------------------------------------
a1516 <- full_join( a15, a16 )
a567 <- full_join( a1516, a17 )
abiotic <- full_join( a567, a18 )
# rename sites
#abiotic$site <- gsub(" ","_",abiotic$site)
#abiotic$site <- gsub("choked_interior","choked_inner",abiotic$site)
#abiotic$site <- gsub("sandspit","choked_sandspit",abiotic$site)
#abiotic$site <- gsub("mcmullin","mcmullins",abiotic$site)
sort(unique(abiotic$site))
abiotic$date
##
##
# make a separate column for year
abiotic <- abiotic %>%
  dplyr::mutate(year = lubridate::year(date))

#### WRITE THE MASTER DATA TO FILE -------------------------------------------------------------
write_csv( abiotic, "~/PostDoc/projects/Calvert_O-Connor_eelgrass/R Code and Analysis/output data/Bia_reviewed_O'Connor_hakai_seagrass_MASTER_abiotic_20200214.csv" )

