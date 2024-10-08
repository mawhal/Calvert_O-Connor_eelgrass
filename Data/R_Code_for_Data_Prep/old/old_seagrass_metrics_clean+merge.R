###################################################################
###  Hakai + O'Connor Seagrass Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will clean and merge raw data from different years 
###  code by Matt Whalen
###  started on   22 January 2018
###  
###################################################################

## Script goals and notes
# read, clean, and merge biometric data across years

# libraries
library(tidyverse)
library(lubridate)
library(reshape2) 

#### Read data-------------------------------------------------------------------------------
### 
# person of record: Nicole Knight
d14 <- NULL
##
##
##


## 2015
##
##
##
# metadata linking od quadrat letters to quadrat numbers
quadnames15 <- read_csv("../../metadata/2015_quadrat_sample_match.csv")
quadnames15 <- quadnames15 %>% select( site, quad_number, id=sample )
# 
# quadrats
q15 <- read_csv( "../data_oconnor/seagrass+epiphytes/Hakai_2015_eelgrass_quadrat.csv" )
q15 <- q15 %>%
  select( site, id=`sample ID`, shoot.count=`# of shoots in quadrat`, 
          biomass=`total dry weight`, unknown.weight.shoots=`weight.#shoots`,
          microepiphyte.bulk=microepiphyte, macroepiphyte.quad=macroepiphyte, 
          drift=`drift seaweed`, anchored=`rooted seaweed` ) %>% 
  mutate( site= gsub(" ","_", site))
q15$site[ q15$site=="mcmullins_nroth"] <- "mcmullins_north"
q15$site[ q15$site=="sandspit"] <- "choked_sandspit"
q15$year <- 2015
q15$id  <- recode(q15$id, woof = "a", hiss = "b", chirp = "c", quack = "d", moo = "e", meow = "f")

# single shoots
s15 <- read_csv( "../data_oconnor/seagrass+epiphytes/Hakai_2015_eelgrass_single_shoots.csv" )
s15 <- s15 %>%
  mutate( date=dmy(date) ) %>%
  select( date, site, id=`sample ID`, shoot.length=`shoot length`, shoot.width=`shoot width`, 
          blades, shoot.biomass=`total dry weight`,microepiphyte.filter=microepiphyte, 
          macroepiphyte.shoot=macroepiphyte, smithora )  # microepiphytes were measured after filtering through GFC filters?
s15$year <- 2015
s15$id  <- recode(s15$id, woof = "a", hiss = "b", chirp = "c", quack = "d", moo = "e", meow = "f")


# merge 2015 data
d15 <- left_join( s15, q15 )

##
##
##
## 2016
##
##
##
# quadrat biomass
b16 <- read_csv( "../data_oconnor/seagrass+epiphytes/Hakai_2016_eelgrass_quadrat_biomass.csv" )
b16 <- b16 %>%
  mutate( date=dmy(date),
          biomass = `foil+ dry weight` - `eelgrass foil weight`,
          detritus = `detritus foil+ dry weight` - `detritus foil weight`,
          epibiont = `epiphyte bag wt.` - `epiphyte bag+ dry wt.` ) %>%
  select( date, site, id=`sample ID`, shoot.count=`# of shoots in quadrat`,  
          biomass, detritus, epibiont.taxon=`epiphyte species`, epibiont )
# make all site names lowercase
b16$site <- tolower(b16$site)
b16$site <- gsub(" ","_",b16$site)
b16$site[ b16$site=="sandspit" ] <- "choked_sandspit"
b16$id <- as.character(b16$id)
b16$year <- 2016
# quadrat shoot measurements
q16 <- read_csv( "../data_oconnor/seagrass+epiphytes/Hakai_2016_eelgrass_quadrat_shoots.csv" )
q16 <- q16 %>%
  mutate( date=dmy(date) ) %>%
  select( date, site, id=`sample ID`, sheath.length=`mean sheath length`, 
          shoot.lengthX1=`shoot length 1`, shoot.lengthX2=`shoot length 2`,
          shoot.widthX1=`shoot width (cm) 1`, shoot.widthX2=`shoot width (cm) 2`,
          bladesX1=`# blades shoot 1`, bladesX2=`# blades shoot 2`  )
q16.melt <- melt( q16, id.vars = 1:4 )
shoot.split <- as.data.frame( do.call( rbind, strsplit( as.character(q16.melt$variable), split = "X" ) ) )
  names( shoot.split ) <- c( "metric","shoot.rep" )
q16.split <- cbind( q16.melt, shoot.split )
q16 <- q16.split %>% 
  select(-variable) %>% 
  group_by(date, site, sheath.length, shoot.rep) %>% 
  spread( metric, value, fill = NA)
# make all site names lowercase
q16$site <- tolower(q16$site)
q16$site <- gsub(" ","_",q16$site)
q16$site[q16$site=="sandspit"] <- "choked_sandspit"
q16$id <- as.character( q16$id )
q16$year <- 2016
# single shoot measurements
s16 <- read_csv( "../data_oconnor/seagrass+epiphytes/Hakai_2016_eelgrass_single.shoots.csv" )
names(s16)[1] <- "date"
s16 <- s16 %>%
  mutate( date=dmy(date),
          microepiphyte.filter = `GF/C with dry sample` - `GF/C filter weight`,  # microepiphyte has negative values
          shoot.biomass = `foil+dry weight` - `foil weight`, 
          bryozoan = `bryozoan dry weight` - `bryozoan foil wt.`, 
          smithora = `smithora dry weight` - `smithora foil weight` ,
          bulk.diatom = `bulk diatoms dry weight` - `bulk diatoms foil weight` ) %>%
  select( date, site, id=`sample ID`, blades, shoot.length=`shoot length`, 
          shoot.width=`shoot width`, shoot.biomass,
          microepiphyte.filter, bryozoan, smithora )  # bulk diatoms removed because only appears once and calculated dry mass is negative
# make all site names lowercase
s16$site <- tolower(s16$site)
s16$site <- gsub(" ","_",s16$site)
s16$site[s16$site=="sandspit"] <- "choked_sandspit"
# think of shoot lengths and widths here as a third replicate from each quadrat
s16$shoot.rep <- factor(3)
# FOIL WEIGHT == EELGRASS
s16$id <- as.character( s16$id )
s16$year <- 2016
# drift seaweed
d16 <- read_csv(  "../data_oconnor/seagrass+epiphytes/Hakai_2016_eelgrass_quadrat_driftseaweed.csv" )
d16 <- d16 %>%
  mutate( date=dmy(date), dry=as.numeric(`drift algae dry wt. +bag wt.`),
          bryozoan.drift=as.numeric(factor(notes))-1,
          drift.mass = dry-`drift algae foil wt.` ) %>%
  select( date, site, id=`sample ID`, taxon=`drift algae species`, 
          drift.mass, bryozoan.drift )
d16$drift.mass[ is.na(d16$drift.mass) ] <- 0.02   
# make all site names lowercase
d16$site <- tolower(d16$site)
d16$site <- gsub(" ","_",d16$site)
d16$site[d16$site=="sandspit"] <- "choked_sandspit"
d16$id <- as.character(d16$id)
d16$year <- 2016
# spread out the taxon columns
d16 <- d16 %>%
  group_by(date,site,id) %>%
  spread( taxon, drift.mass )
# merge 2016 data
# bd16 <- full_join( b16, d16 )
# qs16 <- full_join(q16,s16)
# bds16 <- left_join( bd16, s16  )
bq16 <- left_join(b16,q16)
bqd16 <- left_join(bq16,d16)
##
##
##
## 2017
##
##
##
# quadrat level data
q17 <- read_csv( "../data_oconnor/seagrass+epiphytes/corrected_seagrass_metrics_abiotic_2017_2018/corrected_Hakai_2017_Eelgrass_Biometrics.csv" ) 
q17 <- q17 %>%
  mutate( date=dmy(date), shoot.count=`regular shoot #`+`flowering shoot #`,
          biomass = `live dry weight (g)` - `live paper bag weight (g)`,
          detritus = `detritus dry weight (g)` - `detritus tin weight (g)`,
          smithora = `epi algae dry weight (g)` - `epiphyte bag weight (g)` ) %>%
  select( date, site, id=`quadrat #`, shoot.count, veg.count=`regular shoot #`, 
          flowering.count=`flowering shoot #`, biomass, detritus, smithora,
          shoot.lengthX1=`shoot 1 length`, shoot.lengthX2=`shoot 2 length`, shoot.lengthX3=`shoot 3 length`, 
          shoot.lengthX4=`shoot 4 length`, shoot.lengthX5=`shoot 5 length`,
          sheath.lengthX1=`shoot 1 sheath length`, sheath.lengthX2=`shoot 2 sheath length`, sheath.lengthX3=`shoot 3 sheath length`,
          sheath.lengthX4=`shoot 4 sheath length`, sheath.lengthX5=`shoot 5 sheath length`,
          shoot.widthX1=`shoot 1 width (cm)`, shoot.widthX2=`shoot 2 width (cm)`, shoot.widthX3=`shoot 3 width (cm)`,
          shoot.widthX4=`shoot 4 width (cm)`, shoot.widthX5=`shoot 5 width (cm)`,
          bladesX1=`shoot 1 # blades`, bladesX2=`shoot 2 # blades`, bladesX3=`shoot 3 # blades`,
          bladesX4=`shoot 4 # blades`, bladesX5=`shoot 5 # blades` )
q17.melt <- melt( q17, id.vars = 1:9 )
shoot.split <- as.data.frame( do.call( rbind, strsplit( as.character(q17.melt$variable), split = "X" ) ) )
names( shoot.split ) <- c( "metric","shoot.rep" )
q17.split <- cbind( q17.melt, shoot.split )
q17 <- dcast( q17.split, date+site+id+shoot.count+veg.count+flowering.count+biomass+detritus+smithora+shoot.rep~metric )
q17$id <- as.character( q17$id )
q17$year <- 2017
# shoot level epiphyte data
s17 <- read_csv( "../data_oconnor/seagrass+epiphytes/corrected_seagrass_metrics_abiotic_2017_2018/corrected_Hakai_2017_Epiphyte_Data.csv" ) 
s17 <- s17 %>%
  mutate( date=dmy(date), 
          shoot.biomass=`epiphyte dry weight (g)` - `epi foil weight (g)`,
          microepiphyte=`GF/C with dry sample (g)` - `GF/C weight (g)`) %>%
  select( date, site, id=`quadrat #`, shoot.length=`epi blade length (cm)`, 
          shoot.width=`epi blade width (cm)`,sheath.length=`sheath length (cm)`,
          shoot.biomass, microepiphyte  )
s17$id <- as.character( s17$id )
s17$year <- 2017
# 
# # perhaps deal with 'drift' algae separately
# c(taxonX1=macroalgae.species.1, taxonX2=macroalgae.species.2, taxonX3=macroalgae.species.3, taxonX4=macroalgae.species.4,
#   macroalgae.1.wet.weight..g.,
# macroalgae.1.foil..g.,
# macroalgae.1.dry..g.,
# 
# macroalgae.2.wet.weight..g.,
# macroalgae.2.foil..g.,
# macroalgae.2.dry..g.,
# 
# macroalgae.3.wet.weight..g.,
# macroalgae.3.foil..g.,
# macroalgae.3.dry..g.,
# 
# macroalgae.4.wet.weight..g.,
# macroalgae.4.foil..g.,
# macroalgae.4.dry..g. )
##
##
##
## 2018
##
##
##
# quadrat-level biomass, morphology, drift algae
q18 <- read_csv( "../data_oconnor/seagrass+epiphytes/Hakai_2018_Eelgrass_Biometrics.csv" ) 

q18 <- q18 %>%
  mutate( date=dmy(date), shoot.count=`regular shoot #`+`flowering shoot #`,
          biomass = `live dry weight (g)` - `live paper bag weight (g)`,
          detritus = `detritus dry weight (g)` - `detritus tin weight (g)` ) %>%
  select( date, site, id=`quadrat #`, shoot.count, veg.count=`regular shoot #`, flowering.count=`flowering shoot #`, 
          biomass, detritus, 
          shoot.lengthX1=`shoot 1 length`, shoot.lengthX2=`shoot 2 length`, shoot.lengthX3=`shoot 3 length`, 
          shoot.lengthX4=`shoot 4 length`, shoot.lengthX5=`shoot 5 length`,
          sheath.lengthX1=`shoot 1 sheath length`, sheath.lengthX2=`shoot 2 sheath length`, sheath.lengthX3=`shoot 3 sheath length`,
          sheath.lengthX4=`shoot 4 sheath length`, sheath.lengthX5=`shoot 5 sheath length`,
          shoot.widthX1=`shoot 1 width (cm)`, shoot.widthX2=`shoot 2 width (cm)`, shoot.widthX3=`shoot 3 width (cm)`,
          shoot.widthX4=`shoot 4 width (cm)`, shoot.widthX5=`shoot 5 width (cm)`,
          bladesX1=`shoot 1 # blades`, bladesX2=`shoot 2 # blades`, bladesX3=`shoot 3 # blades`,
          bladesX4=`shoot 4 # blades`, bladesX5=`shoot 5 # blades` )
q18.melt <- melt( q18, id.vars = 1:8 )
shoot.split <- as.data.frame( do.call( rbind, strsplit( as.character(q18.melt$variable), split = "X" ) ) )
names( shoot.split ) <- c( "metric","shoot.rep" )
q18.split <- cbind( q18.melt, shoot.split )
q18 <- dcast( q18.split, date+site+id+shoot.count+veg.count+flowering.count+biomass+detritus+shoot.rep~metric )
q18$id <- as.character( q18$id )
q18$year <- 2018

# shoot level epiphyte data
s18 <- read_csv( "../data_oconnor/seagrass+epiphytes/corrected_seagrass_metrics_abiotic_2017_2018/corrected_Hakai_2018_Epiphyte_Data.csv" ) 
s18 <- s18 %>%
  mutate( date=dmy(date),
          shoot.biomass=`epiphyte dry weight (g)` - `epi foil weight (g)`,
          microepiphyte=`GF/C with dry sample (g)` - `GF/C weight (g)` ) %>%
  select( date, site, id=`quadrat #`, shoot.length=`epi blade length (cm)`,
          shoot.width=`epi blade width (cm)`, sheath.length=`sheath length (cm)`, 
          microepiphyte ) 
s18$id <- as.character( s18$id )
 s18$year <- 2018
# epiphyte.dry.weight refers to the dry mass of the eelgrass blade from which epiphytes were scrapped and quantified
# GF/C with dry sample is the mass of filtered microepiphytes



#####
## MERGE THE DIFFERENT DATA SOURCES
# shoot level
s <- full_join(full_join(full_join(s15,s16),s17),s18)
# quadrat level
q <- full_join(full_join(full_join(q15,bqd16),q17),q18)
#
# d <- full_join(s,q, by=c("site","id","year"))


# rename sites and quadrat ID (sample)
s$site[ s$site == "choked_south_pigu"] <- "choked_inner"
s$site[ s$site == "flat_island"]       <- "choked_inner"
s$site[ s$site == "inner choked i5"]   <- "choked_inner"
s$site[ s$site == "sand spit"]         <- "choked_sandspit"
s$site[ s$site == "triquet_bay"]       <- "triquet_south"
s$site[ s$site == "triquet bay"]       <- "triquet_south"
s$site[ s$site == "triquet north"]     <- "triquet_north"
s$site[ s$site == "goose_east"]        <- "goose_south_east"
s$site[ s$site == "goose_west"]        <- "goose_south_west"
s$site[ s$site == "pruth bay"]         <- "pruth_bay"
s$site[ s$site == "pruth pocket"]      <- "pruth_pocket"
unique(s$site)
with(s, table(site,year))
quadnames15$year <- 2015
s <- left_join( s, quadnames15 )
s$quad_number[ is.na(s$quad_number) ] <- s$id[ is.na(s$quad_number) ]
s <- s %>% 
  mutate( id = as.numeric( s$quad_number ) ) %>% 
  select( -quad_number )

q$site <- gsub(" ","_",q$site)
q$site[ q$site == "goose_east"]        <- "goose_south_east"
q$site[ q$site == "goose_west"]        <- "goose_south_west"
q$site[ q$site == "inner_choked_i5"]   <- "choked_inner"
q$site[ q$site == "flat_island"]       <- "choked_inner"
q$site[ q$site == "sand_spit"]         <- "choked_sandspit"
q$site[ q$site == "mcmullins_nroth"]   <- "mcmullins_north"
q$site[ q$site == "triquet_bay"]       <- "triquet_south"
unique(q$site)
with(q, table(site,year))
unique(q$id)

#####
## Write shoot level and quadrat level data to disk
write_csv( s, "master data/O'Connor_hakai_seagrass_MASTER_shoots.csv" )
write_csv( q, "master data/O'Connor_hakai_seagrass_MASTER_quadrats.csv" )
