####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will clean and merge raw data from different years 
###  code by Matt Whalen
###  started on   22 January 2018
###  updated by Coreen Forbes 22 Feb 2020- updated relative paths 
###  and translated quadrat letters into numbers 2015
####################################################################

## Script goals and notes
# handle each year separately
# determine which taxa are present in each dataset (year), generate taxonomy to deal with tax uncertainty
# determine which sites were sampled in each year, and standardize names across years
# merge data into single long and wide formats and save the output

# libraries
library(tidyverse)
library(reshape2)

#### 2014 DATA--------------------------------------------------------------------------------------------------
### 
# person of record: Nicole Knight

# read data
d14 <- read.csv( "Data/macro_eukaryotes/hakai_grazers_2014.csv", stringsAsFactors = FALSE )

# several extra columns labelled "X","X.1","X.2",etc.
# get rid of these columns
d14 <- d14[,-grep("[X*]",names(d14))]

# several extraneous rows in the dataset
# one shows column sums for each taxon
# others are tacked on NA columns
d14 <- d14[ d14$Sample.number != "", ]

# site "Lower" is actually "Lower Choked"
d14$Site[ d14$Site=="Lower" ] <- "Lower Choked"


# convert to long form for easier manipulation
d14long <- melt( d14, id.vars=1:4, variable.name="taxon", value.name="count" )

# get rid of zeros
l14 <- d14long[ d14long$count!=0, ] 

# add a column for year
l14$year <- 2014

# rename and organize columns
l14 <- l14 %>%
  select( year, date=Date, site=Site, sample=Sample.number, sieve=Sieve.size..mm., taxon, count )
#### END OF 2014 




#### 2015 DATA----------------------------------------------------------------------------------------------------
### 
# person of record: Allison Dennert

# read data
d15 <- read.csv( "Data/macro_eukaryotes/hakai_grazers_2015.csv", stringsAsFactors = FALSE )

##############samples don't match quadrat numbers in 2015, translate between them

##Data with quadrat to sample name

quadrat_sample_match <- read.csv("Data/macro_eukaryotes/2015_quadrat_sample_match_old_names.csv")

sample_match <- quadrat_sample_match %>%
  select(site, quad_number, sample)

sample_match$sample <- toupper(sample_match$sample)

d15 <- left_join(d15, sample_match, by=c("Site"= "site","Sample"= "sample"))

d15 <- d15 %>% 
  mutate(Sample = as.numeric(d15$quad_number))%>%
  select(-c("quad_number"))
  


# lots of NA values in this wide format
# convert these all to zeros
d15[ is.na(d15) ] <- 0

# convert to long form for easier manipulation
d15long <- melt( d15, id.vars=1:4, variable.name="taxon", value.name="count" )

# get rid of zeros
l15 <- d15long[ d15long$count!=0, ] 

# add a column for year
l15$year <- 2015

# rename and organize columns
l15 <- l15 %>%
  select( year, date=Date, site=Site, sample=Sample, sieve=Sieve.size..mm., taxon, count )

# merge with table describing quadrats (samples) as numbers
##This was done above 
#quadid <- read_csv("mapping_18S_quadrat_id_number_letters.csv")
#names(quadid)
#quadid <- quadid %>% 
#  select( site, quad=meso_quadrat_id, sample=meso_shoot_id) %>% 
#  mutate( year=2015 )
#### END OF 2015



#### 2016 DATA----------------------------------------------------------------------------------------------------
### 
# person of record: Tanya Prinzig

# read data
d16 <- read.csv( "Data/macro_eukaryotes/hakai_grazers_2016.csv", stringsAsFactors = FALSE )
# data is already in long format  


# add a column for year
d16$year <- 2016

# rename and organize columns
# Note that sieves were not used after 2015
l16 <- d16 %>%
  select( year, date=date_collected, site, sample=quadrat, 
          taxon=final_id, size=total_length_mm )

str(l16)
# body size (total_length_mm) is not numeric. Coerce to numeric
l16$size <- as.numeric( l16$size )

# store taxonomy in a different place than the count and size data
t16 <- d16 %>%
  select( class, order, suborder, family, genus, taxon=final_id ) %>%
  distinct( )
#### END OF 2016



#### 2017 DATA----------------------------------------------------------------------------------------------------
### 
# person of record: Tanya Prinzig

# read data
d17 <- read.csv( "Data/macro_eukaryotes/hakai_grazers_2017.csv", stringsAsFactors = FALSE )
# data is already in long format

# add a column for year
d17$year <- 2017

# rename and organize columns
# Note that sieves were not used after 2015
l17 <- d17 %>%
  select( year, date=date_collected, site, sample=quadrat, 
          taxon=final_id, size=total_length_mm ) 

# store taxonomy in a different place than the count and size data
t17 <- d17 %>%
  select( class, order, suborder, family, genus, taxon=final_id ) %>%
  distinct( )
#### END OF 2017





# NOTES: Both 2014 and 2015 data started in wide format. These are now in long format
# 2015 currently does not have date information
# 2014 data listed sample numbers as integers, 2015 listed samples as letters (capitalized) <- this is now fixed
# Site names are likely to be different between 2014 and 2015


# site names
sort(unique(l14$site))
sort(unique(l15$site))
sort(unique(l16$site))
sort(unique(l17$site))

with( l14, table(site,sample) )
sort(unique(l14$sample))


with( l15, table(site,sample) )
with( l16, table(site,sample) )
with( l17, table(site,sample) )




# Data starting in 2016 are a little different than previous years
#   because we have length measurements for each individual
# We can either summarize these data, or expand the 2014 and 2015 data to repeat 
#   rows based on the count in each sieve size


#### expand rows for 2014, 2015----------------

## 2014
l14
# this is slow, but works. loop over all rows and replicate them
expand14 <- list()
for( i in 1:nrow(l14) ){
  expand14[[i]] <- l14[i,] %>% slice(rep(1:n(), each = count))  # slice will return a tibble
}
# put list elements back together
expand14 <- do.call( rbind, expand14 )
# we can get rid of the count column now
f14 <- expand14 %>%
  select( year, date, site, sample, size=sieve, taxon )

## 2015
l15
# this is slow, but works. loop over all rows and replicate them
expand15 <- list()
for( i in 1:nrow(l15) ){
  expand15[[i]] <- l15[i,] %>% slice(rep(1:n(), each = count))  # slice will return a tibble
}
# put list elements back together
expand15 <- do.call( rbind, expand15 )
# we can get rid of the count column now
f15 <- expand15 %>%
  select( year, date, site, sample, size=sieve, taxon )

## 2016 + 2017, nothing to do, but go ahead and rename it
f16 <- l16
f17 <- l17

## histograms
windows(4,4)
par( mfrow=c(2,2), mar=c(2,4,2,0)+0.1 )
ylimit <- c(0,1)
hist( f14$size, breaks=8, freq = FALSE, main=2014, ylim=ylimit, col="moccasin" )
hist( f15$size, breaks=8, freq = FALSE, main=2015, ylim=ylimit, col="moccasin" )
# filter out the biggest stuff, which makes the histogram look crazy
hist( l16$size[ l16$size<=8 ], breaks = 8, freq = FALSE, main=2016, ylim=ylimit, col="dodgerblue" ) 
hist( l17$size[ l17$size<=8 ], breaks = 8, freq = FALSE, main=2017, ylim=ylimit, col="dodgerblue" ) 
# measuring size in 2016 on seems to pull things into bigger size classes than those captured by sieves
par( mfrow=c(2,1) )
hist( log(l16$size), breaks = 50, freq = FALSE, main=2016, ylim=ylimit, col="dodgerblue" ) 
hist( log(l17$size), breaks = 50, freq = FALSE, main=2017, ylim=ylimit, col="dodgerblue" ) 


# How many individuals total?
nrow(f14)
nrow(f15)
nrow(f16)
nrow(f17)

# how many quadrats total?
( nq <- c(
length(unique(with( f14, paste(site,sample,sep=".") ))),
length(unique(with( f15, paste(site,sample,sep=".") ))),
length(unique(with( f16, paste(site,sample,sep=".") ))),
length(unique(with( f17, paste(site,sample,sep=".") )))  ) )

# individuals per quadrat
( ipq <- c(
nrow(f14) / length(unique(with( f14, paste(site,sample,sep=".") ))),
nrow(f15) / length(unique(with( f15, paste(site,sample,sep=".") ))),
nrow(f16) / length(unique(with( f16, paste(site,sample,sep=".") ))),
nrow(f17) / length(unique(with( f17, paste(site,sample,sep=".") )))  ) )

data.frame( year=2014:2017, quad.num=nq,ind.per.quad=ipq )

windows(4,4)
par( mar=c(3,4,3,1)+0.1 )
plot( x=factor(2014:2017), y=ipq, las=1, ylab="Number of individuals per quadrat",
      main= "O'Connor Hakai Mesograzers")
text( x=factor(2014:2017), y=ipq, paste( "n =", nq ), pos=c(1,3,1,1) )
# text( x=factor(2014:2017), y=ipq, labels = c("Nicole\nKnight", "Allison\nDennert", 
#                                              "Tanya\nPrinzig","Tanya\nPrinzig"), 
#       pos=c(1,3,1,1))


############ Sort those darn sites out ------------------

# notes about sites from Coreen Forbes and from looking at data sources
# site Goose (2014) is called Goose West (2015)
# site McMullin (2014) is McMullin S (2015) 
# Triquet and Triquet Bay are difficult to pin down
# 2014: Lower and Choked Lower are the same
sort(unique( f14$site ))
sort(unique( f15$site ))
sort(unique( f16$site ))
sort(unique( f17$site ))

f14$site[ f14$site=="Choked" ] <- "choked sandspit"  # not sure about this one
f14$site[ f14$site=="Goose" ] <- "goose south west"
f14$site[ f14$site=="Goose East" ] <- "goose south east"
f14$site[ f14$site=="Lower Choked" ] <- "choked inner"
f14$site[ f14$site=="McMullin" ] <- "mcmullins south"
f14$site[ f14$site=="McMullin North" ] <- "mcmullins north"
f14$site[ f14$site=="Triquet" ] <- "triquet north"
f14$site[ f14$site=="Triquet/No Name Cove" ] <- "triquet south"
f14 <- f14[ !is.na(f14$site), ]

f15$site[ f15$site=="Choked Pass, S. Pigu" ] <- "choked inner"  
f15$site[ f15$site=="Goose E" ] <- "goose south east"
f15$site[ f15$site=="Goose W" ] <- "goose south west"
f15$site[ f15$site=="Goose N" ] <- "goose north"
f15$site[ f15$site=="Sandspit" ] <- "choked sandspit"
f15$site[ f15$site=="McMullins S" ] <- "mcmullins south"
f15$site[ f15$site=="McMullins N" ] <- "mcmullins north"
f15$site[ f15$site=="Triquet N" ] <- "triquet north"
f15$site[ f15$site=="Triquet S" ] <- "triquet south"
# merge with quadrat information- this is already done above in 2015 section
#f15$sample
#right_join( quadid, f15 )

f16$site[ f16$site=="Choked Lower" ] <- "choked inner"  # I am totally not sure about this one either
f16$site[ f16$site=="Goose SW" ] <- "goose south west"
f16$site[ f16$site=="Choked Sandspit" ] <- "choked sandspit"
f16$site[ f16$site=="Pruth Pocket" ] <- "pruth pocket"
f16$site[ f16$site=="Triquet N" ] <- "triquet north"
f16$site[ f16$site=="Triquet S" ] <- "triquet south"

f17$site[ f17$site=="Choked I5" ] <- "choked inner"  # I am totally not sure about this one either
f17$site[ f17$site=="Choked Sandspit" ] <- "choked sandspit"
f17$site[ f17$site=="Pruth Pocket" ] <- "pruth pocket"
f17$site[ f17$site=="Triquet N" ] <- "triquet north"
f17$site[ f17$site=="Triquet S" ] <- "triquet south"



#### merge the data
f15$date <- as.character(f15$date)
f15$sample <- as.character(f15$sample)
f16$sample <- as.character(f16$sample)
f17$sample <- as.character(f17$sample)

m <- full_join(full_join(full_join(f14,f15),f16),f17)

unique(m$site)
m$site <- gsub(" ","_",m$site)

nrow(f14)+nrow(f15)+nrow(f16)+nrow(f17)
dim(m)

####
## write mesograzers to disk
write.csv(m, "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv" )


with(m, table(year,site ))

with(m, table(year,sample ))

m %>%
  group_by(year,site) %>%
  summarise( samples = length(unique(sample))) %>%
  spread( site, samples )



# numbers of individuals per sample across years
mind <- m %>%
  group_by(year,site,sample) %>%
  summarise( individuals = length(size) )

ggplot( mind, aes(x=year,y=individuals)) + facet_wrap(~site) + geom_point() +
  xlim(c(2013.5,2017.5))







##### 
## Taxonomy
# replace all periods with spaces
m$taxon <- gsub( "[.]", " ", m$taxon )
# get rid of double spaces and trailing white space
m$taxon <- gsub( "  ", " ", m$taxon )
m$taxon <- trimws( m$taxon )
unitax <- sort(unique(m$taxon))
write_csv( data.frame(unitax), "Data/R_Code_for_Data_Prep/master_data/unique_taxa_raw.csv" )
m[m$taxon=="",]
#2017 one is nereidae

