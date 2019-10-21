###
# O'Connor Seagrass Communities neat Hakai Institute's Calvert Island Ecological Observatory
# revisiting questionable specimen identification 
# specimens viewed by Matt Whalen and Felipe Amadeo
# specimens viewed on 25 July 2019


## libraries
library(tidyverse)


## read data
# questionable specimens to revisit are located in an annotated list of unique taxa from the dataset
unique_taxa <- read_csv( "unique_taxa_edited_20190423.csv" )
# define questionable specimens
qs <- unique_taxa %>% filter(revisit_20190725==1)
# the data portion of the dataset
d <- read_csv( "O'Connor_hakai_seagrass_MASTER_grazers.csv" )


## replace all spaces with periods
qs$Original <- gsub( " ", ".", qs$Original )
d$taxon     <- gsub( " ", ".", d$taxon )


## identify the taxa we want in the data
sort(unique(d$taxon))
unique(d$taxon)[unique(d$taxon) %in% qs$Original]
length(which(unique(d$taxon) %in% qs$Original)) # looks like it finds all of them


## grab all of the samples for which we find the specimens
dq <- d[ d$taxon %in% qs$Original, ]


## clean up the dataset so we can look at it more easily
samples <- dq %>%
  group_by( year, date, site, sample, taxon ) %>%
  summarize( abundance = length(size) ) %>%  # grab abundance using the number of times a taxon appears
  spread( key=taxon, value=abundance, fill = 0 )  # spread out the dataset so each row is a unqie sample and columns contain questionable taxa


## write to disk
write_csv( samples, "dubious_specimens_Matt+Felipe_20190725.csv" )


## in which years did the questionable specimen appear?
years <- dq %>% 
  group_by( year, taxon ) %>%
  summarize( abundance=length(size) ) %>%
  spread( key=taxon, value=abundance, fill=0 )

## write to disk
write_csv( years, "dubious_years_Matt+Felipe_20190725.csv" )



## pull out particular samples with Felipe
d %>% filter( year==2017, site == "inner choked", sample==5 ) %>%
  arrange( taxon )

d %>% filter( year==2016, site == "sandspit", sample==2 ) %>%
  arrange( taxon ) %>%
  group_by(taxon) %>%
  summarize( abundance=length(size) )

d %>% filter( year==2016, site == "goose south west", sample==2 ) %>%
  arrange( taxon )

d %>% filter( year==2016, site == "goose south west", sample==2 ) %>%
  arrange( taxon )

d %>% filter( taxon=="Pentidotea.wosnesenskii")

d %>% filter( taxon=="Urothoidae")
d %>% filter( year==2016, site == "triquet north", sample==3 ) %>%
  arrange( taxon ) %>%
  group_by(taxon) %>%
  summarize( abundance=length(size) )
