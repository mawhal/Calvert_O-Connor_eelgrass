####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will compare epifaunal communities in basic ways
###      by comparing summary patterns of richness, diversity, composition
###  code by Matt Whalen
###  started on   21 November 2019
###  
####################################################################



# load libraries
library( tidyverse )
library( vegan )
library( cowplot )

# read data
m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan

### select sites
# filter taxa and sites
mfilt <- m %>%
  filter( !is.na(taxon2) )


# summarize taxon counts per sample
m.sum.fine <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon=taxon2 ) %>% 
  summarize( abundance=length(size) )
m.spread.fine <- m.sum.fine %>% 
  spread( taxon, abundance, fill=0 )
meta.fine <- data.frame(m.spread.fine[,c(1:3)])
comm.fine <- data.frame(m.spread.fine[,-c(1:3)])
names(comm.fine) <- make.cepnames( names(comm.fine) )
write_csv( m.spread.fine,"R_Code_and_Analysis/output_data/macroeuk_community_finest.csv" )
m.sum.family <- mfilt %>% 
  group_by( year, site, sample, taxon=taxon4 ) %>% 
  summarize( abundance=length(size) )
m.spread.family <- m.sum.family %>% 
  spread( taxon, abundance, fill=0 )
meta.family <- data.frame(m.spread.family[,c(1:3)])
comm.family <- data.frame(m.spread.family[,-c(1:3)])
comm.family <- comm.family[,names(comm.family)!="X.NA."]
names(comm.family) <- make.cepnames( names(comm.family) )
write_csv( m.spread.family,"R_Code_and_Analysis/output_data/macroeuk_community_family.csv" )
m.sum.coarse <- mfilt %>% 
  group_by( year, site, sample, taxon=taxon5 ) %>% 
  summarize( abundance=length(size) )
m.spread.coarse <- m.sum.coarse %>% 
  spread( taxon, abundance, fill=0 )
meta.coarse <- data.frame(m.spread.coarse[,c(1:3)])
comm.coarse <- data.frame(m.spread.coarse[,-c(1:3)])
comm.coarse <- comm.coarse[,names(comm.coarse)!="X.NA."]
names(comm.coarse) <- make.cepnames( names(comm.coarse) )
write_csv( m.spread.coarse,"R_Code_and_Analysis/output_data/macroeuk_community_coarse.csv" )


# - save the bray-curtis distances
for( i in 2014:2017 ){
  metai <- meta.fine %>%
    filter(year==i) %>%
    unite(sample, site, sample, sep="_")
  write_csv( metai, paste0("R_Code_and_Analysis/output_data/",i,"_macroeuk_metadata_finest.csv") )
  sample.names <- make.cepnames(metai$sample)
  commdisti <- vegdist( comm.fine[meta.fine$year==i,], method = "bray" )
  commdisti <- as.matrix(commdisti)
  rownames(commdisti) <- sample.names
  colnames(commdisti) <- sample.names
  write_csv( data.frame(commdisti), paste0("R_Code_and_Analysis/output_data/",i,"_macroeuk_braycurtis_finest.csv") )
}
for( i in 2014:2017 ){
  metai <- meta.family %>%
    filter(year==i) %>%
    unite(sample, site, sample, sep="_")
  write_csv( metai, paste0("R_Code_and_Analysis/output_data/",i,"_macroeuk_metadata_family.csv") )
  sample.names <- make.cepnames(metai$sample)
  commdisti <- vegdist( comm.fine[meta.fine$year==i,], method = "bray" )
  commdisti <- as.matrix(commdisti)
  rownames(commdisti) <- sample.names
  colnames(commdisti) <- sample.names
  write_csv( data.frame(commdisti), paste0("R_Code_and_Analysis/output_data/",i,"_macroeuk_braycurtis_family.csv") )
}
