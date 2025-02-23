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
m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv" )
# # replace spaces with periods for consistency and merging names
# m$taxon <- gsub( " ", ".", m$taxon )
# # bring in the updated data 
# mtaxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
# # merge
# m <- left_join( m,  mtaxa_update )
# # replace periods with spaces for making simple names in vegan
# 
# ### select sites
# # filter taxa and sites
# mfilt <- m %>%
#   filter( !is.na(taxon2) )
# 
# 
# # summarize taxon counts per sample
# m.sum.fine <- mfilt %>% 
#   # unite( "ID", year,site,sample, remove=FALSE ) %>% 
#   group_by( year, site, sample, taxon=taxon2 ) %>% 
#   summarize( abundance=length(size) )
# m.spread.fine <- m.sum.fine %>% 
#   spread( taxon, abundance, fill=0 )
# meta.fine <- data.frame(m.spread.fine[,c(1:3)])
# comm.fine <- data.frame(m.spread.fine[,-c(1:3)])
# names(comm.fine) <- make.cepnames( names(comm.fine) )
# write_csv( m.spread.fine,"R_Code_and_Analysis/output_data/macroeuk_community_finest.csv" )
# m.sum.family <- mfilt %>% 
#   group_by( year, site, sample, taxon=taxon4 ) %>% 
#   summarize( abundance=length(size) )
# m.spread.family <- m.sum.family %>% 
#   spread( taxon, abundance, fill=0 )
# meta.family <- data.frame(m.spread.family[,c(1:3)])
# comm.family <- data.frame(m.spread.family[,-c(1:3)])
# comm.family <- comm.family[,names(comm.family)!="X.NA."]
# names(comm.family) <- make.cepnames( names(comm.family) )
# write_csv( m.spread.family,"R_Code_and_Analysis/output_data/macroeuk_community_family.csv" )
# m.sum.coarse <- mfilt %>% 
#   group_by( year, site, sample, taxon=taxon5 ) %>% 
#   summarize( abundance=length(size) )
# m.spread.coarse <- m.sum.coarse %>% 
#   spread( taxon, abundance, fill=0 )
# meta.coarse <- data.frame(m.spread.coarse[,c(1:3)])
# comm.coarse <- data.frame(m.spread.coarse[,-c(1:3)])
# comm.coarse <- comm.coarse[,names(comm.coarse)!="X.NA."]
# names(comm.coarse) <- make.cepnames( names(comm.coarse) )
# write_csv( m.spread.coarse,"R_Code_and_Analysis/output_data/macroeuk_community_coarse.csv" )

# metadata
m.meta <- m[1:7]

# community data
m.comm <- m[8:ncol(m)]

# - save the bray-curtis distances
for( i in 2015:2017 ){
  path <- "R_Code_and_Analysis/betadiversity/Bray-Curtis/"
  metai <- m.meta %>%
    filter(year==i) #%>%
    # unite(sample, site, sample, sep="_")
  # filter zero columns and rows
  comm.tmp <- m.comm[m.meta$year==i,]
  keep.col <- which(colSums(comm.tmp)>0)
  keep.row <- which(rowSums(comm.tmp)>0)
  meta.use <- metai[ keep.row, ]
  comm.use <- comm.tmp[ keep.row,keep.col ]
  write_csv( meta.use, paste0(path,i,"_macroeuk_metadata.csv") )
  sample.names <- make.cepnames(meta.use$ID)
  commdisti <- vegdist( comm.use, method = "bray" )
  commdisti <- as.matrix(commdisti)
  rownames(commdisti) <- sample.names
  colnames(commdisti) <- sample.names
  write_csv( data.frame(commdisti), paste0(path,i,"_macroeuk_braycurtis_finest.csv") )
  # commdisti <- vegdist( comm.family[meta.family$year==i,], method = "bray" )
  # commdisti <- as.matrix(commdisti)
  # rownames(commdisti) <- sample.names
  # colnames(commdisti) <- sample.names
  # write_csv( data.frame(commdisti), paste0(path,i,"_macroeuk_braycurtis_family.csv") )
  # commdisti <- vegdist( comm.coarse[meta.coarse$year==i,], method = "bray" )
  # commdisti <- as.matrix(commdisti)
  # rownames(commdisti) <- sample.names
  # colnames(commdisti) <- sample.names
  # write_csv( data.frame(commdisti), paste0(path,i,"_macroeuk_braycurtis_coarse.csv") )
}
