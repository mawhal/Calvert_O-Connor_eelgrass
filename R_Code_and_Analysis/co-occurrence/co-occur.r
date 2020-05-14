####################################################################
###  Hakai - O'Connor - Parfrey Seagrass Associates
### 
### Co-occurrence using correlations (across 16S, 18S, and macroinverts)
###  
####################################################################

# load libraries
library( tidyverse )

# read data 
comm18 <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_family_level.csv" )
# remove choked samples with no quadrat info from 18S dataset
comm18 <- comm18 %>% 
  filter(!SampleID == "ZosCSPoldA", !SampleID == "ZosCSPoldF", !SampleID == "ZosCSPoldM")
# keep only year, site and meso_quadrat_id metadata from 18S dataset
comm18 <- comm18 %>% 
  select( year, site, site_quadrat_id)
comm16 <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_family_level.csv" )
# remove choked samples with no quadrat info from 16S dataset
comm16 <- comm16 %>% 
  filter(!SampleID == "ZosCSPoldA", !SampleID == "ZosCSPoldF", !SampleID == "ZosCSPoldG", !SampleID == "ZosCSPoldH", !SampleID == "ZosCSPoldL")
# keep only year, site and quadrat_id metadata from 16S dataset
comm16 <- comm16 %>% 
  select( year, site, site_quadrat_id)
commin <- read_csv( "R_Code_and_Analysis/output_data/macroeuk_community_family.csv" )
commin <- commin %>% 
  unite( site_quadrat_id, site, sample, remove = F ) %>% 
  separate( site, into = c("region","sub"), sep="_", remove = F) %>% 
  select( -sub, -sample, -"<NA>" )

# merge
dim(comm18);dim(comm16);dim(commin)
commall <- full_join( full_join( comm18, comm16 ), commin )
commall <- commall %>% filter( year != 2014 )
commall0 <- commall[ , c(rep(T,6), colSums(commall[,-c(1:6)],na.rm = T) > 2) ]
commallna  <- commall[apply(commall,1,function(z) !any(is.na(z)) ),]
commallna0 <- commallna[ , c(rep(T,6), colSums(commallna[,-c(1:6)],na.rm = T) > 2) ]

comm.cor <- commallna0[ , -c(1:6) ]

cc <- cor(comm.cor)
cc[lower.tri(cc,diag=T)] <- NA
cc <- as.data.frame(cc)
cc$taxon2 <- rownames(cc)
ccspread <- cc %>%  
  gather( "taxon", "corr",-taxon2 ) %>% 
  filter( !is.na(corr) )
# ccspread$taxon2 <- rep( names, 111)
head( ccspread %>% arrange(-corr), 20 )
head( ccspread %>% arrange(corr), 20 )
