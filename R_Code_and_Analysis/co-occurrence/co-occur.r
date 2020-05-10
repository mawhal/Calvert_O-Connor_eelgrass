####################################################################
###  Hakai - O'Connor - Parfrey Seagrass Associates
### 
### Co-occurrence using correlations (across 16S, 18S, and macroinverts)
###  
####################################################################

# load libraries
library( tidyverse )

# read data 
comm18 <- read_csv( "Data/micro_eukaryotes/family_18S_core_microbes.csv" )
comm16 <- read_csv( "Data/prokaryotes/family_16S_core_microbes.csv" )
comm16 <- comm16 %>% 
  select( -quadrat_id, -meso_shoot_id )
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
