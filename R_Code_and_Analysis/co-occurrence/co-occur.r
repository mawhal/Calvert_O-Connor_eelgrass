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
comm18 <- comm18[!is.na(comm18$meso_quadrat_id),]
# keep only year, site and meso_quadrat_id metadata from 18S dataset
comm18 <- comm18 %>% 
<<<<<<< HEAD
  select( -SampleID, -host_type, -survey_type, -sample_type, -meso_quadrat_id )
=======
  select( year, site, site_quadrat_id)

>>>>>>> 20f558da6469d6e6806fd6bf7d757a8c8dd4f7b8
comm16 <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_family_level.csv" )
# remove choked samples with no quadrat info from 16S dataset
comm16 <- comm16[!is.na(comm16$quadrat_id),]
# keep only year, site and quadrat_id metadata from 16S dataset
comm16 <- comm16 %>% 
<<<<<<< HEAD
  select( -SampleID, -swab_id, -barcode_plate,-barcode_well,-region_year,
          -host_species, -host_type, -survey_type, -sample_type, -meso_shoot_id, -quadrat_id )
=======
  select( year, site, site_quadrat_id)

>>>>>>> 20f558da6469d6e6806fd6bf7d757a8c8dd4f7b8
commin <- read_csv( "R_Code_and_Analysis/output_data/macroeuk_community_family.csv" )
commin <- commin %>% 
  unite( site_quadrat_id, site, sample, remove = F ) %>% 
  separate( site, into = c("region","sub"), sep="_", remove = F) %>% 
  select( -sub, -sample, -"<NA>" )
commin %>% select( year, site, Ischyroceridae) %>% 
  filter( year!=2014) %>% 
  group_by(year, site) %>% 
  summarize( sum=sum(Ischyroceridae))
  

# merge
dim(comm18);dim(comm16);dim(commin)
# several of the same ASVs in 16S and 18S
commall <- full_join( full_join( comm18, comm16, by = c("year", "region", "site_quadrat_id", "site") ),
                      commin )
commall2 <- full_join( comm18, comm16, by = c("year", "region", "site_quadrat_id", "site") )
commalluse <- commall
commall <- commalluse %>% filter( year != 2014 )
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
range(ccspread$corr)
head( ccspread %>% arrange(-corr), 20 )
head( ccspread %>% arrange(corr), 40 )
ccspread %>% filter( corr<0 & corr>-0.75)
nrow(ccspread %>% filter( corr > 0.9))/nrow(ccspread)
ccspread %>% filter( corr > 0.9) %>% arrange(corr)


ggplot( ccspread, aes( y=corr, x=1)) + 
  geom_violin() +
  geom_hline(yintercept=0, lty=2) #+   geom_point(alpha=0.1, position=position_jitter(width=0.25) )
