### Inverts data frame for creating phyloseq ###
### Author: Bianca Trevizan Segovia ### 
### Date created: October 29th, 2020 ###

#### Finest level ####

### Read table metadata and abundances
m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan
# define regions as first word of site
m$region <- unlist(lapply( strsplit(m$site,split = "_"), function(z) z[1]))
# unique sample ID to differential samples from different sites
m$ID <- with(m, paste(site,sample,sep = "_"))
# change year to character
m$group <- paste0( "year",m$year )

# filter taxa and sites # taxon 2 = finest
mfilt_finest <- m %>%
  dplyr::select( year, group, region, site, sample, ID, taxon = taxon2, remove, size ) %>% 
  dplyr::filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum_finest <- mfilt_finest %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean_finest <- m.sum_finest %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance)) 

# make a community dataset
m.meta_finest <- m.sum_finest %>% 
  spread( taxon, abundance, fill=0 ) 

names(m.meta_finest)[1:16]

m.meta_finest <- m.meta_finest %>% 
  ungroup( sample, ID,year, group, region, site )

write.csv(m.meta_finest, file="Data/R_Code_for_Data_Prep/master_data/MASTER_grazers_finest_to_phyloseq.csv", quote=F, row.names=F)

#### Family level ####

### Read table metadata and abundances
m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan
# define regions as first word of site
m$region <- unlist(lapply( strsplit(m$site,split = "_"), function(z) z[1]))
# unique sample ID to differential samples from different sites
m$ID <- with(m, paste(site,sample,sep = "_"))
# change year to character
m$group <- paste0( "year",m$year )

# filter taxa and sites # taxon 2 = family
mfilt_family <- m %>%
  dplyr::select( year, group, region, site, sample, ID, taxon = taxon4, remove, size ) %>% 
  dplyr::filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum_family <- mfilt_family %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean_family <- m.sum_family %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance)) 

# make a community dataset
m.meta_family <- m.sum_family %>% 
  spread( taxon, abundance, fill=0 ) 

names(m.meta_family)[1:16]

m.meta_family <- m.meta_family %>% 
  ungroup( sample, ID,year, group, region, site )

write.csv(m.meta_family, file="Data/R_Code_for_Data_Prep/master_data/MASTER_grazers_family_to_phyloseq.csv", quote=F, row.names=F)
