### Prokaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 18, 2019 ###
### Date last modified: March 04, 2020 ###

### This code is now updated to remove the contaminants found in the 2016 dataset ###
### Data from 2015, 2017 and 2018 is rarefied to 3,000 reads/sample, and 2016 is not rarefied due to contamination

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)

#### Importing files ####
all_years_16S_unfiltered <- readRDS("Data/prokaryotes/seagrass_16s.full_dataset.unfiltered.phyloseq_format.RDS")

#### QUALITY FILTERING TAXA DATA ####
# 1. look at minimum, mean, and maximum sample counts, if desired
smin <- 
  min(sample_sums(all_years_16S_unfiltered))
smean <- 
  mean(sample_sums(all_years_16S_unfiltered))
smax <- 
  max(sample_sums(all_years_16S_unfiltered))

# 2. Remove samples with less than N reads. (N = 1000 in example) wholw dataset
all_years_16S_filtered <- prune_samples(sample_sums(all_years_16S_unfiltered) >= 1000, all_years_16S_unfiltered)

# 3. Remove OTUs with less than N total reads. (N = 250 in example) whole dataset
all_years_16S_filtered <- prune_taxa(taxa_sums(all_years_16S_filtered) >= 250, all_years_16S_filtered) 

# 4. Remove mitochondrial and chloroplast ASVs
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank5 != "Mitochondria") %>%
  subset_taxa(Rank3 != "Chloroplastida") %>% 
  subset_taxa(Rank4 != "Chloroplast") %>%
  subset_taxa(Rank5 != "Chloroplast") %>% 
  subset_taxa(Rank1 != "Unassigned") 

# 5. Remove contaminants (those were on the 2016 data)
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank7 != "Pseudomonas_sp._ANT7125") %>% 
  subset_taxa(Rank7 != "Alcaligenes_faecalis") %>% 
  subset_taxa(Rank7 != "Pseudomonas_sp._ZJY-246")

# 6. Remove ASVs with less than ~ 2-5 reads in a given sample - PER SAMPLE FILTERING
otu <- as.data.frame(otu_table(all_years_16S_filtered))
otu_table(all_years_16S_filtered)[otu <= 3] <- 0 #free of noise, I set to 3 asvs/sample
otu2 <- as.data.frame(otu_table(all_years_16S_filtered)) #free of noise

all_years_16S_filtered

###  include metadata (year column and sample_type_growth), and add it to phyloseq object
date_year_column <- read.csv("Data/prokaryotes/date_year_column_16S_ALL_YEARS_FILTERED.csv")
nrow(date_year_column)

sample_data(all_years_16S_filtered)$year<- date_year_column$year 

sample_type_growth <- read.csv("Data/prokaryotes/sample_type_growth_column_ALL_YEARS_FILTERED.csv")
nrow(sample_type_growth)

sample_data(all_years_16S_filtered)$sample_type_growth <- sample_type_growth$sample_type_growth 

# subset samples from 2015, 2017 and 2018 to rarefy only those to 3,000 * 2016 had lower sequencing depth and will be rarefied to a lower level
all_years_16S_filtered_no_2016 <- all_years_16S_filtered %>% subset_samples(!year=="2016") 
all_years_16S_filtered_ONLY_2016 <- all_years_16S_filtered %>% subset_samples(year=="2016") 
# as.data.frame(sample_data(all_years_16S_filtered_ONLY_2016))[["year"]]


#### RAREFY DATA ####

all_years_16S_3000_no_2016 <- rarefy_even_depth(all_years_16S_filtered_no_2016,
                                      sample.size = 3000, # Estimated from rarefaction plot
                                      rngseed = 7, # set seed for reproducibility
                                      replace = FALSE)# sample without replacement; slower but more accurate

##merging 16S years rarefied to 3,000 but 2016 NOT RAREFIED
all_years_16S_2016_NOT_RAREFIED <- merge_phyloseq(all_years_16S_3000_no_2016, all_years_16S_filtered_ONLY_2016)

##################################################################
############ 16S RAREFIED TO 3,000 AND 2016  NOT RAREFIED ##########
####################################################################

all_years_16S_2016_NOT_RAREFIED.otu <- as.data.frame(otu_table(all_years_16S_2016_NOT_RAREFIED ))

all_years_16S_2016_NOT_RAREFIED.tax <- as.data.frame(tax_table(all_years_16S_2016_NOT_RAREFIED ))

all_years_16S_2016_NOT_RAREFIED.sam <- as.data.frame(sample_data(all_years_16S_2016_NOT_RAREFIED ))

write.csv(all_years_16S_2016_NOT_RAREFIED.otu, file="Data/prokaryotes/16S_ASV_level_otu_table.csv", row.names=T)

write.csv(all_years_16S_2016_NOT_RAREFIED.tax, file="Data/prokaryotes/16S_ASV_level_taxonomy_table.csv", row.names=T)

write.csv(all_years_16S_2016_NOT_RAREFIED.sam , file="Data/prokaryotes/16S_ASV_level_metadata_table.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table <- read.csv(file="Data/prokaryotes/16S_ASV_level_otu_table.csv",header=T )
colnames(otu_table)[1]<-"SampleID"

master_table <- left_join(metadata_sel , otu_table , by = "SampleID")
View(as.data.frame(master_table ))

### Exclude the following samples for analyses: "ZosCSPE", "ZosCSPF" # no info if new or old leaf and ZosCSPoldM and ZosPBSoldD18 which was all NAs
exclude <- c("ZosCSPE", "ZosCSPF", "ZosPBSoldD18", "ZosCSPoldM")
master_table <- master_table %>% 
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
master_table <- master_table %>% 
  dplyr::mutate(site=recode(site,
                            "choked_south_pigu" = "choked_inner",
                            "choked_flat_island" = "choked_inner",
                            "mcmullin_north" = "mcmullins_north",
                            "mcmullin_south" = "mcmullins_south",
                            "goose_southwest" = "goose_south_west",
                            "goose_southeast" = "goose_south_east",
                            "pruth_bay_south" = "pruth_bay",
                            "pruth_baysouth" = "pruth_bay"))

master_table <- master_table %>% 
  dplyr::mutate(region=recode(region,
                            "mcmullin" = "mcmullins"))

# For mastel final table, get only leaf_old 
master_table_final <- master_table %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_final <- master_table_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")

# create a region_year column so can remove only mcmullin 2016 samples
master_table_final <- master_table_final %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

# reorganize column orders (get region_year to first columns together with metadata)
master_table_final <- master_table_final %>%
  dplyr::select(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, everything())

# remove mcmullin 2016 samples
master_table_final <- master_table_final %>% 
  dplyr::filter(!region_year == "mcmullin_2016")

#create a unique site_quadrat_id column
master_table_final <- master_table_final %>% 
  unite(site_quadrat_id, site, quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined
View(master_table_final)
# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", quote=F, row.names=F)

####### COLLAPSE ######

###Collapse data at the genus level
genus_level_16S <- all_years_16S_2016_NOT_RAREFIED  %>%
  tax_glom(taxrank = "Rank6") 

### SAVING TABLES ###

genus_level_16S

genus_level_16S.otu <- as.data.frame(otu_table(genus_level_16S))

genus_level_16S.tax <- as.data.frame(unclass(tax_table(genus_level_16S)))

genus_level_16S.sam <- as.data.frame(sample_data(genus_level_16S))

write.csv(genus_level_16S.otu, file="Data/prokaryotes/16S_genus_level_otu_table.csv", row.names=T)

write.csv(genus_level_16S.tax, file="Data/prokaryotes/16S_genus_level_taxonomy_table.csv", row.names=T)

write.csv(genus_level_16S.sam , file="Data/prokaryotes/16S_genus_level_metadata_table.csv",row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

#metadata <- metadata %>% dplyr::rename( SampleID = X.SampleID)
metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table_genus <- read.csv(file="Data/prokaryotes/16S_genus_level_otu_table.csv",header=T )
colnames(otu_table_genus)[1]<-"SampleID"

master_table_genus <- left_join(metadata_sel , otu_table_genus , by = "SampleID")

### Exclude the following samples for analyses: "ZosCSPE", "ZosCSPF" # no info if new or old leaf and ZosPBSoldD18 which was all NAs
exclude <-  c("ZosCSPE", "ZosCSPF", "ZosPBSoldD18", "ZosCSPoldM")
master_table_genus <- master_table_genus %>% 
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
master_table_genus <- master_table_genus %>% 
  dplyr::mutate(site=recode(site,
                            "choked_south_pigu" = "choked_inner",
                            "choked_flat_island" = "choked_inner",
                            "mcmullin_north" = "mcmullins_north",
                            "mcmullin_south" = "mcmullins_south",
                            "goose_southwest" = "goose_south_west",
                            "goose_southeast" = "goose_south_east",
                            "pruth_bay_south" = "pruth_bay",
                            "pruth_baysouth" = "pruth_bay"))

master_table_genus <- master_table_genus %>% 
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

# For mastel final table, get only leaf_old 
master_table_final_genus <- master_table_genus %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_final_genus <- master_table_final_genus %>% 
  dplyr::filter(survey_type == "meso_quadrat")

# create a region_year column so can remove only mcmullin 2016 samples
master_table_final_genus <- master_table_final_genus %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

# reorganize column orders (get region_year to first columns together with metadata)
master_table_final_genus <- master_table_final_genus %>%
  dplyr::select(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, everything())

# remove mcmullin 2016 samples
master_table_final_genus <- master_table_final_genus %>% 
  dplyr::filter(!region_year == "mcmullin_2016")

#create a unique site_quadrat_id column
master_table_final_genus <- master_table_final_genus %>% 
  unite(site_quadrat_id, site, quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined

# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_final_genus, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level.csv", quote=F, row.names=F)


###Collapse data at the family level
family_level_16S <- all_years_16S_2016_NOT_RAREFIED  %>%
  tax_glom(taxrank = "Rank5") 

### SAVING TABLES ###
family_level_16S

family_level_16S.otu <- as.data.frame(otu_table(family_level_16S))

family_level_16S.tax <- as.data.frame(unclass(tax_table(family_level_16S)))

family_level_16S.sam <- as.data.frame(sample_data(family_level_16S))

write.csv(family_level_16S.otu, file="Data/prokaryotes/16S_family_level_otu_table.csv", row.names=T)

write.csv(family_level_16S.tax, file="Data/prokaryotes/16S_family_level_taxonomy_table.csv", row.names=T)

write.csv(family_level_16S.sam , file="Data/prokaryotes/16S_family_level_metadata_table.csv",row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

#metadata <- metadata %>% dplyr::rename( SampleID = X.SampleID)
metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table_family <- read.csv(file="Data/prokaryotes/16S_family_level_otu_table.csv",header=T )

colnames(otu_table_family)[1]<-"SampleID"

master_table_family <- left_join(metadata_sel , otu_table_family , by = "SampleID")

### Exclude the following samples for analyses: "ZosCSPE", "ZosCSPF" # no info if new or old leaf and ZosPBSoldD18 which was all NAs
exclude <-  c("ZosCSPE", "ZosCSPF", "ZosPBSoldD18", "ZosCSPoldM")
master_table_family <- master_table_family %>% 
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
master_table_family <- master_table_family %>% 
  dplyr::mutate(site=recode(site,
                            "choked_south_pigu" = "choked_inner",
                            "choked_flat_island" = "choked_inner",
                            "mcmullin_north" = "mcmullins_north",
                            "mcmullin_south" = "mcmullins_south",
                            "goose_southwest" = "goose_south_west",
                            "goose_southeast" = "goose_south_east",
                            "pruth_bay_south" = "pruth_bay",
                            "pruth_baysouth" = "pruth_bay"))

master_table_family <- master_table_family %>% 
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

# For mastel final table, get only leaf_old 
master_table_final_family <- master_table_family %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_final_family <- master_table_final_family %>% 
  dplyr::filter(survey_type == "meso_quadrat")

# create a region_year column so can remove only mcmullin 2016 samples
master_table_final_family <- master_table_final_family %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

# reorganize column orders (get region_year to first columns together with metadata)
master_table_final_family <- master_table_final_family %>%
  dplyr::select(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, everything())

# remove mcmullin 2016 samples
master_table_final_family <- master_table_final_family %>% 
  dplyr::filter(!region_year == "mcmullin_2016")

#create a unique site_quadrat_id column
master_table_final_family <- master_table_final_family %>% 
  unite(site_quadrat_id, site, quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined

# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_final_family, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_family_level.csv", quote=F, row.names=F)

