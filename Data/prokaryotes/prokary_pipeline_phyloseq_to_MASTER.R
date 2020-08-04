### Prokaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 18, 2019 ###
### Date last modified: July 07, 2020 ###

### This code is now updated to remove the contaminants found in the 2016 dataset ###
### Data from 2015, 2017 and 2018 is rarefied to 3,000 reads/sample, and 2016 is not rarefied due to contamination
### Latest update refers to changes in the pipeline in taxa filtering steps to avoid removal of other taxa in that rank (i.e. | is.na(Rank5)) and change in the ordering of filtering

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)

#### Importing files ####
all_years_16S_unfiltered <- readRDS("Data/prokaryotes/seagrass_16s.full_dataset.unfiltered.phyloseq_format.RDS")

#### QUALITY FILTERING TAXA DATA ####

# 1. Remove mitochondrial and chloroplast ASVs
all_years_16S_filtered <- all_years_16S_unfiltered %>%
  subset_taxa(Rank5 != "Mitochondria" | is.na(Rank5)) %>%
  subset_taxa(Rank3 != "Chloroplastida" | is.na(Rank3)) %>% 
  subset_taxa(Rank4 != "Chloroplast" | is.na(Rank4)) %>%
  subset_taxa(Rank5 != "Chloroplast"| is.na(Rank5))  %>% 
  subset_taxa(Rank1 != "Unassigned"| is.na(Rank1))

# 2. Remove contaminants (those were on the 2016 data)
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank7 != "Pseudomonas_sp._ANT7125"| is.na(Rank7)) %>% 
  subset_taxa(Rank7 != "Alcaligenes_faecalis"| is.na(Rank7)) %>% 
  subset_taxa(Rank7 != "Pseudomonas_sp._ZJY-246"| is.na(Rank7))

# # FILTERING per sample
# # 3. Remove ASVs with less than ~ 2-5 reads in a given sample
# otu <- as.data.frame(otu_table(all_years_16S_filtered))
# otu_table(all_years_16S_filtered)[otu <= 3] <- 0 #free of noise, I set to 3 asvs/sample
# otu2 <- as.data.frame(otu_table(all_years_16S_filtered)) #free of noise

# FILTERING overall
# 4. Remove OTUs with less than N total reads. (N = 250 in example) whole dataset
all_years_16S_filtered <- prune_taxa(taxa_sums(all_years_16S_filtered) >= 250, all_years_16S_filtered) 

# 5. Remove samples with less than N reads. (N = 1000 in example) wholw dataset
all_years_16S_filtered <- prune_samples(sample_sums(all_years_16S_filtered) >= 1000, all_years_16S_filtered)

all_years_16S_filtered

# 6. look at minimum, mean, and maximum sample counts, if desired
smin <- 
  min(sample_sums(all_years_16S_filtered))
meanreads <- 
  mean(sample_sums(all_years_16S_filtered))
smax <- 
  max(sample_sums(all_years_16S_filtered))
totalreads <- 
  sum(sample_sums(all_years_16S_filtered))

get_sample(all_years_16S_filtered)
sample_sums(all_years_16S_filtered)

###  include metadata (year column and sample_type_growth), and add it to phyloseq object
year_growth_column <- read.csv("Data/prokaryotes/year_growth_column_16S_ALL_YEARS_FILTERED.csv")
nrow(year_growth_column)

sample_data(all_years_16S_filtered)$year<- year_growth_column$year 

sample_data(all_years_16S_filtered)$growth <- year_growth_column$growth

all_years_16S_filtered_meso <- all_years_16S_filtered %>% subset_samples(survey_type == "meso_quadrat" | survey_type == "meso_quadrats" ) 

all_years_16S_filtered_meso_Zos <- all_years_16S_filtered_meso %>% subset_samples(growth =="old") 

#View(as.data.frame(otu_table(all_years_16S_filtered_meso_Zos)))

# subset samples from 2015, 2017 and 2018 to rarefy only those to 3,000 * 2016 had lower sequencing depth and will be rarefied to a lower level
all_years_16S_filtered_no_2016 <- all_years_16S_filtered_meso_Zos %>% subset_samples(!year=="2016") 
all_years_16S_filtered_ONLY_2016 <- all_years_16S_filtered_meso_Zos %>% subset_samples(year=="2016") 
# as.data.frame(sample_data(all_years_16S_filtered_ONLY_2016))[["year"]]


#### RAREFY DATA ####

all_years_16S_3000_no_2016 <- rarefy_even_depth(all_years_16S_filtered_no_2016,
                                      sample.size = 3000, # Estimated from rarefaction plot
                                      rngseed = 7, # set seed for reproducibility
                                      replace = FALSE)# sample without replacement; slower but more accurate
# transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU) )

#normalize 2016 
all_years_16S_filtered_ONLY_2016_NORMALIZED <- all_years_16S_filtered_ONLY_2016
##merging 16S years rarefied to 3,000 *** 2016 NOT RAREFIED *BUT* NORMALIZED
all_years_16S_2016_NOT_RAREFIED <- merge_phyloseq(all_years_16S_3000_no_2016, all_years_16S_filtered_ONLY_2016_NORMALIZED)

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

###########################################
### TEST saving non-rarefied data for iNEXT
###########################################

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)

#### Importing files ####
all_years_16S_unfiltered <- readRDS("Data/prokaryotes/seagrass_16s.full_dataset.unfiltered.phyloseq_format.RDS")

#### QUALITY FILTERING TAXA DATA ####

# 1. Remove mitochondrial and chloroplast ASVs
all_years_16S_filtered <- all_years_16S_unfiltered %>%
  subset_taxa(Rank5 != "Mitochondria" | is.na(Rank5)) %>%
  subset_taxa(Rank3 != "Chloroplastida" | is.na(Rank3)) %>% 
  subset_taxa(Rank4 != "Chloroplast" | is.na(Rank4)) %>%
  subset_taxa(Rank5 != "Chloroplast"| is.na(Rank5))  %>% 
  subset_taxa(Rank1 != "Unassigned"| is.na(Rank1))

# 2. Remove contaminants (those were on the 2016 data)
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank7 != "Pseudomonas_sp._ANT7125"| is.na(Rank7)) %>% 
  subset_taxa(Rank7 != "Alcaligenes_faecalis"| is.na(Rank7)) %>% 
  subset_taxa(Rank7 != "Pseudomonas_sp._ZJY-246"| is.na(Rank7))

# FILTERING per sample
# 3. Remove ASVs with less than ~ 2-5 reads in a given sample
otu <- as.data.frame(otu_table(all_years_16S_filtered))
otu_table(all_years_16S_filtered)[otu <= 3] <- 0 #free of noise, I set to 3 asvs/sample
otu2 <- as.data.frame(otu_table(all_years_16S_filtered)) #free of noise

# FILTERING overall
# 4. Remove OTUs with less than N total reads. (N = 250 in example) whole dataset
all_years_16S_filtered <- prune_taxa(taxa_sums(all_years_16S_filtered) >= 250, all_years_16S_filtered) 

# 5. Remove samples with less than N reads. (N = 1000 in example) wholw dataset
all_years_16S_filtered <- prune_samples(sample_sums(all_years_16S_filtered) >= 1000, all_years_16S_filtered)

all_years_16S_filtered


all_years_16S_NONE_RAREFIED.otu <- as.data.frame(otu_table(all_years_16S_filtered ))

all_years_16S_NONE_RAREFIED.tax <- as.data.frame(tax_table(all_years_16S_filtered ))

all_years_16S_NONE_RAREFIED.sam <- as.data.frame(sample_data(all_years_16S_filtered ))

write.csv(all_years_16S_NONE_RAREFIED.otu, file="Data/prokaryotes/16S_ASV_level_otu_table_NONE_RAREFIED.csv", row.names=T)

write.csv(all_years_16S_NONE_RAREFIED.tax, file="Data/prokaryotes/16S_ASV_level_taxonomy_table_NONE_RAREFIED.csv", row.names=T)

write.csv(all_years_16S_NONE_RAREFIED.sam , file="Data/prokaryotes/16S_ASV_level_metadata_table_NONE_RAREFIED.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table <- read.csv(file="Data/prokaryotes/16S_ASV_level_otu_table_NONE_RAREFIED.csv",header=T )
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
write.csv(master_table_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_NONE_RAREFIED_SING_DOUB.csv", quote=F, row.names=F)

##################################################
### TEST rarefying data using coverage based iNEXT
##################################################
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
library(metagMisc)

# phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

# Samples standardized by size will have different degrees of completness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).

# i.e. a seagrass sample with 10,000 total reads will have a different coverage than a seawater sample with 10,000 reads if seawater samples have many more species
all_years_16S_filtered_meso_Zos

taxa_are_rows(all_years_16S_filtered_meso_Zos)

# transpose so taxa are rows
otu_table(all_years_16S_filtered_meso_Zos) <- t(otu_table(all_years_16S_filtered_meso_Zos))

taxa_are_rows(all_years_16S_filtered_meso_Zos)

x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(all_years_16S_filtered_meso_Zos)),
  correct_singletons = T)

SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

# run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)
all_16S_COVERAGE_RAREF <- phyloseq_coverage_raref(physeq=all_years_16S_filtered_meso_Zos, coverage = 0.8, iter = 1, replace = F, correct_singletons = TRUE, drop_lowcoverage = F) 
# I used 1 iteration but see Vladimir's answer below:
#Due to the stochasticity introduced in random subsampling results could be slightly different.
#So you may average diversity estimates or sample dissimilarities across multiple rarefactions.
#I think that it should be more robust in comparison with the results obtained with a single subsampling (ofcourse, it depends on the data).
# I NEED TO RUN WITH MULTIPLE ITERATIONS AND AVERAGE THE RESULTS TO BE MORE ACCURATE
# WHEN RUNNING MULTIPLE ITERATIONS, OUTPUT IS A LIST WITH MULTIPLE PHYLOSEQ OBJECTS

# check sample coverage
attr(all_16S_COVERAGE_RAREF, "SampleCoverage")

all_16S_COVERAGE_RAREF
# to save results from coverage-based rarefaction, need to transpose the otu table again
otu_table(all_16S_COVERAGE_RAREF) <- t(otu_table(all_16S_COVERAGE_RAREF))

all_16S_COVERAGE_RAREF.otu <- as.data.frame(otu_table(all_16S_COVERAGE_RAREF ))

#all_16S_COVERAGE_RAREF.tax <- as.data.frame(tax_table(all_16S_COVERAGE_RAREF ))

all_16S_COVERAGE_RAREF.sam <- as.data.frame(sample_data(all_16S_COVERAGE_RAREF ))

write.csv(all_16S_COVERAGE_RAREF.otu, file="Data/prokaryotes/16S_ASV_level_otu_table_COVERAGE_RAREF.csv", row.names=T)

#write.csv(all_16S_COVERAGE_RAREF.tax, file="Data/prokaryotes/16S_ASV_level_taxonomy_table_COVERAGE_RAREF.csv", row.names=T)

write.csv(all_16S_COVERAGE_RAREF.sam , file="Data/prokaryotes/16S_ASV_level_metadata_table_COVERAGE_RAREF.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table <- read.csv(file="Data/prokaryotes/16S_ASV_level_otu_table_COVERAGE_RAREF.csv",header=T )
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
write.csv(master_table_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_COVERAGE_RAREF.csv", quote=F, row.names=F)

