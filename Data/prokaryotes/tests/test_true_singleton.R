### Prokaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 18, 2019 ###
### Date modified: July 07, 2020 ###
### Date last modified: September 02, 2020 ###

### This code is now updated to remove the contaminants found in the 2016 dataset ###
### Data from 2015, 2017 and 2018 is rarefied to 3,000 reads/sample, and 2016 is not rarefied due to contamination
### Changed the pipeline in taxa filtering steps to avoid removal of other taxa in that rank (i.e. | is.na(Rank5)) and change in the ordering of filtering
### Added coverage-based rarefaction and saved tables to be used in all analyses

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)
library(dplyr)
library(data.table)

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
#View(as.data.frame(tax_table(all_years_16S_filtered)))
# FILTERING per sample
# 3. Remove ASVs with less than ~ 2-5 reads in a given sample
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

all_years_16S_filtered_meso_Zos.otu <- as.data.frame(otu_table(all_years_16S_filtered_meso_Zos ))

write.csv(all_years_16S_filtered_meso_Zos.otu, file="Data/prokaryotes/16S_ASV_level_otu_table_NOT_RAREFIED.csv", row.names=T)

## add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table <- read.csv(file="Data/prokaryotes/16S_ASV_level_otu_table_NOT_RAREFIED.csv",header=T )
colnames(otu_table)[1]<-"SampleID"

master_table <- right_join(metadata_sel , otu_table , by = "SampleID")
#View(as.data.frame(master_table ))
### Exclude the following samples for analyses: "ZosCSPE", "ZosCSPF" # no info if new or old leaf 
exclude <- c("ZosCSPE", "ZosCSPF")
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
master_table_final <- master_table %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

# reorganize column orders (get region_year to first columns together with metadata)
master_table_final <- master_table_final %>%
  dplyr::select(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, everything())

# remove mcmullin 2016 samples
master_table_final <- master_table_final %>% 
  dplyr::filter(!region_year == "mcmullins_2016")

#create a unique site_quadrat_id column
master_table_final <- master_table_final %>% 
  unite(site_quadrat_id, site, quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined
View(master_table_final)

# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_NOT_RAREFIED.csv", quote=F, row.names=F)

###########################################################
# run the chiu_chao_2016_singleton_count_estimator script #
###########################################################
names(master_table_final)[1:20]
data_16S <- master_table_final %>% 
  select(-c(1:15))

data_16S_num <- sapply(data_16S, as.numeric)
data_16S_integer <- as.integer(data_16S_num)

sing_est_16S <- singleton.Est(data_16S_integer,"abundance")$corrected.data

data_16S_df <- sapply(sing_est_16S, as.data.frame)
