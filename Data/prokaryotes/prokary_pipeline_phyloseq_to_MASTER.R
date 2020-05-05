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

# 2. Remove samples with less than N reads. (N = 1000 in example) 
all_years_16S_filtered <- prune_samples(sample_sums(all_years_16S_unfiltered) >= 1000, all_years_16S_unfiltered)

# 3. Remove OTUs with less than N total reads. (N = 250 in example) 
all_years_16S_filtered <- prune_taxa(taxa_sums(all_years_16S_filtered) >= 250, all_years_16S_filtered) 

# 4. Remove mitochondrial and chloroplast OTUs
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank5 != "Mitochondria") %>%
  subset_taxa(Rank3 != "Chloroplastida") %>% 
  subset_taxa(Rank4 != "Chloroplast") %>%
  subset_taxa(Rank5 != "Chloroplast") %>% 
  subset_taxa(Rank1 != "Unassigned") 

# 5. Remove contaminants (those were on the 2016 data)
### WHEN I REMOVE CONTAMINANTS, MCMULLIN 2015 DISSAPEAR
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank7 != "Pseudomonas_sp._ANT7125") %>% 
  subset_taxa(Rank7 != "Alcaligenes_faecalis") %>% 
  subset_taxa(Rank7 != "Pseudomonas_sp._ZJY-246")

# Remove counts that represent less than specified percentage of the total for each sample
otu <- as.data.frame(otu_table(all_years_16S_filtered))
otu_table(all_years_16S_filtered)[otu <= 3] <- 0 #free of noise, I set to 3 otus/sample
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

write.csv(all_years_16S_2016_NOT_RAREFIED.otu, file="Data/prokaryotes/16S_3000_final.otu_REMOVED_CONT_2016_NOT_RARE.csv", row.names=T)

write.csv(all_years_16S_2016_NOT_RAREFIED.tax, file="Data/prokaryotes/16S_3000_final.tax_REMOVED_CONT_2016_NOT_RARE.csv", row.names=T)

write.csv(all_years_16S_2016_NOT_RAREFIED.sam , file="Data/prokaryotes/16S_3000_final.sam_REMOVED_CONT_2016_NOT_RARE.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )

metadata_sel <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))

otu_table <- read.csv(file="Data/prokaryotes/16S_3000_final.otu_REMOVED_CONT_2016_NOT_RARE.csv",header=T )
colnames(otu_table)[1]<-"SampleID"

master_table <- left_join(metadata_sel , otu_table , by = "SampleID")
View(as.data.frame(master_table ))

### Exclude the following samples for analyses that DON’T require metadata:
### choked_exclude –> c(ZosCSPE, ZosCSPF, ZosCSPnewE, ZosCSPnewG, ZosCSPnewH, ZosCSPnewL, ZosCSPnewM, ZosCSPoldL, ZosCSPoldM, ZosCSPnewD, ZosCSPnewC, , ZosCSPnewB, ZosCSPnewB2, waterCSPa, waterCSPb, waterCSPc, ZosCSPnewA)
# For now, I'll just exclude the OLD ones and ones that don't have old or new - we don't know what those are (others will be filtered out anyway)
choked_exclude <- c("ZosCSPE", "ZosCSPF", "ZosCSPoldL", "ZosCSPoldM")
master_table <- master_table %>% 
  dplyr::filter(!SampleID %in% choked_exclude)

###recode to site names used by grazers
master_table <- master_table %>% 
  dplyr::mutate(site=recode(site,
                            "choked_south_pigu" = "choked_inner",
                            "choked_flat_island" = "choked_inner",
                            "goose_southwest" = "goose_south_west",
                            "goose_southeast" = "goose_south_east",
                            "pruth_bay_south" = "pruth_bay",
                            "pruth_baysouth" = "pruth_bay"))

levels(master_table$site)


write.csv(master_table, file="Data/prokaryotes/16S_3000_MASTER_REMOVED_CONT_2016_NOT_RARE.csv", quote=F, row.names=F)

# For mastel final table, get only leaf_old of meso_quadrat survey and remove mcmullin 2016 samples
master_table_final <- master_table %>% 
  dplyr::filter(sample_type =="leaf_old")

master_table_final <- master_table_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")

master_table_final <- master_table_final %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

master_table_final <- master_table_final %>% 
  dplyr::filter(!region_year == "mcmullin_2016")

write.csv(master_table_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", quote=F, row.names=F)
