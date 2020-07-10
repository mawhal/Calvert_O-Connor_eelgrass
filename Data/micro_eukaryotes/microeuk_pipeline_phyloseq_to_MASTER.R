### Microeukaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: October 03, 2019 ###
### Date last modified: July 07, 2020 ###

### Data from 2015, 2016, 2017 and 2018 is rarefied to 1,000 reads/sample
### Latest update refers to changes in the pipeline in taxa filtering steps to avoid removal of other taxa in that rank (i.e. | is.na(Rank5)) and change in the ordering of filtering

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)


#### Importing files ####
all_years_18S_unfiltered <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/micro_eukaryotes/18S_allyears_unfiltered.RDS")

#### QUALITY FILTERING TAXA DATA ####

# 1. Remove unassigned taxa and zostera OTUs #these are from the host plant, we want only microbial community
all_years_18S_filtered <- all_years_18S_unfiltered %>%
  subset_taxa(Rank1 != "Unassigned"| is.na(Rank1))  %>% 
  subset_taxa(Rank6 != "Zostera"| is.na(Rank6)) 

# 2. Remove seaweeds and plants
all_years_18S_filtered  <- all_years_18S_filtered %>%
  subset_taxa(Rank3 != "Rhodophyceae"| is.na(Rank3))  %>% 
  subset_taxa(Rank5 != "Ulvophyceae"| is.na(Rank5))  %>% 
  subset_taxa(Rank5 != "Ulvales"| is.na(Rank5))  %>% 
  subset_taxa(Rank5 != "Embryophyta"| is.na(Rank5))  %>% 
  subset_taxa(Rank5 != "Phaeophyceae"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Desmarestiales"| is.na(Rank5))  %>%
  subset_taxa(Rank4 != "Charophyta"| is.na(Rank4))  

# 3. Remove large-bodied Metazoans
all_years_18S_filtered  <- all_years_18S_filtered  %>%
  subset_taxa(Rank4 != "Annelida"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Brachiopoda"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Bryozoa"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Cnidaria"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Cyclocoela"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Echinodermata"| is.na(Rank4))   %>%
  subset_taxa(Rank6 != "Mammalia"| is.na(Rank6))   %>%
  subset_taxa(Rank4 != "Mollusca"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Nematoda"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Nemertea"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Platyhelminthes"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Porifera"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Tunicata"| is.na(Rank4))   %>%
  subset_taxa(Rank4 != "Typhlocoela"| is.na(Rank4))  

# 4. Remove large-bodied Metazoans inside Arthropoda
all_years_18S_filtered <- all_years_18S_filtered  %>%
  subset_taxa(Rank5 != "Insecta"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Arachnida"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Ostracoda"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Malacostraca"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Hydroidolina"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Podocopa"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Chromadorea"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Palpata"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Hydroidolina"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Heteroconchia"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Rhabditophora"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Hexacorallia"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Gymnolaemata"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Chromadorea"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Eumalacostraca"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Gastropoda"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Demospongiae"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Brachiopoda"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Scolecida"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Copelata"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Seriata"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Tetrapoda"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Seriata"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Pteriomorphia"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Ascidiacea"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Gastrotricha"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Dorylaimia"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Pteriomorphia"| is.na(Rank5))   %>%
  subset_taxa(Rank5 != "Enopla"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Enoplia"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Enoplea"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Rhabdocoela"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Anopla"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Tentaculata"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Echinodermata"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Polyplacophora"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Tentaculata"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Scyphozoa"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Sipuncula"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Thecostraca"| is.na(Rank5))  %>%
  subset_taxa(Rank6 != "Sessilia"| is.na(Rank6))  %>%
  subset_taxa(Rank5 != "Ectocarpales"| is.na(Rank5))  #removing it inside Maxillopoda

#View(as.data.frame(tax_table(seagrass_set_unfiltered )))

# FILTERING per sample
# 5. Remove ASVs with less than ~ 2-5 reads in a given sample - PER SAMPLE FILTERING
otu <- as.data.frame(otu_table(all_years_18S_filtered)) #get ASV table
otu_table(all_years_18S_filtered)[otu <= 3] <- 0 # for entries where the raw abundance of an ASV in a sample is less than 3 reads, set the raw read count to 0 - free of noise

# FILTERING overall
# 6. Remove OTUs with less than N total reads. (N = 250 in example) whole dataset
all_years_18S_filtered <- prune_taxa(taxa_sums(all_years_18S_filtered) >= 250, all_years_18S_filtered) 

# 7. Remove samples with less than N reads. (N = 1000 in example) whole dataset
all_years_18S_filtered <- prune_samples(sample_sums(all_years_18S_filtered) >= 1000, all_years_18S_filtered)

# 8. look at minimum, mean, and maximum sample counts, if desired
smin <- 
  min(sample_sums(all_years_18S_filtered))
meanreads <- 
  mean(sample_sums(all_years_18S_filtered))
smax <- 
  max(sample_sums(all_years_18S_filtered))
totalreads <- 
  sum(sample_sums(all_years_18S_filtered))


all_years_18S_filtered

#### RAREFY DATA ####
all_years_18S_1000 <- rarefy_even_depth(all_years_18S_filtered,
                                        sample.size = 1000, # Estimated from rarefaction plot
                                        rngseed = 7, # set seed for reproducibility
                                        replace = FALSE)# sample without replacement; slower but more accurate

### Getting number of unique sequences, total reads and mean reads per sample - rarefied data
all_years_18S_1000 # 928 unique sequences (taxa)
reads <- sum(sample_sums(all_years_18S_1000))
reads # 388000 total reads
smean <- mean(sample_sums(all_years_18S_1000))
smean # 1,000 because it is rarefied

##################################################################
############ 18S RAREFIED TO 1,000 READS PER SAMPLE ##############
###################################################################

all_years_18S_1000

all_years_18S_1000.otu <- as.data.frame(otu_table(all_years_18S_1000))

all_years_18S_1000.tax <- as.data.frame(tax_table(all_years_18S_1000))

all_years_18S_1000.sam <- as.data.frame(sample_data(all_years_18S_1000))

write.csv(all_years_18S_1000.otu, file="Data/micro_eukaryotes/18S_ASV_level_otu_table.csv", row.names=T)

write.csv(all_years_18S_1000.tax, file="Data/micro_eukaryotes/18S_ASV_level_taxonomy_table.csv", row.names=T)

write.csv(all_years_18S_1000.sam , file="Data/micro_eukaryotes/18S_ASV_level_metadata_table.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/micro_eukaryotes/18S_ASV_level_metadata_table.csv",header=T )
View(metadata)

otu_table <- read.csv(file="Data/micro_eukaryotes/18S_ASV_level_otu_table.csv",header=T )
colnames(otu_table)[1]<-"SampleID"

master_table <- left_join(metadata , otu_table , by = "SampleID")
View(as.data.frame(master_table ))

### Exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
master_table <- master_table %>% 
  dplyr::filter(!SampleID == "ZosCSPtrans3Amb3")

###recode to site names used by grazers
master_table <- master_table %>%
  dplyr::mutate(site = recode(site,
                              "choked_south_pigu" = "choked_inner",
                              "choked_flat_island" = "choked_inner",
                              "goose_southwest" = "goose_south_west",
                              "mcmullin_north" = "mcmullins_north",
                              "mcmullin_south" = "mcmullins_south",
                              "goose_southeast" = "goose_south_east",
                              "pruth_bay_south" = "pruth_bay"))

master_table <- master_table %>%
  dplyr::mutate(region = recode(region,
                              "mcmullin" = "mcmullins"))

# For mastel final table, get only leaf_old 
master_table_final <- master_table %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_final <- master_table_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")
View(master_table_final)

master_table_final$meso_quadrat_id <- replace(master_table_final$meso_quadrat_id, master_table_final$meso_quadrat_id == "na", NA)

master_table_final$meso_quadrat_id <- replace(master_table_final$meso_quadrat_id, master_table_final$meso_quadrat_id == "0", NA)

#create a unique site_quadrat_id column
master_table_final <- master_table_final %>% 
  unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined

# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", quote=F, row.names=F)


####### COLLAPSE ######

###Collapse data at the genus level
genus_level_18S <- all_years_18S_1000 %>%
  tax_glom(taxrank = "Rank6") 

### SAVING TABLES ###

genus_level_18S

genus_level_18S.otu <- as.data.frame(otu_table(genus_level_18S))

genus_level_18S.tax <- as.data.frame(unclass(tax_table(genus_level_18S)))

genus_level_18S.sam <- as.data.frame(sample_data(genus_level_18S))

write.csv(genus_level_18S.otu, file="Data/micro_eukaryotes/18S_genus_level_otu_table.csv", row.names=T)

write.csv(genus_level_18S.tax, file="Data/micro_eukaryotes/18S_genus_level_taxonomy_table.csv", row.names=T)

write.csv(genus_level_18S.sam , file="Data/micro_eukaryotes/18S_genus_level_metadata_table.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/micro_eukaryotes/18S_genus_level_metadata_table.csv",header=T )

otu_table_genus <- read.csv(file="Data/micro_eukaryotes/18S_genus_level_otu_table.csv",header=T )
colnames(otu_table_genus)[1]<-"SampleID"

master_table_genus <- left_join(metadata , otu_table_genus , by = "SampleID")
View(as.data.frame(master_table_genus ))

### Exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
exclude <- c("ZosCSPtrans3Amb3")
master_table_genus <- master_table_genus %>% 
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
master_table_genus <- master_table_genus %>%
  dplyr::mutate(site = recode(site,
                              "choked_south_pigu" = "choked_inner",
                              "choked_flat_island" = "choked_inner",
                              "mcmullin_north" = "mcmullins_north",
                              "mcmullin_south" = "mcmullins_south",
                              "goose_southwest" = "goose_south_west",
                              "goose_southeast" = "goose_south_east",
                              "pruth_bay_south" = "pruth_bay"))

master_table_genus <- master_table_genus %>%
  dplyr::mutate(region = recode(region,
                                "mcmullin" = "mcmullins"))

# For mastel final table, get only leaf_old 
master_table_genus_final <- master_table_genus %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_genus_final <- master_table_genus_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")
#View(master_table_genus_final)

master_table_genus_final$meso_quadrat_id <- replace(master_table_genus_final$meso_quadrat_id, master_table_genus_final$meso_quadrat_id == "na", NA)

master_table_genus_final$meso_quadrat_id <- replace(master_table_genus_final$meso_quadrat_id, master_table_genus_final$meso_quadrat_id == "0", NA)

#create a unique site_quadrat_id column
master_table_genus_final <- master_table_genus_final %>% 
  unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined

# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_genus_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_genus_level.csv", quote=F, row.names=F) 


####### COLLAPSE ######

###Collapse data at the family level for HMSC
family_level_18S <- all_years_18S_1000 %>%
  tax_glom(taxrank = "Rank5") 

### SAVING TABLES ###

family_level_18S

family_level_18S.otu <- as.data.frame(otu_table(family_level_18S))

family_level_18S.tax <- as.data.frame(unclass(tax_table(family_level_18S)))

family_level_18S.sam <- as.data.frame(sample_data(family_level_18S))

write.csv(family_level_18S.otu, file="Data/micro_eukaryotes/18S_family_level_otu_table.csv", row.names=T)

write.csv(family_level_18S.tax, file="Data/micro_eukaryotes/18S_family_level_taxonomy_table.csv", row.names=T)

write.csv(family_level_18S.sam , file="Data/micro_eukaryotes/18S_family_level_metadata_table.csv", row.names=F)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/micro_eukaryotes/18S_family_level_metadata_table.csv",header=T )

otu_table_family <- read.csv(file="Data/micro_eukaryotes/18S_family_level_otu_table.csv",header=T )
colnames(otu_table_family)[1]<-"SampleID"

master_table_family <- left_join(metadata , otu_table_family , by = "SampleID")
#View(as.data.frame(master_table_family ))

### Exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
exclude <- c("ZosCSPtrans3Amb3")
master_table_family <- master_table_family %>% 
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
master_table_family <- master_table_family %>%
  dplyr::mutate(site = recode(site,
                              "choked_south_pigu" = "choked_inner",
                              "choked_flat_island" = "choked_inner",
                              "mcmullin_north" = "mcmullins_north",
                              "mcmullin_south" = "mcmullins_south",
                              "goose_southwest" = "goose_south_west",
                              "goose_southeast" = "goose_south_east",
                              "pruth_bay_south" = "pruth_bay"))

master_table_family <- master_table_family %>%
  dplyr::mutate(region = recode(region,
                                "mcmullin" = "mcmullins"))

# For mastel final table, get only leaf_old 
master_table_family_final <- master_table_family %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_family_final <- master_table_family_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")

master_table_family_final$meso_quadrat_id <- replace(master_table_family_final$meso_quadrat_id, master_table_family_final$meso_quadrat_id == "na", NA)

master_table_family_final$meso_quadrat_id <- replace(master_table_family_final$meso_quadrat_id, master_table_family_final$meso_quadrat_id == "0", NA)

View(master_table_family_final)
#create a unique site_quadrat_id column
master_table_family_final <- master_table_family_final %>% 
  unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined
#View(as.data.frame(master_table_family ))
# This MASTER table contains samples from choked which we don't have info on quadrat_id on, but we can use those in all analysis that don't require environmental data
write.csv(master_table_family_final, file="Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_family_level.csv", quote=F, row.names=F) 

