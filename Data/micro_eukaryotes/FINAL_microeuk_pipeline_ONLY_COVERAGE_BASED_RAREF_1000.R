### Microeukaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: October 03, 2019 ###
### Date last modified: July 07, 2020 ###
### Date last modified: October 21st, 2020 ###

### Data from 2015, 2016, 2017 and 2018 is rarefied to 1,000 reads/sample
### Latest update refers to changes in the pipeline in taxa filtering steps to avoid removal of other taxa in that rank (i.e. | is.na(Rank5)) and change in the ordering of filtering
### Added coverage-based rarefaction and saved tables to be used in all analyses
### Organized, cleaned and reviewed annotation in the script to upload to git ###

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)
library(data.table)

#### Importing files ####
all_years_18S_unfiltered <- readRDS("Data/micro_eukaryotes/18S_allyears_unfiltered.RDS")

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

# # FILTERING per sample
# # 5. Remove ASVs with less than ~ 2-5 reads in a given sample - PER SAMPLE FILTERING
# otu <- as.data.frame(otu_table(all_years_18S_filtered)) #get ASV table
# otu_table(all_years_18S_filtered)[otu <= 3] <- 0 # for entries where the raw abundance of an ASV in a sample is less than 3 reads, set the raw read count to 0 - free of noise

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

all_years_18S_filtered_meso <- all_years_18S_filtered %>% subset_samples(survey_type == "meso_quadrat") 

all_years_18S_filtered_meso_Zos <- all_years_18S_filtered_meso %>% subset_samples(sample_type =="leaf_old") 

saveRDS(all_years_18S_filtered_meso_Zos, "Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds")

###################################################################
### normalizing data using coverage based rarefaction ASV LEVEL ###
###################################################################
library(devtools)
###install_github("vmikk/metagMisc")
###install_github('JohnsonHsieh/iNEXT')
library(metagMisc)
library(iNEXT)
### phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

### Samples standardized by size will have different degrees of completeness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).
### i.e. a seagrass sample with 10,000 total reads will have a different coverage than a seawater sample with 10,000 reads if seawater samples have many more species

### phyloseq object to be used
all_years_18S_filtered_for_CB <- all_years_18S_filtered_meso_Zos

### the following steps are necessary before running the function because it needs the otu table in the phyloseq object to have taxa as rows

# check if taxa are rows in the phyloseq object
taxa_are_rows(all_years_18S_filtered_for_CB)

# transpose so taxa are rows
otu_table(all_years_18S_filtered_for_CB) <- t(otu_table(all_years_18S_filtered_for_CB))

taxa_are_rows(all_years_18S_filtered_for_CB)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(all_years_18S_filtered_for_CB)),
  correct_singletons = T)

# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

saveRDS(all_years_18S_filtered_for_CB, "Data/micro_eukaryotes/all_years_18S_filtered_for_CB.rds")
all_years_18S_filtered_for_CB <- readRDS("Data/micro_eukaryotes/all_years_18S_filtered_for_CB.rds")

### Due to the stochasticity introduced in random subsampling results could be slightly different.
### So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.

### run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)
### 1,000 iterations
all_18S_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_for_CB, coverage = 0.8, iter = 1000, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_1000, "Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_1000.rds")
all_18S_COVERAGE_RAREF_1000 <- readRDS("Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_1000 <- all_18S_COVERAGE_RAREF_1000[c(1:1000)]

### first, extract otu tables from phyloseq objects
### this is how you do it for a single phyloseq object:
### y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
### now do it for the list of phyloseq objects
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

### average all matrices to get the mean abundance across all iterations
average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
### IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
### add SampleID column back
average_otu_tables_1000_round$SampleID <- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(SampleID, everything())
write.csv(average_otu_tables_1000_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_1000.csv", quote=F, row.names=F )

### add metadata according to #SampleID labels
microeuk_average_1000 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_1000.csv", header=T)
metadata <- read.csv(file="Data/micro_eukaryotes/18S_ASV_level_metadata_table.csv",header=T )
master_table_iter <- inner_join(metadata , microeuk_average_1000 , by = "SampleID")
###View(as.data.frame(master_table_iter ))

### exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
exclude <- c("ZosCSPtrans3Amb3")
master_table_iter <- master_table_iter %>%
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
master_table_iter <- master_table_iter %>%
  dplyr::mutate(site = recode(site,
                              "choked_south_pigu" = "choked_inner",
                              "choked_flat_island" = "choked_inner",
                              "mcmullin_north" = "mcmullins_north",
                              "mcmullin_south" = "mcmullins_south",
                              "goose_southwest" = "goose_south_west",
                              "goose_southeast" = "goose_south_east",
                              "pruth_bay_south" = "pruth_bay"))

master_table_iter_final <- master_table_iter %>%
  dplyr::mutate(region = recode(region,
                                "mcmullin" = "mcmullins"))

master_table_iter_final$meso_quadrat_id <- replace(master_table_iter_final$meso_quadrat_id, master_table_iter_final$meso_quadrat_id == "na", NA)

master_table_iter_final$meso_quadrat_id <- replace(master_table_iter_final$meso_quadrat_id, master_table_iter_final$meso_quadrat_id == "0", NA)

#create a unique site_quadrat_id column
master_table_iter_final <- master_table_iter_final %>%
  unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns

write.csv(master_table_iter_final, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", row.names=F)

######################################################################
### normalized data using coverage based rarefaction FAMILY LEVEL ###
######################################################################
### recreate phyloseq object from ASV level 1,000 iterations ###
microeuk_ASV_level_1000_CB <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header = TRUE)

# Separate to metadata and species table
# metadata
microeuk.meta <- microeuk_ASV_level_1000_CB %>% dplyr::select(c(SampleID, year, region, site_quadrat_id, site, host_type, sample_type, survey_type, meso_quadrat_id))
names(microeuk.meta)
# taxa table
microeuk.otus <- microeuk_ASV_level_1000_CB %>% dplyr::select(SampleID,starts_with("ASV"))

# Make sure species data frame is numeric matrix and transpose so taxa are rows
microeuk.otus_t <- setNames(as.data.frame(t(microeuk.otus[,-1])), microeuk.otus[,1])
microeuk.otus_t_matrix <- data.matrix(microeuk.otus_t)

# Desired format for species table
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30

# Make sure metadata has ID as rownames
microeuk.meta.ID <- data.frame(microeuk.meta[,-1], row.names=microeuk.meta[,1])

#  create a taxonomy table if you want (could do in excel?) in the matrix format
# tax.mat <- matrix(taxa.data)

##       Domain Phylum Class Order Family Genus Species
## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o" 
# You don't have to have those specific taxonomic ranks, you could have any number of columns and can be called anything

all_years_18S_filtered_meso_Zos <- readRDS("Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds")
microeuk.tax <- as.data.frame(unclass(tax_table(all_years_18S_filtered_meso_Zos)))
microeuk.tax <- as.matrix(microeuk.tax)

#### Save each object in "phyloseq" format to be combined in a phyloseq object ####
OTU <- otu_table(microeuk.otus_t_matrix, taxa_are_rows = T)
META <- sample_data(microeuk.meta.ID)
TAX <- tax_table(microeuk.tax)

### microeuk.phyloseq <- phyloseq(OTU, TAX, META)
microeuk.phyloseq <- phyloseq(OTU, TAX, META)
sample_names(OTU)
sample_names(META)

### collapse data at the family level
family_level_18S <- microeuk.phyloseq  %>%
  tax_glom(taxrank = "Rank5")

### get OTU, taxonomy and metadata tables in a data frame format
df_family_level_18S <- psmelt(family_level_18S)

### get OTU table (currently in long format) to convert to wide format
otu_family_level_18S <- df_family_level_18S %>% 
  select(c(Sample, OTU, Abundance))

otu_family_level_18S <- otu_family_level_18S %>% 
  pivot_wider(names_from = "OTU", values_from = "Abundance")

### get metadata
metadata_family_level_18S <- df_family_level_18S %>% 
  select(c(Sample,  year, region, site_quadrat_id, site, host_type, sample_type, survey_type, meso_quadrat_id)) %>% 
  unique()

master_family_level_18S <- left_join(metadata_family_level_18S, otu_family_level_18S, by = "Sample")

write.csv(master_family_level_18S, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_family_level_1000_COVERAGE_RAREF.csv", row.names=F)