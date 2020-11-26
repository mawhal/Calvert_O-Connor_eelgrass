### Prokaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 18, 2019 ###
### Date modified: July 07, 2020 ###
### Date last modified: October 19th, 2020 ###

### This code is now updated to remove the contaminants found in the 2016 dataset ###
### Data from 2015, 2017 and 2018 is rarefied to 3,000 reads/sample, and 2016 is not rarefied due to contamination ***this was changed to coverage-based-rarefaction***
### Changed the pipeline in taxa filtering steps to avoid removal of other taxa in that rank (i.e. | is.na(Rank5)) and change in the ordering of filtering
### Added coverage-based rarefaction and saved tables to be used in all analyses
### Organized, cleaned and reviewed annotation in the script to upload to git ###

library(phyloseq)
library(tidyverse)
library(dplyr)
library(data.table)
library(readr)

# Importing files processed with DADA2 #
all_years_16S_unfiltered <- readRDS("Data/prokaryotes/seagrass_16s.full_dataset.unfiltered.phyloseq_format.RDS")

## QUALITY FILTERING TAXA DATA ##

### 1. Remove mitochondrial and chloroplast ASVs
all_years_16S_filtered <- all_years_16S_unfiltered %>%
  subset_taxa(Rank5 != "Mitochondria" | is.na(Rank5)) %>%
  subset_taxa(Rank3 != "Chloroplastida" | is.na(Rank3)) %>% 
  subset_taxa(Rank4 != "Chloroplast" | is.na(Rank4)) %>%
  subset_taxa(Rank5 != "Chloroplast"| is.na(Rank5))  %>% 
  subset_taxa(Rank1 != "Unassigned"| is.na(Rank1))

### 2. Remove contaminants (those were on the 2016 data)
all_years_16S_filtered <- all_years_16S_filtered %>%
  subset_taxa(Rank7 != "Pseudomonas_sp._ANT7125"| is.na(Rank7)) %>% 
  subset_taxa(Rank7 != "Alcaligenes_faecalis"| is.na(Rank7)) %>% 
  subset_taxa(Rank7 != "Pseudomonas_sp._ZJY-246"| is.na(Rank7))
#View(as.data.frame(tax_table(all_years_16S_filtered)))

# THIS STEP (FILTERING NOISE) WAS NOT USED BECAUSE CHAO's PROCEDURE USES LOW FREQUENCY COUNTS TO ESTIMATE SINGLETONS
# FILTERING per sample
# 3. Remove ASVs with less than ~ 2-5 reads in a given sample
# otu <- as.data.frame(otu_table(all_years_16S_filtered))
# otu_table(all_years_16S_filtered)[otu <= 3] <- 0 #free of noise, I set to 3 asvs/sample
# otu2 <- as.data.frame(otu_table(all_years_16S_filtered)) #free of noise

## FILTERING overall
### 4. Remove OTUs with less than N total reads. (N = 250 in example) whole dataset
all_years_16S_filtered <- prune_taxa(taxa_sums(all_years_16S_filtered) >= 250, all_years_16S_filtered) 

### 5. Remove samples with less than N reads. (N = 1000 in example) wholw dataset
all_years_16S_filtered <- prune_samples(sample_sums(all_years_16S_filtered) >= 1000, all_years_16S_filtered)

all_years_16S_filtered

### 6. look at minimum, mean, and maximum sample counts, if desired
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

### include metadata (year column and sample_type_growth), and add it to phyloseq object
year_growth_column <- read.csv("Data/prokaryotes/year_growth_column_16S_ALL_YEARS_FILTERED.csv")
nrow(year_growth_column)

sample_data(all_years_16S_filtered)$year<- year_growth_column$year 

sample_data(all_years_16S_filtered)$growth <- year_growth_column$growth

all_years_16S_filtered_meso <- all_years_16S_filtered %>% subset_samples(survey_type == "meso_quadrat" | survey_type == "meso_quadrats" ) 

all_years_16S_filtered_meso_Zos <- all_years_16S_filtered_meso %>% subset_samples(growth =="old") 

saveRDS(all_years_16S_filtered_meso_Zos, "Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")

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
all_years_16S_filtered_for_CB <- all_years_16S_filtered_meso_Zos
###View(as.data.frame(otu_table(all_years_16S_filtered_for_CB)))
### the following steps are necessary before running the function because it needs the otu table in the phyloseq object to have taxa as rows

### check if taxa are rows in the phyloseq object
taxa_are_rows(all_years_16S_filtered_for_CB)

### transpose so taxa are rows
otu_table(all_years_16S_filtered_for_CB) <- t(otu_table(all_years_16S_filtered_for_CB))
taxa_are_rows(all_years_16S_filtered_for_CB)

### prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(all_years_16S_filtered_for_CB)),
  correct_singletons = T)
View(as.data.frame(otu_table(all_years_16S_filtered_for_CB)))
### check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

saveRDS(all_years_16S_filtered_for_CB, "Data/prokaryotes/all_years_16S_filtered_for_CB.rds")
all_years_16S_filtered_for_CB <- readRDS("Data/prokaryotes/all_years_16S_filtered_for_CB.rds")

### Due to the stochasticity introduced in random subsampling results could be slightly different.
### So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.

### run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)
### 1,000 iterations
all_16S_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=all_years_16S_filtered_for_CB, coverage = 0.8, iter = 1000, replace = F, correct_singletons = TRUE, drop_lowcoverage = F) 
saveRDS(all_16S_COVERAGE_RAREF_1000, "Data/prokaryotes/all_16S_COVERAGE_RAREF_1000.rds")
all_16S_COVERAGE_RAREF_1000 <- readRDS("Data/prokaryotes/all_16S_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_1000 <- all_16S_COVERAGE_RAREF_1000[c(1:1000)]

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
write.csv(average_otu_tables_1000_round, "Data/prokaryotes/prok_average_otu_tables_1000.csv", quote=F, row.names=F )

### add metadata according to #SampleID labels
prok_average_1000 <- read.csv("Data/prokaryotes/prok_average_otu_tables_1000.csv", header=T)
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )
metadata_sel_iter <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))
master_table_iter <- inner_join(metadata_sel_iter , prok_average_1000 , by = "SampleID")
###View(as.data.frame(master_table_iter ))

### exclude the following samples for analyses: "ZosCSPE", "ZosCSPF" # no info if new or old leaf
exclude <- c("ZosCSPE", "ZosCSPF")
master_table_iter <- master_table_iter %>% 
  dplyr::filter(!SampleID %in% exclude)

### recode to site names used by grazers
master_table_iter <- master_table_iter %>% 
  dplyr::mutate(site=recode(site,
                            "choked_south_pigu" = "choked_inner",
                            "choked_flat_island" = "choked_inner",
                            "mcmullin_north" = "mcmullins_north",
                            "mcmullin_south" = "mcmullins_south",
                            "goose_southwest" = "goose_south_west",
                            "goose_southeast" = "goose_south_east",
                            "pruth_bay_south" = "pruth_bay",
                            "pruth_baysouth" = "pruth_bay"))

master_table_iter <- master_table_iter %>% 
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

### for master final table, get only leaf_old 
master_table_iter_final <- master_table_iter %>% 
  dplyr::filter(sample_type =="leaf_old")

### get only meso_quadrat survey 
master_table_iter_final <- master_table_iter_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")

### create a region_year column so can remove only mcmullin 2016 samples
master_table_iter_final <- master_table_iter_final %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

### reorganize column orders (get region_year to first columns together with metadata)
master_table_iter_final <- master_table_iter_final %>%
  dplyr::select(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, everything())

### remove mcmullin 2016 samples
master_table_iter_final <- master_table_iter_final %>% 
  dplyr::filter(!region_year == "mcmullins_2016")

### create a unique site_quadrat_id column
master_table_iter_final <- master_table_iter_final %>% 
  unite(site_quadrat_id, site, quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined
View(master_table_iter_final)

write.csv(master_table_iter_final, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", row.names=F)

######################################################################
### normalized data using coverage based rarefaction FAMILY LEVEL ###
######################################################################

### recreate phyloseq object from ASV level 1,000 iterations ###
prokary_ASV_level_1000_CB <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header = TRUE)

# Separate to metadata and species table
# metadata
prokary.meta <- prokary_ASV_level_1000_CB %>% dplyr::select(c(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, region, site_quadrat_id, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))
names(prokary.meta)
### taxa table
prokary.otus <- prokary_ASV_level_1000_CB %>% dplyr::select(SampleID,starts_with("ASV"))

### Make sure species data frame is numeric matrix and transpose so taxa are rows
prokary.otus_t <- setNames(as.data.frame(t(prokary.otus[,-1])), prokary.otus[,1])
prokary.otus_t_matrix <- data.matrix(prokary.otus_t)

### Desired format for species table
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30

### Make sure metadata has ID as rownames
prokary.meta.ID <- data.frame(prokary.meta[,-1], row.names=prokary.meta[,1])

###  create a taxonomy table if you want (could do in excel?) in the matrix format
# tax.mat <- matrix(taxa.data)

##       Domain Phylum Class Order Family Genus Species
## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o" 
# You don't have to have those specific taxonomic ranks, you could have any number of columns and can be called anything

all_years_16S_filtered_meso_Zos <- readRDS("Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")
prokary.tax <- as.data.frame(unclass(tax_table(all_years_16S_filtered_meso_Zos)))
prokary.tax <- as.matrix(prokary.tax)

#### Save each object in "phyloseq" format to be combined in a phyloseq object ####
OTU <- otu_table(prokary.otus_t_matrix, taxa_are_rows = T)
META <- sample_data(prokary.meta.ID)
TAX <- tax_table(prokary.tax)

### prokary.phyloseq <- phyloseq(OTU, TAX, META)
prokary.phyloseq <- phyloseq(OTU, TAX, META)
sample_names(OTU)
sample_names(META)

### collapse data at the family level
family_level_16S <- prokary.phyloseq  %>%
  tax_glom(taxrank = "Rank5")

### get OTU, taxonomy and metadata tables in a data frame format
df_family_level_16S <- psmelt(family_level_16S)

### get OTU table (currently in long format) to convert to wide format
otu_family_level_16S <- df_family_level_16S %>% 
  select(c(Sample, OTU, Abundance))

otu_family_level_16S <- otu_family_level_16S %>% 
  pivot_wider(names_from = "OTU", values_from = "Abundance")

### get metadata
metadata_family_level_16S <- df_family_level_16S %>% 
  select(c(Sample, year, region_year, region, site_quadrat_id, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id)) %>% 
  unique()

master_family_level_16S <- left_join(metadata_family_level_16S, otu_family_level_16S, by = "Sample")

write.csv(master_family_level_16S, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_family_level_1000_COVERAGE_RAREF.csv", row.names=F)
