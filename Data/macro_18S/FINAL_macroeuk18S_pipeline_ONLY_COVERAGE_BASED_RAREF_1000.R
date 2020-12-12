### Microeukaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 02, 2020 ###

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)
library(data.table)
library(vegan)

#### Importing files ####
all_years_18S_unfiltered <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/18S_allyears_unfiltered.RDS")

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
  subset_taxa(Rank4 != "Charophyta"| is.na(Rank4)) %>% 
  subset_taxa(Rank5 != "Ectocarpales"| is.na(Rank5))

# 3. Remove large-bodied Metazoans
all_years_18S_filtered  <- all_years_18S_filtered  %>%
  subset_taxa(Rank4 == "Metazoa_(Animalia)")  # | is.na(Rank4) removed because it was getting other stuff together
all_years_18S_filtered  <- all_years_18S_filtered  %>%
subset_taxa(Rank6 != "Mammalia"| is.na(Rank6))

#View(as.data.frame(tax_table(all_years_18S_filtered )))

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

all_years_18S_filtered_meso_Zos <- all_years_18S_filtered_meso_Zos %>% subset_samples(SampleID !="ZosCSPtrans3Amb3") 

### replace year 2018 wrongly assigned 2017 for two samples
### get OTU, taxonomy and metadata tables in a data frame format
df_all_years_18S <- psmelt(all_years_18S_filtered_meso_Zos)

### get metadata
metadata_18S <- df_all_years_18S %>% 
  select(c(Sample,  year, region, site, host_type, sample_type, survey_type, meso_quadrat_id)) %>% 
  unique() %>% 
  arrange(year)
metadata_18S <- metadata_18S %>% mutate(row_number = row_number())
metadata_18S.ID <- data.frame(metadata_18S[,-1], row.names=metadata_18S[,1])

#relabel 2018 year that got wrongly assigned to 2017
metadata_18S.ID[75, "year"] <- "2018"

sample_data(all_years_18S_filtered_meso_Zos) <- metadata_18S.ID

# get rid of taxa that sums to zero
minTotRelAbun = 0
x = taxa_sums(all_years_18S_filtered_meso_Zos)
keepTaxa = (x / sum(x)) > minTotRelAbun
macro18S_final = prune_taxa(keepTaxa, all_years_18S_filtered_meso_Zos)
#View(as.data.frame(otu_table(macro18S_final)))
saveRDS(macro18S_final, "Data/micro_eukaryotes/macro18S_filtered.rds")
View(as.data.frame(tax_table(macro18S_final)))
##############################
#### INVESTIGATE OUTLIERS ####
##############################
choked_samples <-  macro18S_final %>% 
  subset_samples(region =="choked")
sample_sums(choked_samples)
smin_choked <-
  min(sample_sums(choked_samples))
meanreads_choked <-
  mean(sample_sums(choked_samples))
# overall
sample_sums(macro18S_final)

### get OTU, taxonomy and metadata tables in a data frame format
df_macro18S <- psmelt(macro18S_final)

### get OTU table (currently in long format) to convert to wide format
otu_macro18S <- df_macro18S %>% 
  select(c(Sample, OTU, Abundance))

otu_macro18S <- otu_macro18S %>% 
  pivot_wider(names_from = "OTU", values_from = "Abundance")

### get metadata
metadata_macro18S <- df_macro18S %>% 
  select(c(Sample,  year, region, site, host_type, sample_type, survey_type, meso_quadrat_id)) %>% 
  unique()

metadata_macro18S_SampleID <- metadata_macro18S %>% 
  rename(SampleID = Sample)
write.csv(metadata_macro18S_SampleID, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/metadata_macro18S.csv", row.names = FALSE )

master_macro18S <- left_join(metadata_macro18S, otu_macro18S, by = "Sample")

write.csv(master_macro18S, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/master_macro18S.csv", row.names = FALSE)

# get taxonomy table
taxonomy_macro18S <- df_macro18S %>% 
  select(c(Rank1:Accession)) %>% 
  arrange(Rank5) %>% 
  unique()

write.csv(taxonomy_macro18S, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/taxonomy_macro18S.csv", row.names = FALSE)

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

#macro18S_final <- macro18S_final %>% subset_samples(year != "2017")

### phyloseq object to be used
macro18S_final_for_CB <- macro18S_final

### the following steps are necessary before running the function because it needs the otu table in the phyloseq object to have taxa as rows

# check if taxa are rows in the phyloseq object
taxa_are_rows(macro18S_final_for_CB)

# transpose so taxa are rows
otu_table(macro18S_final_for_CB) <- t(otu_table(macro18S_final_for_CB))

taxa_are_rows(macro18S_final_for_CB)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(macro18S_final_for_CB)),
  correct_singletons = T)

# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)


macro18S_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq= macro18S_final_for_CB, coverage = 0.8, iter = 1000, replace = F, correct_singletons = F, drop_lowcoverage = T)
saveRDS(macro18S_COVERAGE_RAREF_1000, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/macro_18S_COVERAGE_RAREF_1000.rds")
macro18S_COVERAGE_RAREF_1000 <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/macro_18S_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_1000 <- macro18S_COVERAGE_RAREF_1000[c(1:1000)]
subset_phylo_objects_3 <- macro18S_COVERAGE_RAREF_1000[c(1:3)]
# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_3 <- lapply(subset_phylo_objects_3, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# average all matrices to get the mean abundance across all iterations

average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_1000_round$SampleID <- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(SampleID, everything())

average_otu_tables_3 <- Reduce("+",otu_tables_3)/length(otu_tables_3)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_3_round <- average_otu_tables_3 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_3_round$SampleID <- rownames(average_otu_tables_3) 
average_otu_tables_3_round <- average_otu_tables_3_round %>% 
  select(SampleID, everything())

write.csv(average_otu_tables_1000_round, "Data/macro_18S/macro_18S_average_otu_tables_1000.csv", quote=F, row.names=F )
write.csv(average_otu_tables_3_round, "Data/macro_18S/macro_18S_average_otu_tables_3.csv", quote=F, row.names=F )

### 1000 CORRECTED FOR SINGLETONS
macroeuk18S_average_1000_COR_SING <- read.csv("Data/macro_18S/macro_18S_average_otu_tables_1000_COR_SING.csv", header=T)

### add metadata according to #SampleID labels
metadata <- read.csv(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/metadata_macro18S.csv",header=T )
master_table_iter <- inner_join(metadata , macroeuk18S_average_1000_COR_SING , by = "SampleID")
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
  unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) 
#relabel 2018 year that got worngly assigned to 2017
master_table_iter_final["43", "year"] <- "2018"

write.csv(master_table_iter_final, "Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv", row.names=F)
# 
# ### 2000 NOT CORRECTED FOR SINGLETONS
# macroeuk18S_average_2000_NOT_COR_SING <- read.csv("Data/macro_18S/macro_18S_average_otu_tables_2000_NOT_COR_SING.csv", header=T)
# metadata <- read.csv(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/metadata_macro18S.csv",header=T )
# master_table_iter <- inner_join(metadata , macroeuk18S_average_2000_NOT_COR_SING , by = "SampleID")
# ###View(as.data.frame(master_table_iter ))
# 
# ### exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
# exclude <- c("ZosCSPtrans3Amb3")
# master_table_iter <- master_table_iter %>%
#   dplyr::filter(!SampleID %in% exclude)
# 
# ###recode to site names used by grazers
# master_table_iter <- master_table_iter %>%
#   dplyr::mutate(site = recode(site,
#                               "choked_south_pigu" = "choked_inner",
#                               "choked_flat_island" = "choked_inner",
#                               "mcmullin_north" = "mcmullins_north",
#                               "mcmullin_south" = "mcmullins_south",
#                               "goose_southwest" = "goose_south_west",
#                               "goose_southeast" = "goose_south_east",
#                               "pruth_bay_south" = "pruth_bay"))
# 
# master_table_iter_final <- master_table_iter %>%
#   dplyr::mutate(region = recode(region,
#                                 "mcmullin" = "mcmullins"))
# 
# master_table_iter_final$meso_quadrat_id <- replace(master_table_iter_final$meso_quadrat_id, master_table_iter_final$meso_quadrat_id == "na", NA)
# 
# master_table_iter_final$meso_quadrat_id <- replace(master_table_iter_final$meso_quadrat_id, master_table_iter_final$meso_quadrat_id == "0", NA)
# 
# #create a unique site_quadrat_id column
# master_table_iter_final <- master_table_iter_final %>%
#   unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns
# 
# write.csv(master_table_iter_final, "Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_2000_NOT_COR_SING_COVERAGE_RAREF.csv", row.names=F)
# 
# # 
# # ### 5000 NOT CORRECTED FOR SINGLETONS
# # macroeuk18S_average_5000_NOT_COR_SING <- read.csv("Data/macro_18S/macro_18S_average_otu_tables_5000_NOT_COR_SING.csv", header=T)
# # metadata <- read.csv(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/metadata_macro18S.csv",header=T )
# # master_table_iter <- inner_join(metadata , macroeuk18S_average_5000_NOT_COR_SING , by = "SampleID")
# # ###View(as.data.frame(master_table_iter ))
# # 
# # ### exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
# # exclude <- c("ZosCSPtrans3Amb3")
# # master_table_iter <- master_table_iter %>%
# #   dplyr::filter(!SampleID %in% exclude)
# # 
# # ###recode to site names used by grazers
# # master_table_iter <- master_table_iter %>%
# #   dplyr::mutate(site = recode(site,
# #                               "choked_south_pigu" = "choked_inner",
# #                               "choked_flat_island" = "choked_inner",
# #                               "mcmullin_north" = "mcmullins_north",
# #                               "mcmullin_south" = "mcmullins_south",
# #                               "goose_southwest" = "goose_south_west",
# #                               "goose_southeast" = "goose_south_east",
# #                               "pruth_bay_south" = "pruth_bay"))
# # 
# # master_table_iter_final <- master_table_iter %>%
# #   dplyr::mutate(region = recode(region,
# #                                 "mcmullin" = "mcmullins"))
# # 
# # master_table_iter_final$meso_quadrat_id <- replace(master_table_iter_final$meso_quadrat_id, master_table_iter_final$meso_quadrat_id == "na", NA)
# # 
# # master_table_iter_final$meso_quadrat_id <- replace(master_table_iter_final$meso_quadrat_id, master_table_iter_final$meso_quadrat_id == "0", NA)
# # 
# # #create a unique site_quadrat_id column
# # master_table_iter_final <- master_table_iter_final %>%
# #   unite(site_quadrat_id, site, meso_quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns
# # 
# # write.csv(master_table_iter_final, "Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_5000_NOT_COR_SING_COVERAGE_RAREF.csv", row.names=F)
# # 
# # ######################################################################
# # ### normalized data using coverage based rarefaction FAMILY LEVEL ###
# # ######################################################################
# # ### recreate phyloseq object from ASV level 1,000 iterations ###
# # macroeuk18S_1000_CB <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_level_1000_COVERAGE_RAREF.csv", header = TRUE)
# # 
# # # Separate to metadata and species table
# # # metadata
# # macroeuk18S.meta <- macroeuk18S_1000_CB %>% dplyr::select(c(SampleID, year, region, site_quadrat_id, site, host_type, sample_type, survey_type, meso_quadrat_id))
# # names(macroeuk18S.meta)
# # # taxa table
# # macroeuk18S.otus <- macroeuk18S_1000_CB %>% dplyr::select(SampleID,starts_with("ASV"))
# # 
# # # Make sure species data frame is numeric matrix and transpose so taxa are rows
# # macroeuk18S.otus_t <- setNames(as.data.frame(t(macroeuk18S.otus[,-1])), macroeuk18S.otus[,1])
# # macroeuk18S.otus_t_matrix <- data.matrix(macroeuk18S.otus_t)
# # 
# # # Desired format for species table
# # ##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
# # ## OTU1       96      50      36      35      59      80      83      63
# # ## OTU2       52      67      39      39      37      57      20      15
# # ## OTU3       94      18      15      11      14      75       1      12
# # ## OTU4       27      88      98     100      59      27      30      30
# # 
# # # Make sure metadata has ID as rownames
# # macroeuk18S.meta.ID <- data.frame(macroeuk18S.meta[,-1], row.names=macroeuk18S.meta[,1])
# # 
# # #  create a taxonomy table if you want (could do in excel?) in the matrix format
# # # tax.mat <- matrix(taxa.data)
# # 
# # ##       Domain Phylum Class Order Family Genus Species
# # ## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
# # ## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
# # ## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
# # ## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o" 
# # # You don't have to have those specific taxonomic ranks, you could have any number of columns and can be called anything
# # 
# # macro18S_final <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/macro_18S/macro18S_final.rds")
# # macro18S.tax <- as.data.frame(unclass(tax_table(macro18S_final)))
# macro18S.tax <- as.matrix(macro18S.tax)
# 
# #### Save each object in "phyloseq" format to be combined in a phyloseq object ####
# OTU <- otu_table(macroeuk18S.otus_t_matrix, taxa_are_rows = T)
# META <- sample_data(macroeuk18S.meta.ID)
# TAX <- tax_table(macro18S.tax)
# 
# ### microeuk.phyloseq <- phyloseq(OTU, TAX, META)
# macro18S.phyloseq <- phyloseq(OTU, TAX, META)
# sample_names(OTU)
# sample_names(META)
# 
# ### collapse data at the family level
# family_macro18S <- macro18S.phyloseq  %>%
#   tax_glom(taxrank = "Rank5")
# 
# ### get OTU, taxonomy and metadata tables in a data frame format
# df_family_macro18S <- psmelt(family_macro18S)
# 
# ### get OTU table (currently in long format) to convert to wide format
# otu_family_macro18S <- df_family_macro18S %>% 
#   select(c(Sample, OTU, Abundance))
# 
# otu_family_macro18S <- otu_family_macro18S %>% 
#   pivot_wider(names_from = "OTU", values_from = "Abundance")
# 
# ### get metadata
# metadata_family_macro18S <- df_family_macro18S %>% 
#   select(c(Sample,  year, region, site_quadrat_id, site, host_type, sample_type, survey_type, meso_quadrat_id)) %>% 
#   unique()
# 
# master_family_macro18S <- left_join(metadata_family_macro18S, otu_family_macro18S, by = "Sample")
# 
# write.csv(master_family_macro18S, "Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_family_level_1000_COVERAGE_RAREF.csv", row.names=F)
