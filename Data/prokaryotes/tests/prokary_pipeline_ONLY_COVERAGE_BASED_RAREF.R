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

#### Importing files processed with DADA2 ####
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

# THIS STEP (FILTERING NOISE) WAS NOT USED BECAUSE CHAO's PROCEDURE USES LOW FREQUENCY COUNTS TO ESTIMATE SINGLETONS
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

#View(as.data.frame(otu_table(all_years_16S_filtered_meso_Zos)))

#########################################################
### normalizing data using coverage based rarefaction ###
#########################################################
library(devtools)
#install_github("vmikk/metagMisc")
#install_github('JohnsonHsieh/iNEXT')
library(metagMisc)
library(iNEXT)
# phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

# Samples standardized by size will have different degrees of completeness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).
# i.e. a seagrass sample with 10,000 total reads will have a different coverage than a seawater sample with 10,000 reads if seawater samples have many more species

# phyloseq object to be used
all_years_16S_filtered_meso_Zos
View(as.data.frame(otu_table(all_years_16S_filtered_meso_Zos)))
### the following steps are necessary before running the function because it needs the otu table in the phyloseq object to have taxa as rows

# check if taxa are rows in the phyloseq object
taxa_are_rows(all_years_16S_filtered_meso_Zos)

# transpose so taxa are rows
otu_table(all_years_16S_filtered_meso_Zos) <- t(otu_table(all_years_16S_filtered_meso_Zos))
taxa_are_rows(all_years_16S_filtered_meso_Zos)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(all_years_16S_filtered_meso_Zos)),
  correct_singletons = T)
View(as.data.frame(otu_table(all_years_16S_filtered_meso_Zos)))
# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

#saveRDS(all_years_16S_filtered_meso_Zos, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")
# 
all_years_16S_filtered_meso_Zos <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")

#Due to the stochasticity introduced in random subsampling results could be slightly different.
#So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.
# 
# # run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)
all_16S_COVERAGE_RAREF_5 <- phyloseq_coverage_raref(physeq=all_years_16S_filtered_meso_Zos, coverage = 0.8, iter = 5, replace = F, correct_singletons = TRUE, drop_lowcoverage = F) 
saveRDS(all_16S_COVERAGE_RAREF_5, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_5.rds")
all_16S_COVERAGE_RAREF_5 <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_5.rds")

all_16S_COVERAGE_RAREF_50 <- phyloseq_coverage_raref(physeq=all_years_16S_filtered_meso_Zos, coverage = 0.8, iter = 50, replace = F, correct_singletons = TRUE, drop_lowcoverage = F) 
saveRDS(all_16S_COVERAGE_RAREF_50, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_50.rds")
all_16S_COVERAGE_RAREF_50 <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_50.rds")

all_16S_COVERAGE_RAREF_100 <- phyloseq_coverage_raref(physeq=all_years_16S_filtered_meso_Zos, coverage = 0.8, iter = 100, replace = F, correct_singletons = TRUE, drop_lowcoverage = F) 
saveRDS(all_16S_COVERAGE_RAREF_100, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_100.rds")
all_16S_COVERAGE_RAREF_100 <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_100.rds")

# With 200 iterations *but trying to get 1,000 in Entamoeba*
all_16S_COVERAGE_RAREF_200 <- phyloseq_coverage_raref(physeq=all_years_16S_filtered_meso_Zos, coverage = 0.8, iter = 200, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_16S_COVERAGE_RAREF_200, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_200.rds")
all_16S_COVERAGE_RAREF_200 <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_200.rds")

# read RDS ran with 200 iterations but without removing low frequency counts for comparison
all_16S_COVERAGE_RAREF_200_remove_LFC <- readRDS("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/prokaryotes/all_16S_COVERAGE_RAREF_200_REMOVED_LOW_COUNT.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_5 <- all_16S_COVERAGE_RAREF_5[c(1:5)]
subset_phylo_objects_50 <- all_16S_COVERAGE_RAREF_50[c(1:50)]
subset_phylo_objects_100 <- all_16S_COVERAGE_RAREF_100[c(1:100)]
subset_phylo_objects_200 <- all_16S_COVERAGE_RAREF_200[c(1:200)]
subset_phylo_objects_200_remove_LFC <- all_16S_COVERAGE_RAREF_200_remove_LFC[c(1:200)]
subset_phylo_objects_1000 <- all_16S_COVERAGE_RAREF_1000[c(1:1000)]

# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
otu_tables_5 <- lapply(subset_phylo_objects_5, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_50 <- lapply(subset_phylo_objects_50, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_100 <- lapply(subset_phylo_objects_100, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_200 <- lapply(subset_phylo_objects_200, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_200_remove_LFC <- lapply(subset_phylo_objects_200_remove_LFC, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# # average all matrices to get the mean abundance across all iterations
average_otu_tables_5 <- Reduce("+",otu_tables_5)/length(otu_tables_5)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_5_round <- average_otu_tables_5 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_5_round$SampleID <- rownames(average_otu_tables_5)
average_otu_tables_5_round <- average_otu_tables_5_round %>%
  select(SampleID, everything())

average_otu_tables_50 <- Reduce("+",otu_tables_50)/length(otu_tables_50)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_50_round <- average_otu_tables_50 %>% dplyr::mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_50_round$SampleID <- rownames(average_otu_tables_50)
average_otu_tables_50_round <- average_otu_tables_50_round %>%
  select(SampleID, everything())

average_otu_tables_100 <- Reduce("+",otu_tables_100)/length(otu_tables_100)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_100 <- average_otu_tables_100 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_100_round$SampleID <- rownames(average_otu_tables_100)
average_otu_tables_100_round <- average_otu_tables_100_round %>%
  select(SampleID, everything())

average_otu_tables_200 <- Reduce("+",otu_tables_200)/length(otu_tables_200)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_200_round <- average_otu_tables_200 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_200_round$SampleID <- rownames(average_otu_tables_200)
average_otu_tables_200_round <- average_otu_tables_200_round %>%
  select(SampleID, everything())

average_otu_tables_200_remove_LFC <- Reduce("+",otu_tables_200_remove_LFC)/length(otu_tables_200_remove_LFC)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_200_round_remove_LFC <- average_otu_tables_200_remove_LFC %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_200_round_remove_LFC$SampleID <- rownames(average_otu_tables_200_remove_LFC)
average_otu_tables_200_round_remove_LFC <- average_otu_tables_200_round_remove_LFC %>%
  dplyr::select(SampleID, everything())

average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_1000_round$SampleID <- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(SampleID, everything())

write.csv(average_otu_tables_5_round, "Data/prokaryotes/prok_average_otu_tables_5.csv", quote=F, row.names=F )
write.csv(average_otu_tables_50_round, "Data/prokaryotes/prok_average_otu_tables_50.csv", quote=F, row.names=F )
write.csv(average_otu_tables_100_round, "Data/prokaryotes/prok_average_otu_tables_100.csv", quote=F, row.names=F )
write.csv(average_otu_tables_200_round, "Data/prokaryotes/prok_average_otu_tables_200.csv", quote=F, row.names=F )
write.csv(average_otu_tables_200_round_remove_LFC, "Data/prokaryotes/prok_average_otu_tables_200_remove_LFC.csv", quote=F, row.names=F )
write.csv(average_otu_tables_1000_round, "/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/prokaryotes/prok_average_otu_tables_1000.csv", quote=F, row.names=F )

# Now check if different number of iterations yield the same pattern
average_otu_tables_5 <- read.csv("Data/prokaryotes/prok_average_otu_tables_5.csv", header=T)
average_otu_tables_5$cover_based_iterations <- "five"
average_otu_tables_5 <- average_otu_tables_5 %>%
  dplyr::select(cover_based_iterations, everything())

average_otu_tables_50 <- read.csv("Data/prokaryotes/prok_average_otu_tables_50.csv", header=T)
average_otu_tables_50$cover_based_iterations <- "fifty"
average_otu_tables_50 <- average_otu_tables_50 %>%
  dplyr::select(cover_based_iterations, everything())

average_otu_tables_100 <- read.csv("Data/prokaryotes/prok_average_otu_tables_100.csv", header=T)
average_otu_tables_100$cover_based_iterations <- "onehundred"
average_otu_tables_100 <- average_otu_tables_100 %>%
  dplyr::select(cover_based_iterations, everything())

average_otu_tables_200 <- read.csv("Data/prokaryotes/prok_average_otu_tables_200.csv", header=T)
average_otu_tables_200$cover_based_iterations <- "twohundred"
average_otu_tables_200 <- average_otu_tables_200 %>%
  dplyr::select(cover_based_iterations, everything())

# average_otu_tables_200_remove_LFC <- read.csv("Data/prokaryotes/prok_average_otu_tables_200_remove_LFC.csv", header=T)
# average_otu_tables_200_remove_LFC$cover_based_iterations <- "twohundred_remove_LFC"
# average_otu_tables_200_remove_LFC <- average_otu_tables_200_remove_LFC %>%
#   dplyr::select(cover_based_iterations, everything())

average_otu_tables_1000 <- read.csv("Data/prokaryotes/prok_average_otu_tables_1000.csv", header=T)
average_otu_tables_1000$cover_based_iterations <- "onethousand"
average_otu_tables_1000 <- average_otu_tables_1000 %>%
  dplyr::select(cover_based_iterations, everything())

# # join all in a single data frame to make boxplot for alpha diversity
otu_tables_iterations <- bind_rows(average_otu_tables_5, average_otu_tables_50,average_otu_tables_100, average_otu_tables_200, average_otu_tables_1000)

### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/prokaryotes/EDITED_16S_final_metadata.csv",header=T )
metadata_sel_iter <- metadata %>% 
  dplyr::select(c(SampleID,swab_id, barcode_plate, barcode_well, year ,region, site, host_species, host_type, sample_type, survey_type, quadrat_id, meso_shoot_id))
master_table_iter <- inner_join(metadata_sel_iter , otu_tables_iterations , by = "SampleID")
View(as.data.frame(master_table_iter ))

############################################################
### Compare alpha and beta diversities across iterations ###
############################################################

library(vegan)
### Creating an object to store abundances only
abundances_16S <- master_table_iter %>% 
  dplyr::select(-(1:14))
names(master_table_iter)[1:17]
View(as.data.frame(abundances_16S))
# Calculate alpha diversity metrics
shannon <- diversity(abundances_16S, index = "shannon")
View(as.data.frame(shannon))
## creating data frame with chao and metadata
alpha_16S <- data.frame(shannon, master_table_iter$region, master_table_iter$year,master_table_iter$cover_based_iterations)

### renaming columns new name = old name
alpha_16S_metrics <- alpha_16S %>%
  dplyr::rename(region= master_table_iter.region, year=  master_table_iter.year, cover_based_iterations =master_table_iter.cover_based_iterations) 

alpha_16S_metrics$year = factor(alpha_16S_metrics$year, levels=c("2015","2016", "2017", "2018"))
alpha_16S_metrics$region <- factor(alpha_16S_metrics$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
alpha_16S_metrics$cover_based_iterations <- factor(alpha_16S_metrics$cover_based_iterations, levels=c("five", "fifty", "onehundred","twohundred", "onethousand"))

boxplot_cover_based_iterations <- ggplot(alpha_16S_metrics, aes(x = cover_based_iterations, y = shannon, fill = year)) + #fill allows to set different colors 
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_brewer(palette = "Dark2")
boxplot_cover_based_iterations
ggsave("R_Code_and_Analysis/alphadiversity/test_cover_based_iterations_16S_COVER_BAS_RAR_years.png", plot = boxplot_cover_based_iterations, width=250, height=200, units="mm",dpi=300)

# Calculate beta diversity
###LOG-transformation *** log1p = log(1+x)
spe.log <- log1p(abundances_16S)
#attach(master_table_iter)
spe.dist<-vegdist(spe.log,method='bray')
spe.dist
beta.spe.host <-betadisper(spe.dist, cover_based_iterations, type = c("median","centroid")) ##always attach so it will remember the "host_type" column
plot(beta.spe.host)

boxplot(beta.spe.host)

##############################################################
### saving Coverage-based rarefaction with more iterations ###
#############################################################
master_table_iter

### Exclude the following samples for analyses: "ZosCSPE", "ZosCSPF" # no info if new or old leaf and ZosCSPoldM and ZosPBSoldD18 which was all NAs
exclude <- c("ZosCSPE", "ZosCSPF")
master_table_iter <- master_table_iter %>% 
  dplyr::filter(!SampleID %in% exclude)

###recode to site names used by grazers
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

# For mastel final table, get only leaf_old 
master_table_iter_final <- master_table_iter %>% 
  dplyr::filter(sample_type =="leaf_old")

# get only meso_quadrat survey 
master_table_iter_final <- master_table_iter_final %>% 
  dplyr::filter(survey_type == "meso_quadrat")

# create a region_year column so can remove only mcmullin 2016 samples
master_table_iter_final <- master_table_iter_final %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

# reorganize column orders (get region_year to first columns together with metadata)
master_table_iter_final <- master_table_iter_final %>%
  dplyr::select(SampleID, swab_id, barcode_plate, barcode_well, year, region_year, everything())

# remove mcmullin 2016 samples
master_table_iter_final <- master_table_iter_final %>% 
  dplyr::filter(!region_year == "mcmullins_2016")

#create a unique site_quadrat_id column
master_table_iter_final <- master_table_iter_final %>% 
  unite(site_quadrat_id, site, quadrat_id, sep = "_" , remove = FALSE) #remove F so it doesn't remove the columns that were combined
View(master_table_iter_final)


master_table_5_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "five")
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
master_table_5_iter <- master_table_5_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_5_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_5_COVERAGE_RAREF.csv", row.names=F)

master_table_50_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "fifty")
master_table_50_iter <- master_table_50_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_50_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_50_COVERAGE_RAREF.csv", row.names=F)

master_table_100_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "onehundred")
master_table_100_iter <- master_table_100_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_100_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_100_COVERAGE_RAREF.csv", row.names=F)

master_table_200_iter <- master_table_iter_final %>% 
  filter(cover_based_iterations == "twohundred")
master_table_200_iter <- master_table_200_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_200_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_200_COVERAGE_RAREF.csv", row.names=F)

master_table_200_remove_LFC_iter <- master_table_iter_final %>% 
  filter(cover_based_iterations == "twohundred_remove_LFC")
master_table_200_remove_LFC_iter <- master_table_200_remove_LFC_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_200_remove_LFC_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_200_remove_LFC_COVERAGE_RAREF.csv", row.names=F)

master_table_1000_iter <- master_table_iter_final %>% 
  filter(cover_based_iterations == "onethousand")
master_table_1000_iter <- master_table_1000_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_1000_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", row.names=F)
