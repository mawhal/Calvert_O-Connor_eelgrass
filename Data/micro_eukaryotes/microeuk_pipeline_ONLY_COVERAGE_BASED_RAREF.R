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

all_years_18S_filtered_meso_Zos

##################################################
### rarefying data using coverage based iNEXT
##################################################
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
library(metagMisc)

# phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

# Samples standardized by size will have different degrees of completness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).

# i.e. a seagrass sample with 10,000 total reads will have a different coverage than a seawater sample with 10,000 reads if seawater samples have many more species

# phyloseq object to be used
all_years_18S_filtered_meso_Zos

### the following steps are necessary before running the function because it needs the otu table in the phyloseq object to have taxa as rows

# check if taxa are rows in the phyloseq object
taxa_are_rows(all_years_18S_filtered_meso_Zos)

# transpose so taxa are rows
otu_table(all_years_18S_filtered_meso_Zos) <- t(otu_table(all_years_18S_filtered_meso_Zos))

taxa_are_rows(all_years_18S_filtered_meso_Zos)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(all_years_18S_filtered_meso_Zos)),
  correct_singletons = T)

# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

#Due to the stochasticity introduced in random subsampling results could be slightly different.
#So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.

# run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)

all_18S_COVERAGE_RAREF_5 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_meso_Zos, coverage = 0.8, iter = 5, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_5, "Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_5.rds")
all_18S_COVERAGE_RAREF_5 <- readRDS("Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_5.rds")

all_18S_COVERAGE_RAREF_50 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_meso_Zos, coverage = 0.8, iter = 50, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_50, "Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_50.rds")
all_18S_COVERAGE_RAREF_50 <- readRDS("Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_50.rds")

all_18S_COVERAGE_RAREF_100 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_meso_Zos, coverage = 0.8, iter = 100, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_5=100, "Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_100.rds")
all_18S_COVERAGE_RAREF_100 <- readRDS("Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_100.rds")

all_18S_COVERAGE_RAREF_200 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_meso_Zos, coverage = 0.8, iter = 200, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_5=200, "Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_200.rds")
all_18S_COVERAGE_RAREF_200 <- readRDS("Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_200.rds")

all_18S_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_meso_Zos, coverage = 0.8, iter = 1000, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_1000, "Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_1000.rds")
all_18S_COVERAGE_RAREF_1000 <- readRDS("Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_5 <- all_18S_COVERAGE_RAREF_5[c(1:5)]
subset_phylo_objects_50 <- all_18S_COVERAGE_RAREF_50[c(1:50)]
subset_phylo_objects_100 <- all_18S_COVERAGE_RAREF_100[c(1:100)]
subset_phylo_objects_200 <- all_18S_COVERAGE_RAREF_200[c(1:200)]
subset_phylo_objects_1000 <- all_18S_COVERAGE_RAREF_1000[c(1:1000)]
# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
otu_tables_5 <- lapply(subset_phylo_objects_5, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_50 <- lapply(subset_phylo_objects_50, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_100 <- lapply(subset_phylo_objects_100, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_200 <- lapply(subset_phylo_objects_200, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# average all matrices to get the mean abundance across all iterations
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
average_otu_tables_100_round <- average_otu_tables_100 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
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

average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_1000_round$SampleID <- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(SampleID, everything())

write.csv(average_otu_tables_5_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_5.csv", quote=F, row.names=F )
write.csv(average_otu_tables_50_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_50.csv", quote=F, row.names=F )
write.csv(average_otu_tables_100_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_100.csv", quote=F, row.names=F )
write.csv(average_otu_tables_200_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_200.csv", quote=F, row.names=F )
write.csv(average_otu_tables_1000_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_1000.csv", quote=F, row.names=F )

# Now check if different number of iterations yield the same pattern
average_otu_tables_5 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_5.csv", header=T)
average_otu_tables_5$cover_based_iterations <- "five"
average_otu_tables_5 <- average_otu_tables_5 %>%
  select(cover_based_iterations, everything())

average_otu_tables_50 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_50.csv", header=T)
average_otu_tables_50$cover_based_iterations <- "fifty"
average_otu_tables_50 <- average_otu_tables_50 %>%
  select(cover_based_iterations, everything())

average_otu_tables_100 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_100.csv", header=T)
average_otu_tables_100$cover_based_iterations <- "onehundred"
average_otu_tables_100 <- average_otu_tables_100 %>%
  select(cover_based_iterations, everything())

average_otu_tables_200 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_200.csv", header=T)
average_otu_tables_200$cover_based_iterations <- "twohundred"
average_otu_tables_200 <- average_otu_tables_200 %>%
  select(cover_based_iterations, everything())

average_otu_tables_1000 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_1000.csv", header=T)
average_otu_tables_1000$cover_based_iterations <- "onethousand"
average_otu_tables_1000 <- average_otu_tables_1000 %>%
  select(cover_based_iterations, everything())

# # join all in a single data frame to make boxplot for alpha diversity
otu_tables_iterations <- bind_rows(average_otu_tables_5, average_otu_tables_50,average_otu_tables_100, average_otu_tables_200, average_otu_tables_1000)

# # ### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/micro_eukaryotes/18S_ASV_level_metadata_table.csv",header=T )
# View(metadata)
master_table_iter <- inner_join(metadata , otu_tables_iterations , by = "SampleID")
View(as.data.frame(master_table_iter ))

# # ############################################################
# # ### Compare alpha and beta diversities across iterations ###
# # ############################################################
# #
library(vegan)
# ### Creating an object to store abundances only
abundances_18S <- master_table_iter %>%
  dplyr::select(-(1:9))
names(master_table_iter)
# # Calculate alpha diversity metrics
shannon <- diversity(abundances_18S, index = "shannon")

# ## creating data frame with chao and metadata
alpha_18S <- data.frame(shannon, master_table_iter$region, master_table_iter$year,master_table_iter$cover_based_iterations)

# ### renaming columns new name = old name
alpha_18S_metrics <- alpha_18S %>%
  dplyr::rename(region= master_table_iter.region, year=  master_table_iter.year, cover_based_iterations =master_table_iter.cover_based_iterations)

alpha_18S_metrics$year = factor(alpha_18S_metrics$year, levels=c("2015","2016", "2017", "2018"))
alpha_18S_metrics$region <- factor(alpha_18S_metrics$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
alpha_18S_metrics$cover_based_iterations <- factor(alpha_18S_metrics$cover_based_iterations, levels=c("five", "fifty", "onehundred","twohundred", "onethousand"))

boxplot_cover_based_iterations <- ggplot(alpha_18S_metrics, aes(x = cover_based_iterations, y = shannon, fill = year)) + #fill allows to set different colors
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_brewer(palette = "Dark2")
boxplot_cover_based_iterations
ggsave("R_Code_and_Analysis/alphadiversity/test_iterations_18S_COVER_BAS_RAR_years.png", plot = boxplot_cover_based_iterations, width=250, height=200, units="mm",dpi=300)

# # Calculate beta diversity
# ###LOG-transformation *** log1p = log(1+x)
spe.log <- log1p(abundances_18S)
attach(master_table_iter)
spe.dist<-vegdist(spe.log,method='bray')
spe.dist
beta.spe.host <-betadisper(spe.dist, cover_based_iterations, type = c("median","centroid")) ##always attach so it will remember the "host_type" column
plot(beta.spe.host)

boxplot(beta.spe.host)
#=
# ##############################################################
# ### saving Coverage-based rarefaction with more iterations ###
# #############################################################
#
master_table_iter
#
# ### Exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
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


master_table_5_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "five")
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
master_table_5_iter <- master_table_5_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_5_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_5_COVERAGE_RAREF.csv", row.names=F)

master_table_50_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "fifty")
master_table_50_iter <- master_table_50_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_50_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_50_COVERAGE_RAREF.csv", row.names=F)

master_table_100_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "onehundred")
master_table_100_iter <- master_table_100_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_100_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_100_COVERAGE_RAREF.csv", row.names=F)

master_table_200_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "twohundred")
master_table_200_iter <- master_table_200_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_200_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_200_COVERAGE_RAREF.csv", row.names=F)

master_table_1000_iter <- master_table_iter_final %>%
  filter(cover_based_iterations == "onethousand")
master_table_1000_iter <- master_table_1000_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
write.csv(master_table_1000_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", row.names=F)
