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

library(devtools)
library(metagMisc)
library(iNEXT)

#### Importing files ####
all_years_18S_filtered_for_CB <- readRDS("/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/micro_eukaryotes/all_years_18S_filtered_for_CB.rds")

all_18S_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=all_years_18S_filtered_for_CB, coverage = 0.8, iter = 1000, replace = F, correct_singletons = TRUE, drop_lowcoverage = F)
saveRDS(all_18S_COVERAGE_RAREF_1000, "/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_1000.rds")
all_18S_COVERAGE_RAREF_1000 <- readRDS("/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/micro_eukaryotes/all_18S_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
# subset_phylo_objects_5 <- all_18S_COVERAGE_RAREF_5[c(1:5)]
# subset_phylo_objects_50 <- all_18S_COVERAGE_RAREF_50[c(1:50)]
# subset_phylo_objects_100 <- all_18S_COVERAGE_RAREF_100[c(1:100)]
# subset_phylo_objects_200 <- all_18S_COVERAGE_RAREF_200[c(1:200)]
subset_phylo_objects_1000 <- all_18S_COVERAGE_RAREF_1000[c(1:1000)]
# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
# otu_tables_5 <- lapply(subset_phylo_objects_5, function(z) as.data.frame(t(phyloseq::otu_table(z))))
# otu_tables_50 <- lapply(subset_phylo_objects_50, function(z) as.data.frame(t(phyloseq::otu_table(z))))
# otu_tables_100 <- lapply(subset_phylo_objects_100, function(z) as.data.frame(t(phyloseq::otu_table(z))))
# otu_tables_200 <- lapply(subset_phylo_objects_200, function(z) as.data.frame(t(phyloseq::otu_table(z))))
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# average all matrices to get the mean abundance across all iterations
# average_otu_tables_5 <- Reduce("+",otu_tables_5)/length(otu_tables_5)
# # IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
# average_otu_tables_5_round <- average_otu_tables_5 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # add SampleID column back
# average_otu_tables_5_round$SampleID <- rownames(average_otu_tables_5)
# average_otu_tables_5_round <- average_otu_tables_5_round %>%
#   select(SampleID, everything())
# 
# average_otu_tables_50 <- Reduce("+",otu_tables_50)/length(otu_tables_50)
# # IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
# average_otu_tables_50_round <- average_otu_tables_50 %>% dplyr::mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # add SampleID column back
# average_otu_tables_50_round$SampleID <- rownames(average_otu_tables_50)
# average_otu_tables_50_round <- average_otu_tables_50_round %>%
#   select(SampleID, everything())
# 
# average_otu_tables_100 <- Reduce("+",otu_tables_100)/length(otu_tables_100)
# # IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
# average_otu_tables_100_round <- average_otu_tables_100 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # add SampleID column back
# average_otu_tables_100_round$SampleID <- rownames(average_otu_tables_100)
# average_otu_tables_100_round <- average_otu_tables_100_round %>%
#   select(SampleID, everything())
# 
# average_otu_tables_200 <- Reduce("+",otu_tables_200)/length(otu_tables_200)
# # IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
# average_otu_tables_200_round <- average_otu_tables_200 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # add SampleID column back
# average_otu_tables_200_round$SampleID <- rownames(average_otu_tables_200)
# average_otu_tables_200_round <- average_otu_tables_200_round %>%
#   select(SampleID, everything())

average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_1000_round$SampleID <- rownames(average_otu_tables_1000) 
average_otu_tables_1000_round <- average_otu_tables_1000_round %>% 
  select(SampleID, everything())
# 
# write.csv(average_otu_tables_5_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_5.csv", quote=F, row.names=F )
# write.csv(average_otu_tables_50_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_50.csv", quote=F, row.names=F )
# write.csv(average_otu_tables_100_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_100.csv", quote=F, row.names=F )
# write.csv(average_otu_tables_200_round, "Data/micro_eukaryotes/microeuk_average_otu_tables_200.csv", quote=F, row.names=F )
write.csv(average_otu_tables_1000_round, "/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/micro_eukaryotes/microeuk_average_otu_tables_1000.csv", quote=F, row.names=F )

# Now check if different number of iterations yield the same pattern
# average_otu_tables_5 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_5.csv", header=T)
# average_otu_tables_5$cover_based_iterations <- "five"
# average_otu_tables_5 <- average_otu_tables_5 %>%
#   select(cover_based_iterations, everything())
# 
# average_otu_tables_50 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_50.csv", header=T)
# average_otu_tables_50$cover_based_iterations <- "fifty"
# average_otu_tables_50 <- average_otu_tables_50 %>%
#   select(cover_based_iterations, everything())
# 
# average_otu_tables_100 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_100.csv", header=T)
# average_otu_tables_100$cover_based_iterations <- "onehundred"
# average_otu_tables_100 <- average_otu_tables_100 %>%
#   select(cover_based_iterations, everything())
# 
# average_otu_tables_200 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_200.csv", header=T)
# average_otu_tables_200$cover_based_iterations <- "twohundred"
# average_otu_tables_200 <- average_otu_tables_200 %>%
#   select(cover_based_iterations, everything())

average_otu_tables_1000 <- read.csv("Data/micro_eukaryotes/microeuk_average_otu_tables_1000.csv", header=T)
average_otu_tables_1000$cover_based_iterations <- "onethousand"
average_otu_tables_1000 <- average_otu_tables_1000 %>%
  select(cover_based_iterations, everything())

# # join all in a single data frame to make boxplot for alpha diversity
otu_tables_iterations <- bind_rows(average_otu_tables_5, average_otu_tables_50,average_otu_tables_100, average_otu_tables_200)

# # ### add metadata according to #SampleID labels
metadata <- read.csv(file="Data/micro_eukaryotes/18S_ASV_level_metadata_table.csv",header=T )
# View(metadata)
master_table_iter <- inner_join(metadata , otu_tables_iterations , by = "SampleID")
View(as.data.frame(master_table_iter ))
# #
# # ############################################################
# # ### Compare alpha and beta diversities across iterations ###
# # ############################################################
# #
# # library(vegan)
# # ### Creating an object to store abundances only
# # abundances_18S <- master_table_iter %>%
# #   dplyr::select(-(1:9))
# # names(master_table_iter)
# # # Calculate alpha diversity metrics
# # shannon <- diversity(abundances_18S, index = "shannon")
# #
# # ## creating data frame with chao and metadata
# # alpha_18S <- data.frame(shannon, master_table_iter$region, master_table_iter$year,master_table_iter$cover_based_iterations)
# #
# # ### renaming columns new name = old name
# # alpha_18S_metrics <- alpha_18S %>%
# #   dplyr::rename(region= master_table_iter.region, year=  master_table_iter.year, cover_based_iterations =master_table_iter.cover_based_iterations)
# #
# # alpha_18S_metrics$year = factor(alpha_18S_metrics$year, levels=c("2015","2016", "2017", "2018"))
# # alpha_18S_metrics$region <- factor(alpha_18S_metrics$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# # alpha_18S_metrics$cover_based_iterations <- factor(alpha_18S_metrics$cover_based_iterations, levels=c("five", "fifty", "onehundred","twohundred"))
# #
# # boxplot_cover_based_iterations <- ggplot(alpha_18S_metrics, aes(x = cover_based_iterations, y = shannon, fill = year)) + #fill allows to set different colors
# #   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
# #   scale_fill_brewer(palette = "Dark2")
# # boxplot_cover_based_iterations
# # ggsave("R_Code_and_Analysis/alphadiversity/test_iterations_18S_COVER_BAS_RAR_years.png", plot = boxplot_cover_based_iterations, width=250, height=200, units="mm",dpi=300)
# #
# # # Calculate beta diversity
# # ###LOG-transformation *** log1p = log(1+x)
# # spe.log <- log1p(abundances_18S)
# # attach(master_table_iter)
# # spe.dist<-vegdist(spe.log,method='bray')
# # spe.dist
# # beta.spe.host <-betadisper(spe.dist, cover_based_iterations, type = c("median","centroid")) ##always attach so it will remember the "host_type" column
# # plot(beta.spe.host)
# #
# # boxplot(beta.spe.host)
# #
# # ##############################################################
# # ### saving Coverage-based rarefaction with more iterations ###
# # #############################################################
# #
# # master_table_iter
# #
# # ### Exclude the following samples for analyses: ZosCSPtrans3Amb3 (not meso_quadrat but wrongly assigned as such)
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
# #
# # master_table_5_iter <- master_table_iter_final %>%
# #   filter(cover_based_iterations == "five")
# # # IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
# # master_table_5_iter <- master_table_5_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # write.csv(master_table_5_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_5_COVERAGE_RAREF.csv", row.names=F)
# #
# # master_table_50_iter <- master_table_iter_final %>%
# #   filter(cover_based_iterations == "fifty")
# # master_table_50_iter <- master_table_50_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # write.csv(master_table_50_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_50_COVERAGE_RAREF.csv", row.names=F)
# #
# # master_table_100_iter <- master_table_iter_final %>%
# #   filter(cover_based_iterations == "onehundred")
# # master_table_100_iter <- master_table_100_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # write.csv(master_table_100_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_100_COVERAGE_RAREF.csv", row.names=F)
# #
# # master_table_200_iter <- master_table_iter_final %>%
# #   filter(cover_based_iterations == "twohundred")
# # master_table_200_iter <- master_table_200_iter %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# # write.csv(master_table_200_iter, "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_200_COVERAGE_RAREF.csv", row.names=F)
