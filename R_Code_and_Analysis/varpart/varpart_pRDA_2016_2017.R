### varpart pRDA 2016 2017 ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 03, 2020 ###

# environmental variables used: temperature, salinity, water depth, leaf area index (lai),  seagrass dry weight (biomass), microepiphyte biomass and bed area 

library(vegan)
library(stats)
library(readr)
library(reshape)
library(dplyr)
library(adespatial)
library(tidyverse)

#########################################
############ 16S prokaryotes ############
#########################################
# 
# ### Read table metadata and abundances
# microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# 
# # keep only 2016 and 2017
# keep_2016_2017 <- c("2016", "2017")
# microbes_16S_ASV <- microbes_16S_ASV %>% 
#   filter( year %in% keep_2016_2017) %>% 
#   arrange(year)
# 
# microbes_16S_ASV <- microbes_16S_ASV %>%
#   dplyr::mutate(year_label = recode(year,
#                                     "2016"="A",
#                                     "2017"="B"))
# 
# microbes_16S_ASV <- microbes_16S_ASV %>% 
#   dplyr::select(SampleID, year_label, everything())
# 
# ### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16" # no abiot data
# ### remove "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17" # no lai, microepiph data
# remove_no_env <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16", "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17", "ZosTQSoldD17")
# microbes_16S_ASV <- microbes_16S_ASV %>% 
#   dplyr::filter(!SampleID %in% remove_no_env )
# 
# #create a unique site_quadrat_id_year column
# microbes_16S_ASV <- microbes_16S_ASV %>%
#   unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
#   select(SampleID, quadrat_year, everything())
# 
# ### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
# environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
# names(environmental_hakai)
# 
# # environmental_hakai <- environmental_hakai %>% 
# #   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))
# 
# prok_quad_envir <- inner_join(environmental_hakai, microbes_16S_ASV, by = "quadrat_year")
# 
# prok_quad_envir <- prok_quad_envir %>% 
#   arrange(site.y, region_year)
# 
# # check sample order to create distance matrix
# sample_order <- as.data.frame(prok_quad_envir$SampleID)
# write.csv(sample_order, "Data/prokaryotes/sample_order_for_geogr_dist_matrix_2016_2017.csv", row.names = F)
# 
# ### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
# prokary_dist<-read.csv("Data/prokaryotes/prokary_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values
# 
# #view sample order to select only 2016 and 2017 from the distance matrix 
# samples_2016_2017 <- c("ZosCFIoldB16", "ZosCFIoldF16", "ZosCIoldA17", "ZosCIoldB17", "ZosCIoldC17", "ZosCIoldD17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldB16", "ZosCSSoldA16","ZosCSSoldG16", "ZosCSSoldA17", "ZosCSSoldB17", "ZosCSSolFF17", "ZosGSWoldF16", "ZosGSWoldD16", "ZosGSWoldA16", "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldH16", "ZosPBSoldA16", "ZosPBSoldB16", "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16", "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldA17", "ZosPBSoldC17", "ZosPBSolEE17", "ZosPBSoldB17", "ZosPBSoldD17", "ZosPBPoldH16", "ZosPBPoldB16", "ZosPBPoldA16", "ZosPBPoldE16", "ZosPBPoldF16", "ZosPBPoldC16", "ZosPBPoldD16", "ZosPBPoldA17", "ZosPBPoldB17", "ZosPBPoldC17", "ZosPBPoldD17", "ZosPBPolEE17", "ZosTQNoldG16", "ZosTQNoldH16", "ZosTQNoldF16", "ZosTQNoldD16", "ZosTQNoldB16", "ZosTQNoldE16", "ZosTQNoldC16", "ZosTQNoldA16", "ZosTQNoldB17", "ZosTQNoldA17", "ZosTQNoldC17", "ZosTQNolEE17", "ZosTQNoldD17", "ZosTQSoldG16", "ZosTQSoldA16", "ZosTQSoldH16", "ZosTQSoldB16", "ZosTQSoldF16", "ZosTQSoldC16", "ZosTQSoldD16", "ZosTQSoldB17", "ZosTQSoldA17", "ZosTQSoldC17")
# 
# # select sometimes bug if using a vector to select columns, so I'll use the whole list of samples I want to keep here
# prokary_dist_16_17 <- prokary_dist %>% 
#   dplyr::select(c("SampleID","ZosCFIoldB16", "ZosCFIoldF16", "ZosCIoldA17", "ZosCIoldB17", "ZosCIoldC17", "ZosCIoldD17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldB16", "ZosCSSoldA16","ZosCSSoldG16", "ZosCSSoldA17", "ZosCSSoldB17", "ZosCSSolFF17", "ZosGSWoldF16", "ZosGSWoldD16", "ZosGSWoldA16", "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldH16", "ZosPBSoldA16", "ZosPBSoldB16", "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16", "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldA17", "ZosPBSoldC17", "ZosPBSolEE17", "ZosPBSoldB17", "ZosPBSoldD17", "ZosPBPoldH16", "ZosPBPoldB16", "ZosPBPoldA16", "ZosPBPoldE16", "ZosPBPoldF16", "ZosPBPoldC16", "ZosPBPoldD16", "ZosPBPoldA17", "ZosPBPoldB17", "ZosPBPoldC17", "ZosPBPoldD17", "ZosPBPolEE17", "ZosTQNoldG16", "ZosTQNoldH16", "ZosTQNoldF16", "ZosTQNoldD16", "ZosTQNoldB16", "ZosTQNoldE16", "ZosTQNoldC16", "ZosTQNoldA16", "ZosTQNoldB17", "ZosTQNoldA17", "ZosTQNoldC17", "ZosTQNolEE17", "ZosTQNoldD17", "ZosTQSoldG16", "ZosTQSoldA16", "ZosTQSoldH16", "ZosTQSoldB16", "ZosTQSoldF16", "ZosTQSoldC16", "ZosTQSoldD16", "ZosTQSoldB17", "ZosTQSoldA17", "ZosTQSoldC17"))
# 
# prokary_dist_16_17 <- prokary_dist_16_17 %>% 
#   dplyr::filter(SampleID %in% samples_2016_2017 )
# nrow(prokary_dist_16_17)
# 
# labels <- as.data.frame(prokary_dist_16_17$SampleID)
# names(labels)[1] <- "SampleID"
# 
# prokary_dist_16_17 <- prokary_dist_16_17 %>% dplyr::select(c(-(SampleID)))
# nrow(prokary_dist_16_17)
# prokary_dist_16_17_matrix <- as.dist(prokary_dist_16_17, diag = TRUE)
# 
# ### CREATE EUCLIDEAN DIST MATRIX 
# distgeo_prokary_16_17 <-vegdist(prokary_dist_16_17_matrix, method="eu")
# distgeo_prokary_16_17_matrix <- as.dist(distgeo_prokary_16_17, diag = TRUE)
# 
# ########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# # putting FINAL TABLE in the same order as the distance matrix
# # remove samples with no env data
# prok_quad_envir <- prok_quad_envir %>% 
#   dplyr::filter(SampleID %in% samples_2016_2017 ) %>% 
#   dplyr::filter(!SampleID %in% remove_no_env)
# 
# prokary_right_order <- left_join(labels, prok_quad_envir, by="SampleID")
# # check if it worked
# label_table <- as.data.frame(prokary_right_order$SampleID)
# sanity_check <- cbind(label_table, labels)
# 
# # ### NOW CAN CREATE DIST MATRIX FOR ZOS
# ### ZOSTERA DIST ###
# # remove metadata 
# names(prokary_right_order)[1:40]
# prokary_no_label <- prokary_right_order %>% dplyr::select(-(1:38))
# 
# #Hellinger pre-transformation of the species data
# prokary_hell <- decostand(prokary_no_label, "hellinger")
# prokary_hell
# 
# ### ENVIRONMENTAL VARIABLES ###
# prokary_env <- prokary_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g, quadrat_macroalgae_g,bed_area_m2)
# 
# ###VIF
# names(prokary_env)
# rda(prokary_hell ~.,prokary_env)
# ###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
# rdacomplete_env<-  rda(formula = prokary_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_biomass_g + quadrat_macroalgae_g + bed_area_m2, data = prokary_env)
# vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again
# 
# # Forward.sel
# sel_env_prokary <-forward.sel(prokary_hell, prokary_env)
# sel_env_prokary
# names(prokary_env)
# prokary_env_sel <- prokary_env %>% dplyr::select(c(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_macroalgae_g, quadrat_biomass_g, bed_area_m2))
# 
# ### prokaryTERA USING PCNMs ###
# pcnm_prokary <- pcnm(distgeo_prokary_16_17)
# pcnm_prokary_vectors <- pcnm_prokary$vectors
# 
# sel_pcnm_prokary<-forward.sel(prokary_hell, pcnm_prokary_vectors)
# sel_pcnm_prokary
# ### In this case, all PCNMs were selected
# pcnm_prokary_vectors <- as.data.frame(pcnm_prokary_vectors)
# pcnm_prokary_selected <- pcnm_prokary_vectors %>% dplyr::select(c(PCNM1, PCNM2, PCNM4, PCNM3))
# 
# ###VARPART WITH YEAR
# # include dummy variable year A = 2015, B = 2016
# prokary_right_order$year_label
# prokary_year_factor <-as.data.frame(prokary_right_order$year_label)
# 
# # varpart
# prokary_varpart_separately_YEARS <- varpart(prokary_hell, pcnm_prokary_selected,prokary_env_sel,prokary_year_factor )
# jpeg('R_Code_and_Analysis/varpart/prokary_varpart_16_17.jpg')
# plot(prokary_varpart_separately_YEARS, digits = 1, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
# title(main = "Prokary_16_17")
# dev.off()
# 
# # Component [a], pure spatial
# anova(rda(prokary_hell, pcnm_prokary_vectors, cbind(prokary_env_sel, year_factor)))
# 
# # Component [b], pure abiotic
# anova(rda(prokary_hell, prokary_env_sel, cbind(pcnm_prokary_vectors, year_factor)))
# 
# # Component [c], pure year
# anova(rda(prokary_hell, year_factor, cbind(pcnm_prokary_vectors, prokary_env_sel)))
# 
# 
# #########################################
# ############ 18S microeukaryotes ############
# #########################################
# 
# ### Read table metadata and abundances
# microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# 
# # keep only 2016 and 2017
# keep_2016_2017 <- c("2016", "2017")
# microbes_18S_ASV <- microbes_18S_ASV %>% 
#   filter( year %in% keep_2016_2017) %>% 
#   arrange(year)
# 
# microbes_18S_ASV <- microbes_18S_ASV %>%
#   dplyr::mutate(year_label = recode(year,
#                                     "2016"="A",
#                                     "2017"="B"))
# 
# microbes_18S_ASV <- microbes_18S_ASV %>% 
#   dplyr::select(SampleID, year_label, everything())
# 
# ### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16" # no abiot data
# ### remove "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17" # no lai, microepiph data
# remove_no_env <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16", "ZosGSEoldG16","ZosGSEoldE16","ZosGSEoldA16","ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17", "ZosTQSoldD17", "ZosTQSolEE17")
# microbes_18S_ASV <- microbes_18S_ASV %>% 
#   dplyr::filter(!SampleID %in% remove_no_env )
# 
# #create a unique site_quadrat_id_year column
# microbes_18S_ASV <- microbes_18S_ASV %>%
#   unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
#   select(SampleID, quadrat_year, everything())
# 
# ### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
# environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
# names(environmental_hakai)
# 
# # environmental_hakai <- environmental_hakai %>% 
# #   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))
# 
# micro_quad_envir <- inner_join(environmental_hakai, microbes_18S_ASV, by = "quadrat_year")
# 
# micro_quad_envir <- micro_quad_envir %>% 
#   arrange(site.y, year.y)
# 
# # check sample order to create distance matrix
# sample_order <- as.data.frame(micro_quad_envir$SampleID)
# write.csv(sample_order, "Data/micro_eukaryotes/microeuk_sample_order_for_geogr_dist_matrix_2016_2017.csv", row.names = F)
# 
# ### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
# microeuk_dist<-read.csv("Data/micro_eukaryotes/microeuk_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values
# 
# #view sample order to select only 2016 and 2017 from the distance matrix 
# samples_2016_2017 <- c("ZosCFIoldB16", "ZosCFIoldG16", "ZosCIoldC17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldA16", "ZosCSSoldG16",  "ZosGSWoldF16", "ZosGSWoldD16",  "ZosGSWoldA16",  "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldE16",  "ZosGSWoldH16",  "ZosPBSoldA16",  "ZosPBSoldB16",  "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16",  "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldC17", "ZosPBPoldH16", "ZosPBPoldB16",  "ZosPBPoldG16", "ZosPBPoldE16",   "ZosPBPoldC16",  "ZosPBPoldD16",   "ZosPBPolEE17",  "ZosTQNoldG16",  "ZosTQNoldH16", "ZosTQNoldF16",  "ZosTQNoldD16",  "ZosTQNoldB16",  "ZosTQNoldE16",  "ZosTQNoldC16",  "ZosTQNoldA16",  "ZosTQNoldB17",  "ZosTQNoldC17",  "ZosTQNoldD17",  "ZosTQSoldG16",  "ZosTQSoldA16",  "ZosTQSoldH16",  "ZosTQSoldB16",  "ZosTQSoldF16",  "ZosTQSoldC16",   "ZosTQSoldD16",  "ZosTQSoldB17", "ZosTQSoldC17")
# # select sometimes bug if using a vector to select columns, so I'll use the whole list of samples I want to keep here
# microeuk_dist_16_17 <- microeuk_dist %>% 
#   dplyr::select(c("SampleID","ZosCFIoldB16", "ZosCFIoldG16", "ZosCIoldC17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldA16", "ZosCSSoldG16",  "ZosGSWoldF16", "ZosGSWoldD16",  "ZosGSWoldA16",  "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldE16",  "ZosGSWoldH16",  "ZosPBSoldA16",  "ZosPBSoldB16",  "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16",  "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldC17", "ZosPBPoldH16", "ZosPBPoldB16",  "ZosPBPoldG16", "ZosPBPoldE16",   "ZosPBPoldC16",  "ZosPBPoldD16",   "ZosPBPolEE17",  "ZosTQNoldG16",  "ZosTQNoldH16", "ZosTQNoldF16",  "ZosTQNoldD16",  "ZosTQNoldB16",  "ZosTQNoldE16",  "ZosTQNoldC16",  "ZosTQNoldA16",  "ZosTQNoldB17",  "ZosTQNoldC17",  "ZosTQNoldD17",  "ZosTQSoldG16",  "ZosTQSoldA16",  "ZosTQSoldH16",  "ZosTQSoldB16",  "ZosTQSoldF16",  "ZosTQSoldC16",   "ZosTQSoldD16",  "ZosTQSoldB17", "ZosTQSoldC17"))
# 
# microeuk_dist_16_17 <- microeuk_dist_16_17 %>% 
#   dplyr::filter(SampleID %in% samples_2016_2017 )
# nrow(microeuk_dist_16_17)
# 
# labels <- as.data.frame(microeuk_dist_16_17$SampleID)
# names(labels)[1] <- "SampleID"
# 
# microeuk_dist_16_17 <- microeuk_dist_16_17 %>% dplyr::select(c(-(SampleID)))
# nrow(microeuk_dist_16_17)
# microeuk_dist_16_17_matrix <- as.dist(microeuk_dist_16_17, diag = TRUE)
# 
# ### CREATE EUCLIDEAN DIST MATRIX 
# distgeo_microeuk_16_17 <-vegdist(microeuk_dist_16_17_matrix, method="eu")
# distgeo_microeuk_16_17_matrix <- as.dist(distgeo_microeuk_16_17, diag = TRUE)
# 
# ########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# # putting FINAL TABLE in the same order as the distance matrix
# # remove samples with no env data
# micro_quad_envir <- micro_quad_envir %>% 
#   dplyr::filter(SampleID %in% samples_2016_2017 ) %>% 
#   dplyr::filter(!SampleID %in% remove_no_env)
# 
# microeuk_right_order <- left_join(labels, micro_quad_envir, by="SampleID")
# # check if it worked
# label_table <- as.data.frame(microeuk_right_order$SampleID)
# sanity_check <- cbind(label_table, labels)
# 
# # ### NOW CAN CREATE DIST MATRIX FOR ZOS
# ### ZOSTERA DIST ###
# # remove metadata 
# names(microeuk_right_order)[1:40]
# microeuk_no_label <- microeuk_right_order %>% dplyr::select(-(1:32))
# 
# #Hellinger pre-transformation of the species data
# microeuk_hell <- decostand(microeuk_no_label, "hellinger")
# microeuk_hell
# 
# ### ENVIRONMENTAL VARIABLES ###
# microeuk_env <- microeuk_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g, quadrat_macroalgae_g,bed_area_m2)
# 
# ###VIF
# names(microeuk_env)
# rda(microeuk_hell ~.,microeuk_env)
# ###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
# rdacomplete_env<-  rda(formula = microeuk_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_biomass_g + quadrat_macroalgae_g + bed_area_m2, data = microeuk_env)
# vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again
# 
# # Forward.sel
# sel_env_microeuk <-forward.sel(microeuk_hell, microeuk_env)
# sel_env_microeuk
# names(microeuk_env)
# microeuk_env_sel <- microeuk_env %>% dplyr::select(c(temperature, salinity, depth, quadrat_microepiphyte_mg,  bed_area_m2))
# 
# ### spatial distance using PCNMs ###
# pcnm_microeuk <- pcnm(distgeo_microeuk_16_17)
# pcnm_microeuk_vectors <- pcnm_microeuk$vectors
# 
# sel_pcnm_microeuk<-forward.sel(microeuk_hell, pcnm_microeuk_vectors)
# sel_pcnm_microeuk
# ### In this case, all PCNMs were selected
# pcnm_microeuk_vectors <- as.data.frame(pcnm_microeuk_vectors)
# pcnm_microeuk_selected <- pcnm_microeuk_vectors %>% dplyr::select(c(PCNM1, PCNM3, PCNM2))
# 
# ###VARPART WITH YEAR
# # include dummy variable year A = 2015, B = 2016
# microeuk_right_order$year_label
# microeuk_year_factor <-as.data.frame(microeuk_right_order$year_label)
# 
# # varpart
# microeuk_varpart_separately_YEARS <- varpart(microeuk_hell, pcnm_microeuk_selected, microeuk_env_sel,microeuk_year_factor )
# jpeg('R_Code_and_Analysis/varpart/microeuk_varpart_16_17.jpg')
# plot(microeuk_varpart_separately_YEARS, digits = 1, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
# title(main = "Microeuk_16_17")
# dev.off()
# 
# # Component [a], pure spatial
# anova(rda(microeuk_hell, pcnm_microeuk_vectors, cbind(microeuk_env_sel, year_factor)))
# 
# # Component [b], pure abiotic
# anova(rda(microeuk_hell, microeuk_env_sel, cbind(pcnm_microeuk_vectors, year_factor)))
# 
# # Component [c], pure year
# anova(rda(microeuk_hell, year_factor, cbind(pcnm_microeuk_vectors, microeuk_env_sel)))

#################################################
############ Inverts Macroeukaryotes ############
#################################################

#### Finest level ####

### Read table metadata and abundances
inverts_finest <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")
names(inverts_finest)[1:12]

# keep only 2016 and 2017
keep_2016_2017 <- c("2016", "2017")
inverts_finest <- inverts_finest %>% 
  filter( year %in% keep_2016_2017) %>% 
  arrange(year)

inverts_finest <- inverts_finest %>%
  dplyr::mutate(year_label = recode(year,
                                    "2016"="A",
                                    "2017"="B"))

inverts_finest <- inverts_finest %>% 
  dplyr::select(ID_year, year_label, everything())


### remove no env data:  "triquet_south_5_2016", "triquet_south_4_2017", "triquet_south_5_2017"
remove_no_env_inverts <- c("triquet_south_5_2016", "triquet_south_4_2017", "triquet_south_5_2017")
inverts_finest <- inverts_finest %>% 
  dplyr::filter(!ID_year %in% remove_no_env_inverts )

#create a unique site_quadrat_id_year column
inverts_finest <- inverts_finest %>%
  dplyr::rename(quadrat_year = "ID_year")

### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

# environmental_hakai <- environmental_hakai %>% 
#   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))

inverts_quad_envir <- inner_join(environmental_hakai, inverts_finest, by = "quadrat_year")

inverts_quad_envir <- inverts_quad_envir %>% 
  arrange(site.y, year.y)

# check sample order to create distance matrix
sample_order <- as.data.frame(inverts_quad_envir$quadrat_year)
write.csv(sample_order, "Data/macro_eukaryotes/macroeuk_sample_order_for_geogr_dist_matrix_2016_2017.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
inverts_dist_16_17<-read.csv("Data/macro_eukaryotes/macroeuk_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values

labels <- as.data.frame(inverts_dist_16_17$SampleID)
names(labels)[1] <- "quadrat_year"

inverts_dist_16_17 <- inverts_dist_16_17 %>% dplyr::select(c(-(SampleID)))
nrow(inverts_dist_16_17)
inverts_dist_16_17_matrix <- as.dist(inverts_dist_16_17, diag = TRUE)

### CREATE EUCLIDEAN DIST MATRIX 
distgeo_inverts_16_17 <-vegdist(inverts_dist_16_17_matrix, method="eu")
distgeo_inverts_16_17_matrix <- as.dist(distgeo_inverts_16_17, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
inverts_right_order <- left_join(labels, inverts_quad_envir, by="quadrat_year")
# check if it worked
label_table <- as.data.frame(inverts_right_order$quadrat_year)
sanity_check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
### ZOSTERA DIST ###
# remove metadata 
names(inverts_right_order)[1:35]
inverts_no_label <- inverts_right_order %>% dplyr::select(-(1:29))

#Hellinger pre-transformation of the species data
inverts_hell <- decostand(inverts_no_label, "hellinger")
inverts_hell

### ENVIRONMENTAL VARIABLES ###
inverts_env <- inverts_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g, quadrat_macroalgae_g,bed_area_m2)

###VIF
names(inverts_env)
rda(inverts_hell ~.,inverts_env)
###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
rdacomplete_env<-  rda(formula = inverts_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_biomass_g + quadrat_macroalgae_g + bed_area_m2, data = inverts_env)
vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again

# Forward.sel
sel_env_inverts <-forward.sel(inverts_hell, inverts_env)
sel_env_inverts
names(inverts_env)
inverts_env_sel <- inverts_env %>% dplyr::select(c(temperature, salinity, depth, quadrat_macroalgae_g, quadrat_biomass_g, bed_area_m2))

### invertsTERA USING PCNMs ###
pcnm_inverts <- pcnm(distgeo_inverts_16_17)
pcnm_inverts_vectors <- pcnm_inverts$vectors

sel_pcnm_inverts<-forward.sel(inverts_hell, pcnm_inverts_vectors)
sel_pcnm_inverts
### In this case, all PCNMs were selected
pcnm_inverts_vectors <- as.data.frame(pcnm_inverts_vectors)
pcnm_inverts_sel <- pcnm_inverts_vectors %>% dplyr::select(c(PCNM3, PCNM2, PCNM1, PCNM4))

###VARPART WITH YEAR
# include dummy variable year A = 2016, B = 2017
inverts_right_order$year_label
inverts_year_factor <-as.data.frame(inverts_right_order$year_label)

# varpart
inverts_varpart_separately_YEARS <- varpart(inverts_hell, pcnm_inverts_sel,inverts_env_sel,inverts_year_factor )
jpeg('R_Code_and_Analysis/varpart/inverts_varpart_16_17.jpg')
plot(inverts_varpart_separately_YEARS, digits = 4, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
title(main = "Inverts_16_17")
dev.off()

inverts_varpart_separately_YEARS

# Component [a], pure spatial
anova(rda(inverts_hell, pcnm_inverts_vectors, cbind(inverts_env_sel, inverts_year_factor)))

# Component [b], pure aenvironmental
anova(rda(inverts_hell, inverts_env_sel, cbind(pcnm_inverts_vectors, inverts_year_factor)))

# Component [c], pure year
anova(rda(inverts_hell, inverts_year_factor, cbind(pcnm_inverts_vectors, inverts_env_sel)))
# 
# #####################################
# ############ Macroeuk18S ############
# #####################################
# 
# ### Read table metadata and abundances
# master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv", header=TRUE)
# 
# # keep only 2016 and 2017
# keep_2016_2017 <- c("2016", "2017")
# master_df_Macro18S <- master_df_Macro18S %>% 
#   filter( year %in% keep_2016_2017) %>% 
#   arrange(year)
# 
# master_df_Macro18S <- master_df_Macro18S %>%
#   dplyr::mutate(year_label = recode(year,
#                                     "2016"="A",
#                                     "2017"="B"))
# 
# master_df_Macro18S <- master_df_Macro18S %>% 
#   dplyr::select(SampleID, year_label, everything())
# 
# ### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16" # no abiot data
# ### remove "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17" # no lai, microepiph data
# remove_no_env <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldE16", "ZosGSEoldA16",  "ZosGSEoldH16" ,"ZosTQSoldE16")
# master_df_Macro18S <- master_df_Macro18S %>% 
#   dplyr::filter(!SampleID %in% remove_no_env )
# 
# #create a unique site_quadrat_id_year column
# master_df_Macro18S <- master_df_Macro18S %>%
#   unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
#   select(SampleID, quadrat_year, everything())
# 
# ### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
# environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
# names(environmental_hakai)
# 
# # environmental_hakai <- environmental_hakai %>% 
# #   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))
# 
# macro18S_quad_envir <- inner_join(environmental_hakai, master_df_Macro18S, by = "quadrat_year")
# 
# macro18S_quad_envir <- macro18S_quad_envir %>% 
#   arrange(site.y, year.y)
# 
# # check sample order to create distance matrix
# sample_order <- as.data.frame(macro18S_quad_envir$SampleID)
# write.csv(sample_order, "Data/macro_18S/macro18S_sample_order_for_geogr_dist_matrix_2016_2017.csv", row.names = F)
# 
# ### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
# ### *** used the same as microeuk because they are derived from the same dataset
# macro18S_dist<-read.csv("Data/macro_18S/macro18S_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values
# 
# #view sample order to select only 2016 and 2017 from the distance matrix 
# samples_2016_2017 <- c("ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldD16","ZosCSSoldB16", "ZosCSSoldA16", "ZosCSSoldG16", "ZosGSWoldF16", "ZosGSWoldD16", "ZosGSWoldA16", "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldE16", "ZosGSWoldH16", "ZosPBSoldA16", "ZosPBSoldB16", "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16", "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBPoldH16", "ZosPBPoldB16", "ZosPBPoldG16", "ZosPBPoldE16", "ZosPBPoldC16", "ZosPBPoldA17", "ZosTQNoldG16", "ZosTQNoldH16", "ZosTQNoldF16", "ZosTQNoldB16", "ZosTQNoldE16", "ZosTQNoldA16", "ZosTQNoldB17", "ZosTQNoldD17", "ZosTQSoldG16", "ZosTQSoldA16", "ZosTQSoldH16", "ZosTQSoldB16", "ZosTQSoldF16", "ZosTQSoldC16", "ZosTQSoldD16")
# # select sometimes bug if using a vector to select columns, so I'll use the whole list of samples I want to keep here
# macro18S_dist_16_17 <- macro18S_dist %>% 
#   dplyr::select(c("SampleID", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldD16","ZosCSSoldB16", "ZosCSSoldA16", "ZosCSSoldG16", "ZosGSWoldF16", "ZosGSWoldD16", "ZosGSWoldA16", "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldE16", "ZosGSWoldH16", "ZosPBSoldA16", "ZosPBSoldB16", "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16", "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBPoldH16", "ZosPBPoldB16", "ZosPBPoldG16", "ZosPBPoldE16", "ZosPBPoldC16", "ZosPBPoldA17", "ZosTQNoldG16", "ZosTQNoldH16", "ZosTQNoldF16", "ZosTQNoldB16", "ZosTQNoldE16", "ZosTQNoldA16", "ZosTQNoldB17", "ZosTQNoldD17", "ZosTQSoldG16", "ZosTQSoldA16", "ZosTQSoldH16", "ZosTQSoldB16", "ZosTQSoldF16", "ZosTQSoldC16", "ZosTQSoldD16"))
# 
# macro18S_dist_16_17 <- macro18S_dist_16_17 %>% 
#   dplyr::filter(SampleID %in% samples_2016_2017 )
# nrow(macro18S_dist_16_17)
# 
# labels <- as.data.frame(macro18S_dist_16_17$SampleID)
# names(labels)[1] <- "SampleID"
# 
# macro18S_dist_16_17 <- macro18S_dist_16_17 %>% dplyr::select(c(-(SampleID)))
# nrow(macro18S_dist_16_17)
# macro18S_dist_16_17_matrix <- as.dist(macro18S_dist_16_17, diag = TRUE)
# 
# ### CREATE EUCLIDEAN DIST MATRIX 
# distgeo_macro18S_16_17 <-vegdist(macro18S_dist_16_17_matrix, method="eu")
# distgeo_macro18S_16_17_matrix <- as.dist(distgeo_macro18S_16_17, diag = TRUE)
# 
# ########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# # putting FINAL TABLE in the same order as the distance matrix
# # remove samples with no env data
# macro18S_quad_envir <- macro18S_quad_envir %>% 
#   dplyr::filter(SampleID %in% samples_2016_2017 ) %>% 
#   dplyr::filter(!SampleID %in% remove_no_env)
# 
# macro18S_right_order <- left_join(labels, macro18S_quad_envir, by="SampleID")
# # check if it worked
# label_table <- as.data.frame(macro18S_right_order$SampleID)
# sanity_check <- cbind(label_table, labels)
# 
# # ### NOW CAN CREATE DIST MATRIX FOR ZOS
# ### ZOSTERA DIST ###
# # remove metadata 
# names(macro18S_right_order)[1:40]
# macro18S_no_label <- macro18S_right_order %>% dplyr::select(-(1:32))
# 
# #Hellinger pre-transformation of the species data
# macro18S_hell <- decostand(macro18S_no_label, "hellinger")
# macro18S_hell
# 
# ### ENVIRONMENTAL VARIABLES ###
# macro18S_env <- macro18S_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g, quadrat_macroalgae_g,bed_area_m2)
# 
# ###VIF
# names(macro18S_env)
# rda(macro18S_hell ~.,macro18S_env)
# ###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
# rdacomplete_env<-  rda(formula = macro18S_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_biomass_g + quadrat_macroalgae_g + bed_area_m2, data = macro18S_env)
# vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again
# 
# # Forward.sel
# sel_env_macro18S <-forward.sel(macro18S_hell, macro18S_env)
# sel_env_macro18S
# macro18S_env_sel <- macro18S_env %>% dplyr::select(c(temperature, salinity, depth, quadrat_microepiphyte_mg, quadrat_biomass_g, bed_area_m2))
# 
# ### spatial distance using PCNMs ###
# pcnm_macro18S <- pcnm(distgeo_macro18S_16_17)
# pcnm_macro18S_vectors <- pcnm_macro18S$vectors
# 
# sel_pcnm_macro18S<-forward.sel(macro18S_hell, pcnm_macro18S_vectors)
# sel_pcnm_macro18S
# ### In this case, all PCNMs were selected
# pcnm_macro18S_vectors <- as.data.frame(pcnm_macro18S_vectors)
# pcnm_macro18S_selected <- pcnm_macro18S_vectors %>% dplyr::select(c(PCNM1, PCNM3))
# 
# ###VARPART WITH YEAR
# # include dummy variable year A = 2015, B = 2016
# macro18S_right_order$year_label
# macro18S_year_factor <-as.data.frame(macro18S_right_order$year_label)
# 
# # varpart
# macro18S_varpart_separately_YEARS <- varpart(macro18S_hell, pcnm_macro18S_selected, macro18S_env_sel,macro18S_year_factor )
# jpeg('R_Code_and_Analysis/varpart/macro18S_varpart_16_17.jpg')
# plot(macro18S_varpart_separately_YEARS, digits = 1, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
# title(main = "macro18S_varpart_16_17")
# dev.off()
# 
# # Component [a], pure spatial
# anova(rda(macro18S_hell, pcnm_macro18S_vectors, cbind(macro18S_env_sel, macro18S_year_factor)))
# 
# # Component [b], pure abiotic
# anova(rda(macro18S_hell, macro18S_env_sel, cbind(pcnm_macro18S_vectors, macro18S_year_factor)))
# 
# # Component [c], pure year
# anova(rda(macro18S_hell, macro18S_year_factor, cbind(pcnm_macro18S_vectors, macro18S_env_sel)))

##############################################################
################ EXTRA MICROBES AT GENUS LEVEL ###############
##############################################################


#########################################
############ 16S prokaryotes ############
#########################################

### Read table metadata and abundances
microbes_16S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv", header=T)

# keep only 2016 and 2017
keep_2016_2017 <- c("2016", "2017")
microbes_16S_genus <- microbes_16S_genus %>% 
  filter( year %in% keep_2016_2017) %>% 
  arrange(year)

microbes_16S_genus <- microbes_16S_genus %>%
  dplyr::mutate(year_label = recode(year,
                                    "2016"="A",
                                    "2017"="B"))

microbes_16S_genus <- microbes_16S_genus %>% 
  dplyr::select(SampleID, year_label, everything())

### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16" # no abiot data
### remove "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17" # no lai, microepiph data
remove_no_env <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16", "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17", "ZosTQSoldD17")
microbes_16S_genus <- microbes_16S_genus %>% 
  dplyr::filter(!SampleID %in% remove_no_env )

#create a unique site_quadrat_id_year column
microbes_16S_genus <- microbes_16S_genus %>%
  unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
  select(SampleID, quadrat_year, everything())

### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

# environmental_hakai <- environmental_hakai %>% 
#   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))

prok_quad_envir <- inner_join(environmental_hakai, microbes_16S_genus, by = "quadrat_year")

prok_quad_envir <- prok_quad_envir %>% 
  arrange(site.y, region_year)

# check sample order to create distance matrix
sample_order <- as.data.frame(prok_quad_envir$SampleID)
write.csv(sample_order, "Data/prokaryotes/sample_order_for_geogr_dist_matrix_2016_2017.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
prokary_dist<-read.csv("Data/prokaryotes/prokary_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values

#view sample order to select only 2016 and 2017 from the distance matrix 
samples_2016_2017 <- c("ZosCFIoldB16", "ZosCFIoldF16", "ZosCIoldA17", "ZosCIoldB17", "ZosCIoldC17", "ZosCIoldD17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldB16", "ZosCSSoldA16","ZosCSSoldG16", "ZosCSSoldA17", "ZosCSSoldB17", "ZosCSSolFF17", "ZosGSWoldF16", "ZosGSWoldD16", "ZosGSWoldA16", "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldH16", "ZosPBSoldA16", "ZosPBSoldB16", "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16", "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldA17", "ZosPBSoldC17", "ZosPBSolEE17", "ZosPBSoldB17", "ZosPBSoldD17", "ZosPBPoldH16", "ZosPBPoldB16", "ZosPBPoldA16", "ZosPBPoldE16", "ZosPBPoldF16", "ZosPBPoldC16", "ZosPBPoldD16", "ZosPBPoldA17", "ZosPBPoldB17", "ZosPBPoldC17", "ZosPBPoldD17", "ZosPBPolEE17", "ZosTQNoldG16", "ZosTQNoldH16", "ZosTQNoldF16", "ZosTQNoldD16", "ZosTQNoldB16", "ZosTQNoldE16", "ZosTQNoldC16", "ZosTQNoldA16", "ZosTQNoldB17", "ZosTQNoldA17", "ZosTQNoldC17", "ZosTQNolEE17", "ZosTQNoldD17", "ZosTQSoldG16", "ZosTQSoldA16", "ZosTQSoldH16", "ZosTQSoldB16", "ZosTQSoldF16", "ZosTQSoldC16", "ZosTQSoldD16", "ZosTQSoldB17", "ZosTQSoldA17", "ZosTQSoldC17")

# select sometimes bug if using a vector to select columns, so I'll use the whole list of samples I want to keep here
prokary_dist_16_17 <- prokary_dist %>% 
  dplyr::select(c("SampleID","ZosCFIoldB16", "ZosCFIoldF16", "ZosCIoldA17", "ZosCIoldB17", "ZosCIoldC17", "ZosCIoldD17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldB16", "ZosCSSoldA16","ZosCSSoldG16", "ZosCSSoldA17", "ZosCSSoldB17", "ZosCSSolFF17", "ZosGSWoldF16", "ZosGSWoldD16", "ZosGSWoldA16", "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldH16", "ZosPBSoldA16", "ZosPBSoldB16", "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16", "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldA17", "ZosPBSoldC17", "ZosPBSolEE17", "ZosPBSoldB17", "ZosPBSoldD17", "ZosPBPoldH16", "ZosPBPoldB16", "ZosPBPoldA16", "ZosPBPoldE16", "ZosPBPoldF16", "ZosPBPoldC16", "ZosPBPoldD16", "ZosPBPoldA17", "ZosPBPoldB17", "ZosPBPoldC17", "ZosPBPoldD17", "ZosPBPolEE17", "ZosTQNoldG16", "ZosTQNoldH16", "ZosTQNoldF16", "ZosTQNoldD16", "ZosTQNoldB16", "ZosTQNoldE16", "ZosTQNoldC16", "ZosTQNoldA16", "ZosTQNoldB17", "ZosTQNoldA17", "ZosTQNoldC17", "ZosTQNolEE17", "ZosTQNoldD17", "ZosTQSoldG16", "ZosTQSoldA16", "ZosTQSoldH16", "ZosTQSoldB16", "ZosTQSoldF16", "ZosTQSoldC16", "ZosTQSoldD16", "ZosTQSoldB17", "ZosTQSoldA17", "ZosTQSoldC17"))

prokary_dist_16_17 <- prokary_dist_16_17 %>% 
  dplyr::filter(SampleID %in% samples_2016_2017 )
nrow(prokary_dist_16_17)

labels <- as.data.frame(prokary_dist_16_17$SampleID)
names(labels)[1] <- "SampleID"

prokary_dist_16_17 <- prokary_dist_16_17 %>% dplyr::select(c(-(SampleID)))
nrow(prokary_dist_16_17)
prokary_dist_16_17_matrix <- as.dist(prokary_dist_16_17, diag = TRUE)

### CREATE EUCLIDEAN DIST MATRIX 
distgeo_prokary_16_17 <-vegdist(prokary_dist_16_17_matrix, method="eu")
distgeo_prokary_16_17_matrix <- as.dist(distgeo_prokary_16_17, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
# remove samples with no env data
prok_quad_envir <- prok_quad_envir %>% 
  dplyr::filter(SampleID %in% samples_2016_2017 ) %>% 
  dplyr::filter(!SampleID %in% remove_no_env)

prokary_right_order <- left_join(labels, prok_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(prokary_right_order$SampleID)
sanity_check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
### ZOSTERA DIST ###
# remove metadata 
names(prokary_right_order)[1:40]
prokary_no_label <- prokary_right_order %>% dplyr::select(-(1:38))

#Hellinger pre-transformation of the species data
prokary_hell <- decostand(prokary_no_label, "hellinger")
prokary_hell

### ENVIRONMENTAL VARIABLES ###
prokary_env <- prokary_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g, quadrat_macroalgae_g,bed_area_m2)

###VIF
names(prokary_env)
rda(prokary_hell ~.,prokary_env)
###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
rdacomplete_env<-  rda(formula = prokary_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_biomass_g + quadrat_macroalgae_g + bed_area_m2, data = prokary_env)
vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again

# Forward.sel
sel_env_prokary <-forward.sel(prokary_hell, prokary_env)
sel_env_prokary
names(prokary_env)
prokary_env_sel <- prokary_env %>% dplyr::select(c(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_macroalgae_g, quadrat_biomass_g, bed_area_m2))

### prokaryTERA USING PCNMs ###
pcnm_prokary <- pcnm(distgeo_prokary_16_17)
pcnm_prokary_vectors <- pcnm_prokary$vectors

sel_pcnm_prokary<-forward.sel(prokary_hell, pcnm_prokary_vectors)
sel_pcnm_prokary
### In this case, all PCNMs were selected
pcnm_prokary_vectors <- as.data.frame(pcnm_prokary_vectors)
pcnm_prokary_selected <- pcnm_prokary_vectors %>% dplyr::select(c(PCNM1, PCNM4, PCNM2, PCNM3,  PCNM10))

###VARPART WITH YEAR
# include dummy variable year A = 2015, B = 2016
prokary_right_order$year_label
prokary_year_factor <-as.data.frame(prokary_right_order$year_label)

# varpart
prokary_varpart_separately_YEARS_genus <- varpart(prokary_hell, pcnm_prokary_selected,prokary_env_sel,prokary_year_factor )
jpeg('R_Code_and_Analysis/varpart/GENUS_prokary_varpart_16_17.jpg')
plot(prokary_varpart_separately_YEARS_genus, digits = 4, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
title(main = "Prokary_GENUS_16_17")
dev.off()

# Component [a], pure spatial
anova(rda(prokary_hell, pcnm_prokary_vectors, cbind(prokary_env_sel, prokary_year_factor)))

# Component [b], pure abiotic
anova(rda(prokary_hell, prokary_env_sel, cbind(pcnm_prokary_vectors, prokary_year_factor)))

# Component [c], pure year
anova(rda(prokary_hell, prokary_year_factor, cbind(pcnm_prokary_vectors, prokary_env_sel)))


################################################
############ 18S microeukaryotes GENUS #########
################################################

### Read table metadata and abundances
microbes_18S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv", header=T)

# keep only 2016 and 2017
keep_2016_2017 <- c("2016", "2017")
microbes_18S_genus <- microbes_18S_genus %>% 
  filter( year %in% keep_2016_2017) %>% 
  arrange(year)

microbes_18S_genus <- microbes_18S_genus %>%
  dplyr::mutate(year_label = recode(year,
                                    "2016"="A",
                                    "2017"="B"))

microbes_18S_genus <- microbes_18S_genus %>% 
  dplyr::select(SampleID, year_label, everything())

### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16" # no abiot data
### remove "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17" # no lai, microepiph data
remove_no_env <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16", "ZosGSEoldG16","ZosGSEoldE16","ZosGSEoldA16","ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17", "ZosTQSoldD17", "ZosTQSolEE17")
microbes_18S_genus <- microbes_18S_genus %>% 
  dplyr::filter(!SampleID %in% remove_no_env )

#create a unique site_quadrat_id_year column
microbes_18S_genus <- microbes_18S_genus %>%
  unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
  select(SampleID, quadrat_year, everything())

### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

# environmental_hakai <- environmental_hakai %>% 
#   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))

micro_quad_envir <- inner_join(environmental_hakai, microbes_18S_genus, by = "quadrat_year")

micro_quad_envir <- micro_quad_envir %>% 
  arrange(site.y, year.y)

# check sample order to create distance matrix
sample_order <- as.data.frame(micro_quad_envir$SampleID)
write.csv(sample_order, "Data/micro_eukaryotes/microeuk_sample_order_for_geogr_dist_matrix_2016_2017.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
microeuk_dist<-read.csv("Data/micro_eukaryotes/microeuk_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values

#view sample order to select only 2016 and 2017 from the distance matrix 
samples_2016_2017 <- c("ZosCFIoldB16", "ZosCFIoldG16", "ZosCIoldC17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldA16", "ZosCSSoldG16",  "ZosGSWoldF16", "ZosGSWoldD16",  "ZosGSWoldA16",  "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldE16",  "ZosGSWoldH16",  "ZosPBSoldA16",  "ZosPBSoldB16",  "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16",  "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldC17", "ZosPBPoldH16", "ZosPBPoldB16",  "ZosPBPoldG16", "ZosPBPoldE16",   "ZosPBPoldC16",  "ZosPBPoldD16",   "ZosPBPolEE17",  "ZosTQNoldG16",  "ZosTQNoldH16", "ZosTQNoldF16",  "ZosTQNoldD16",  "ZosTQNoldB16",  "ZosTQNoldE16",  "ZosTQNoldC16",  "ZosTQNoldA16",  "ZosTQNoldB17",  "ZosTQNoldC17",  "ZosTQNoldD17",  "ZosTQSoldG16",  "ZosTQSoldA16",  "ZosTQSoldH16",  "ZosTQSoldB16",  "ZosTQSoldF16",  "ZosTQSoldC16",   "ZosTQSoldD16",  "ZosTQSoldB17", "ZosTQSoldC17")
# select sometimes bug if using a vector to select columns, so I'll use the whole list of samples I want to keep here
microeuk_dist_16_17 <- microeuk_dist %>% 
  dplyr::select(c("SampleID","ZosCFIoldB16", "ZosCFIoldG16", "ZosCIoldC17", "ZosCSSoldC16", "ZosCSSoldF16", "ZosCSSoldE16", "ZosCSSoldA16", "ZosCSSoldG16",  "ZosGSWoldF16", "ZosGSWoldD16",  "ZosGSWoldA16",  "ZosGSWoldB16", "ZosGSWoldG16", "ZosGSWoldC16", "ZosGSWoldE16",  "ZosGSWoldH16",  "ZosPBSoldA16",  "ZosPBSoldB16",  "ZosPBSoldC16", "ZosPBSoldD16", "ZosPBSoldG16", "ZosPBSoldF16",  "ZosPBSoldH16", "ZosPBSoldE16", "ZosPBSoldC17", "ZosPBPoldH16", "ZosPBPoldB16",  "ZosPBPoldG16", "ZosPBPoldE16",   "ZosPBPoldC16",  "ZosPBPoldD16",   "ZosPBPolEE17",  "ZosTQNoldG16",  "ZosTQNoldH16", "ZosTQNoldF16",  "ZosTQNoldD16",  "ZosTQNoldB16",  "ZosTQNoldE16",  "ZosTQNoldC16",  "ZosTQNoldA16",  "ZosTQNoldB17",  "ZosTQNoldC17",  "ZosTQNoldD17",  "ZosTQSoldG16",  "ZosTQSoldA16",  "ZosTQSoldH16",  "ZosTQSoldB16",  "ZosTQSoldF16",  "ZosTQSoldC16",   "ZosTQSoldD16",  "ZosTQSoldB17", "ZosTQSoldC17"))

microeuk_dist_16_17 <- microeuk_dist_16_17 %>% 
  dplyr::filter(SampleID %in% samples_2016_2017 )
nrow(microeuk_dist_16_17)

labels <- as.data.frame(microeuk_dist_16_17$SampleID)
names(labels)[1] <- "SampleID"

microeuk_dist_16_17 <- microeuk_dist_16_17 %>% dplyr::select(c(-(SampleID)))
nrow(microeuk_dist_16_17)
microeuk_dist_16_17_matrix <- as.dist(microeuk_dist_16_17, diag = TRUE)

### CREATE EUCLIDEAN DIST MATRIX 
distgeo_microeuk_16_17 <-vegdist(microeuk_dist_16_17_matrix, method="eu")
distgeo_microeuk_16_17_matrix <- as.dist(distgeo_microeuk_16_17, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
# remove samples with no env data
micro_quad_envir <- micro_quad_envir %>% 
  dplyr::filter(SampleID %in% samples_2016_2017) 

microeuk_right_order <- left_join(labels, micro_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(microeuk_right_order$SampleID)
sanity_check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
### ZOSTERA DIST ###
# remove metadata 
names(microeuk_right_order)[1:40]
microeuk_no_label <- microeuk_right_order %>% dplyr::select(-(1:32))

#Hellinger pre-transformation of the species data
microeuk_hell <- decostand(microeuk_no_label, "hellinger")
microeuk_hell

### ENVIRONMENTAL VARIABLES ###
microeuk_env <- microeuk_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g, quadrat_macroalgae_g,bed_area_m2)

###VIF
names(microeuk_env)
rda(microeuk_hell ~.,microeuk_env)
###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
rdacomplete_env<-  rda(formula = microeuk_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_biomass_g + quadrat_macroalgae_g + bed_area_m2, data = microeuk_env)
vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again

# Forward.sel
sel_env_microeuk <-forward.sel(microeuk_hell, microeuk_env)
sel_env_microeuk
names(microeuk_env)
microeuk_env_sel <- microeuk_env %>% dplyr::select(c(temperature, salinity, quadrat_microepiphyte_mg, bed_area_m2))

### spatial distance using PCNMs ###
pcnm_microeuk <- pcnm(distgeo_microeuk_16_17)
pcnm_microeuk_vectors <- pcnm_microeuk$vectors

sel_pcnm_microeuk<-forward.sel(microeuk_hell, pcnm_microeuk_vectors)
sel_pcnm_microeuk
### In this case, all PCNMs were selected
pcnm_microeuk_vectors <- as.data.frame(pcnm_microeuk_vectors)
pcnm_microeuk_selected <- pcnm_microeuk_vectors %>% dplyr::select(c(PCNM1, PCNM3, PCNM2))

###VARPART WITH YEAR
# include dummy variable year A = 2015, B = 2016
microeuk_right_order$year_label
microeuk_year_factor <-as.data.frame(microeuk_right_order$year_label)

# varpart
microeuk_varpart_separately_YEARS_GENUS <- varpart(microeuk_hell, pcnm_microeuk_selected, microeuk_env_sel,microeuk_year_factor )
jpeg('R_Code_and_Analysis/varpart/GENUS_microeuk_varpart_16_17.jpg')
plot(microeuk_varpart_separately_YEARS_GENUS, digits = 4, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
title(main = "Microeuk_GENUS_16_17")
dev.off()

# Component [a], pure spatial
anova(rda(microeuk_hell, pcnm_microeuk_vectors, cbind(microeuk_env_sel, microeuk_year_factor)))

# Component [b], pure abiotic
anova(rda(microeuk_hell, microeuk_env_sel, cbind(pcnm_microeuk_vectors, microeuk_year_factor)))

# Component [c], pure year
anova(rda(microeuk_hell, microeuk_year_factor, cbind(pcnm_microeuk_vectors, microeuk_env_sel)))

#################################################
############### STACKED HISTOGRAM ###############
#################################################
varpart_results <- read.csv("R_Code_and_Analysis/varpart/varpart_results_stacked_histo.csv", header = TRUE)

varpart_results$host <- factor(varpart_results$host, levels = c("prokaryotes", "microeukaryotes", "macroeukaryotes"))

varpart_results$component <- factor(varpart_results$component, levels = c("pure_spa", "pure_env", "shared_env_spa", "pure_year", "shared_env_year", "shared_spa_year" ,"shared_env_spa_year", "residuals"))

colour_portion <- c("darkorange1", "green3", "brown3", "darkorchid3",  "blue3", "gold1", "maroon1", "grey36")

varpart_histo <- ggplot(varpart_results,aes(host,proportion_explained,fill=component))+
  geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
  facet_grid(. ~ host, drop=TRUE,scale="free",space="free_x")

varpart_histo <- varpart_histo + ggtitle("Variation Partitioning")

varpart_histo <- varpart_histo + scale_fill_manual(values=colour_portion)

varpart_histo <- varpart_histo + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))

varpart_histo <- varpart_histo + 
  theme (strip.text.x = element_text(size = 11),
         axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size = 12), #change font size of numbers
                axis.title.y=element_text(size = 14), #change font size of y title
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                legend.text=element_text(size=12),
                legend.text.align=0,
                plot.title = element_text(hjust = 0.5, size = 16, face="bold")) 

varpart_histo

ggsave("R_Code_and_Analysis/varpart/varpart_histogram_results.png", plot = varpart_histo, width=180, height=150, units="mm",dpi=300)
