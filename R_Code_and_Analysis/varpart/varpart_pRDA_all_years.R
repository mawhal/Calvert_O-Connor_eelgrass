### varpart pRDA ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 03, 2020 ###

# environmental variables used: temperature, salinity, dissolved oxygen, water depth, leaf area index (lai),  seagrass dry weight (biomass) and microepiphyte biomass 

library(vegan)
library(stats)
library(readr)
library(reshape)
library(dplyr)
library(adespatial)
library(tidyverse)

#########################################
############ 16S prokaryotes GENUS ############
#########################################

### Read table metadata and abundances
microbes_16S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv", header=T)

microbes_16S_genus <- microbes_16S_genus %>%
  dplyr::mutate(year_label = recode(year,
                                    "2015"="A",
                                    "2016"="B",
                                    "2017"="C",
                                    "2018" = "D"))

microbes_16S_genus <- microbes_16S_genus %>% 
  dplyr::select(SampleID, year_label, everything())

### remove ZosCSPoldA, 	ZosCSPoldF, ZosCSPoldG, ZosCSPoldH, ZosCSPoldL, ZosCSPoldM
#### no info on quadrat it
remove <- c("ZosCSPoldA", "ZosCSPoldF", "ZosCSPoldG", "ZosCSPoldH", "ZosCSPoldL", "ZosCSPoldM")
microbes_16S_genus <- microbes_16S_genus %>% 
  filter(!SampleID %in% remove)

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
write.csv(sample_order, "Data/prokaryotes/sample_order_for_geogr_dist_matrix.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
# This assumes your distance file is in a tab-delimited table format
# also that it has a header, and that the first column contains the sample names (sites)
prokary_dist<-read.csv("Data/prokaryotes/prokary_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values

### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16" # no abiot data
### remove "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17" # no lai, microepiph data
prokary_dist <- prokary_dist %>% 
  dplyr::select(-c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16", "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17"))

remove_gse16_no_abiot <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldC16", "ZosGSEoldB16", "ZosGSEoldH16", "ZosTQSoldE16", "ZosCSSoldH16", "ZosCSSoldC17", "ZosCSSoldD17")
prokary_dist <- prokary_dist %>% 
  dplyr::filter(!SampleID %in% remove_gse16_no_abiot )
nrow(prokary_dist)

labels <- as.data.frame(prokary_dist$SampleID)
names(labels)[1] <- "SampleID"

prokary_dist <- prokary_dist %>% dplyr::select(c(-(SampleID)))
nrow(prokary_dist)
prokary_dist_matrix <- as.dist(prokary_dist, diag = TRUE)

### CREATE EUCLIDEAN DIST MATRIX 
distgeo_prokary <-vegdist(prokary_dist_matrix, method="eu")
distgeo_prokary_matrix <- as.dist(distgeo_prokary, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
# remove samples with no env data
prok_quad_envir <- prok_quad_envir %>% 
  dplyr::filter(!SampleID %in% remove_gse16_no_abiot )

prokary_right_order <- left_join(labels, prok_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(prokary_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
### ZOSTERA DIST ###
# remove metadata 
names(prokary_right_order)[1:40]
prokary_no_label <- prokary_right_order %>% dplyr::select(-(1:38))

#Hellinger pre-transformation of the species data
prokary_hell <- decostand(prokary_no_label, "hellinger")
prokary_hell

### ENVIRONMENTAL VARIABLES ###
prokary_env <- prokary_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_macroalgae_g,bed_area_m2)

# zos_abiot <- final_zos_right_order %>% dplyr::select(depth, temperature, dissolved_oxygen, salinity)
# 
# zos_host <- final_zos_right_order %>% dplyr::select(quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g)

###VIF
names(prokary_env)
rda(prokary_hell ~.,prokary_env)
###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
rdacomplete_env<-  rda(formula = prokary_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_macroalgae_g + bed_area_m2, data = prokary_env)
vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again

# Forward.sel
sel_env_prokary <-forward.sel(prokary_hell, prokary_env)
sel_env_prokary
names(prokary_env)
prokary_env_sel <- prokary_env %>% dplyr::select(c(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_macroalgae_g, bed_area_m2))

### prokaryTERA USING PCNMs ###
pcnm_prokary <- pcnm(distgeo_prokary)
pcnm_prokary_vectors <- pcnm_prokary$vectors

sel_pcnm_prokary<-forward.sel(prokary_hell, pcnm_prokary_vectors)
sel_pcnm_prokary
### In this case, all PCNMs were selected
pcnm_prokary_vectors <- as.data.frame(pcnm_prokary_vectors)
pcnm_prokary_selected <- pcnm_prokary_vectors %>% dplyr::select(c(PCNM1, PCNM2, PCNM5, PCNM3, PCNM4, PCNM7))

###VARPART WITH YEAR
# include dummy variable year A = 2015, B = 2016
prokary_right_order$year_label
prokary_year_factor <-as.data.frame(prokary_right_order$year_label)

# varpart
prokary_varpart_separately_YEARS <- varpart(prokary_hell, pcnm_prokary_selected,prokary_env_sel,prokary_year_factor )
jpeg('R_Code_and_Analysis/varpart/GENUS_prokary_varpart_all_years.jpg')
plot(prokary_varpart_separately_YEARS, digits = 1, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
title(main = "Prokary_all_years")
dev.off()

# Component [a], pure spatial
anova(rda(prokary_hell, pcnm_prokary_selected, cbind(prokary_env_sel, prokary_year_factor)))

# Component [b], pure abiotic
anova(rda(prokary_hell, prokary_env_sel, cbind(pcnm_prokary_selected, prokary_year_factor)))

# Component [c], pure year
anova(rda(prokary_hell, prokary_year_factor, cbind(pcnm_prokary_selected, prokary_env_sel)))

#########################################
############ 18S microeukaryotes GENUS ########
#########################################

### Read table metadata and abundances
microbes_18S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv", header=T)

microbes_18S_genus <- microbes_18S_genus %>%
  dplyr::mutate(year_label = recode(year,
                                    "2015"="A",
                                    "2016"="B",
                                    "2017"="C",
                                    "2018" = "D"))

microbes_18S_genus <- microbes_18S_genus %>% 
  dplyr::select(SampleID, year_label, everything())

### remove ZosCSPoldA, 	ZosCSPoldF, ZosCSPoldM
#### no info on quadrat it
remove <- c("ZosCSPoldA", "ZosCSPoldF", "ZosCSPoldM", "ZosPBPoldA17")
microbes_18S_genus <- microbes_18S_genus %>% 
  filter(!SampleID %in% remove)

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

check_NAS_in_env_data <- micro_quad_envir %>% 
  select(c("SampleID", "depth", "temperature", "salinity", "quadrat_lai", "quadrat_microepiphyte_mg", "quadrat_macroalgae_g", "bed_area_m2"))

# check sample order to create distance matrix
sample_order_micro <- as.data.frame(micro_quad_envir$SampleID)
write.csv(sample_order_micro, "Data/micro_eukaryotes/microeuk_sample_order_for_geogr_dist_matrix.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
# This assumes your distance file is in a tab-delimited table format
# also that it has a header, and that the first column contains the sample names (sites)
microeuk_dist<-read.csv("Data/micro_eukaryotes/microeuk_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values

### remove  (("ZosCSSoldH16","ZosCSSoldD17","ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldG16", "ZosGSEoldE16", "ZosGSEoldC16", "ZosGSEolddA16", "ZosGSEolddB16",  "ZosGSEoldH16", "ZosCSSoldH16" ,"ZosCSSoldD17","ZosTQSoldE16" # no env data

microeuk_dist <- microeuk_dist %>% 
  dplyr::select(-c("ZosCSSoldH16","ZosCSSoldD17","ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldG16", "ZosGSEoldE16", "ZosGSEoldC16", "ZosGSEoldA16", "ZosGSEoldB16",  "ZosGSEoldH16", "ZosCSSoldH16" ,"ZosCSSoldD17","ZosTQSoldE16", "ZosPBPoldA17"))
                
remove_gse16_no_abiot <- c("ZosCSSoldH16","ZosCSSoldD17","ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldG16", "ZosGSEoldE16", "ZosGSEoldC16", "ZosGSEoldA16", "ZosGSEoldB16",  "ZosGSEoldH16", "ZosCSSoldH16" ,"ZosCSSoldD17","ZosTQSoldE16", "ZosPBPoldA17")

microeuk_dist <- microeuk_dist %>% 
  dplyr::filter(!SampleID %in% remove_gse16_no_abiot )
nrow(microeuk_dist)

labels <- as.data.frame(microeuk_dist$SampleID)
names(labels)[1] <- "SampleID"

microeuk_dist <- microeuk_dist %>% dplyr::select(c(-(SampleID)))
nrow(microeuk_dist)
microeuk_dist_matrix <- as.dist(microeuk_dist, diag = TRUE)

### CREATE EUCLIDEAN DIST MATRIX 
distgeo_microeuk <-vegdist(microeuk_dist_matrix, method="eu")
distgeo_microeuk_matrix <- as.dist(distgeo_microeuk, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
# remove samples with no env data
micro_quad_envir <- micro_quad_envir %>% 
  dplyr::filter(!SampleID %in% remove_gse16_no_abiot )

microeuk_right_order <- left_join(labels, micro_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(microeuk_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
### ZOSTERA DIST ###
# remove metadata 
names(microeuk_right_order)[1:40]
microeuk_no_label <- microeuk_right_order %>% dplyr::select(-(1:32))

#Hellinger pre-transformation of the species data
microeuk_hell <- decostand(microeuk_no_label, "hellinger")
microeuk_hell

### ENVIRONMENTAL VARIABLES ###
microeuk_env <- microeuk_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_macroalgae_g,bed_area_m2)

# zos_abiot <- final_zos_right_order %>% dplyr::select(depth, temperature, dissolved_oxygen, salinity)
# 
# zos_host <- final_zos_right_order %>% dplyr::select(quadrat_lai, quadrat_microepiphyte_mg, quadrat_biomass_g)

###VIF
names(microeuk_env)
rda(microeuk_hell ~.,microeuk_env)
###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
rdacomplete_env<-  rda(formula = microeuk_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_macroalgae_g + bed_area_m2, data = microeuk_env)
vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again

# Forward.sel
sel_env_microeuk <-forward.sel(microeuk_hell, microeuk_env)
sel_env_microeuk
names(microeuk_env)
microeuk_env_sel <- microeuk_env %>% dplyr::select(c(temperature, salinity, quadrat_microepiphyte_mg, quadrat_macroalgae_g, bed_area_m2))

### microeukTERA USING PCNMs ###
pcnm_microeuk <- pcnm(distgeo_microeuk)
pcnm_microeuk_vectors <- pcnm_microeuk$vectors

sel_pcnm_microeuk<-forward.sel(microeuk_hell, pcnm_microeuk_vectors)
sel_pcnm_microeuk
### In this case, all PCNMs were selected
pcnm_microeuk_vectors <- as.data.frame(pcnm_microeuk_vectors)
pcnm_microeuk_selected <- pcnm_microeuk_vectors %>% dplyr::select(c(PCNM1, PCNM2, PCNM4, PCNM5, PCNM7))

###VARPART WITH YEAR
# include dummy variable year A = 2015, B = 2016
microeuk_right_order$year_label
microeuk_year_factor <-as.data.frame(microeuk_right_order$year_label)

# varpart
microeuk_varpart_separately_YEARS_GENUS <- varpart(microeuk_hell, pcnm_microeuk_selected, microeuk_env_sel,microeuk_year_factor )
jpeg('R_Code_and_Analysis/varpart/GENUS_microeuk_varpart_all_years.jpg')
plot(microeuk_varpart_separately_YEARS_GENUS, digits = 1, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
title(main = "Microeuk_all_years")
dev.off()

# Component [a], pure spatial
anova(rda(microeuk_hell, pcnm_microeuk_vectors, cbind(microeuk_env_sel, microeuk_year_factor)))

# Component [b], pure abiotic
anova(rda(microeuk_hell, microeuk_env_sel, cbind(pcnm_microeuk_vectors, microeuk_year_factor)))

# Component [c], pure year
anova(rda(microeuk_hell, microeuk_year_factor, cbind(pcnm_microeuk_vectors, microeuk_env_sel)))


#################################
############ Macroeuk18S ########
#################################

### Read table metadata and abundances
master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_genus_level_1000_COR_SING_COVERAGE_RAREF.csv", header=TRUE)

master_df_Macro18S <- master_df_Macro18S %>%
  dplyr::mutate(year_label = recode(year,
                                    "2015"="A",
                                    "2016"="B",
                                    "2017"="C",
                                    "2018" = "D"))

master_df_Macro18S <- master_df_Macro18S %>% 
  dplyr::select(SampleID, year_label, everything())

### remove ZosCSPoldA, 	ZosCSPoldF, ZosCSPoldM
#### no info on quadrat it
remove <- c("ZosCSPoldA", "ZosCSPoldF", "ZosCSPoldM")
master_df_Macro18S <- master_df_Macro18S %>% 
  filter(!SampleID %in% remove)

#create a unique site_quadrat_id_year column
master_df_Macro18S <- master_df_Macro18S %>%
  unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
  select(SampleID, quadrat_year, everything())


### corrected data, plus fixed choked sandspit 2016 depth from zero to 4
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

# environmental_hakai <- environmental_hakai %>% 
#   dplyr::select(-c( quadrat_year, temperature, salinity, depth, dissolved_oxygen, quadrat_lai,quadrat_biomass_g, quadrat_microepiphyte_mg))

macro18S_quad_envir <- inner_join(environmental_hakai, master_df_Macro18S, by = "quadrat_year")

macro18S_quad_envir <- macro18S_quad_envir %>% 
  arrange(site.y, year.y)

check_NAS_in_env_data <- macro18S_quad_envir %>% 
  select(c("SampleID", "depth", "temperature", "salinity", "quadrat_lai", "quadrat_microepiphyte_mg", "quadrat_macroalgae_g", "bed_area_m2"))

# check sample order to create distance matrix
sample_order_macro18S <- as.data.frame(macro18S_quad_envir$SampleID)
write.csv(sample_order_macro18S, "Data/macro_18S/macro18S_sample_order_for_geogr_dist_matrix.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
# This assumes your distance file is in a tab-delimited table format
# also that it has a header, and that the first column contains the sample names (sites)
macro18S_dist<-read.csv("Data/macro_18S/macro18S_geogr_dist_matrix.csv", header=TRUE) # had to put NA in all blank values

### remove  "ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldE16", "ZosGSEoldA16",  "ZosGSEoldH16" ,"ZosTQSoldE16" # no env data

macro18S_dist <- macro18S_dist %>% 
  dplyr::select(-c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldE16", "ZosGSEoldA16",  "ZosGSEoldH16" ,"ZosTQSoldE16")) 

remove_env <- c("ZosGSEoldD16", "ZosGSEoldF16", "ZosGSEoldE16", "ZosGSEoldA16",  "ZosGSEoldH16","ZosTQSoldE16")

macro18S_dist <- macro18S_dist %>% 
  dplyr::filter(!SampleID %in% remove_env )
nrow(macro18S_dist)

labels <- as.data.frame(macro18S_dist$SampleID)
names(labels)[1] <- "SampleID"

macro18S_dist <- macro18S_dist %>% dplyr::select(c(-(SampleID)))
nrow(macro18S_dist)
macro18S_dist_matrix <- as.dist(macro18S_dist, diag = TRUE)

### CREATE EUCLIDEAN DIST MATRIX 
distgeo_macro18S <-vegdist(macro18S_dist_matrix, method="eu")
distgeo_macro18S_matrix <- as.dist(distgeo_macro18S, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
# remove samples with no env data
macro18S_quad_envir <- macro18S_quad_envir %>% 
  dplyr::filter(!SampleID %in% remove_gse16_no_abiot )

macro18S_right_order <- left_join(labels, macro18S_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(macro18S_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
### ZOSTERA DIST ###
# remove metadata 
names(macro18S_right_order)[1:40]
macro18S_no_label <- macro18S_right_order %>% dplyr::select(-(1:38))

#Hellinger pre-transformation of the species data
macro18S_hell <- decostand(macro18S_no_label, "hellinger")
macro18S_hell

### ENVIRONMENTAL VARIABLES ###
macro18S_env <- macro18S_right_order %>% dplyr::select(depth, temperature, salinity, quadrat_lai, quadrat_microepiphyte_mg, quadrat_macroalgae_g, bed_area_m2)

# zos_abiot <- final_zos_right_order %>% dplyr::select(depth, temperature, dissolved_oxygen, salinity)
# 
# zos_host <- final_zos_right_order %>% dplyr::select(quadrat_lai, quadrat_macro18Sepiphyte_mg, quadrat_biomass_g)

###VIF
names(macro18S_env)
rda(macro18S_hell ~.,macro18S_env)
###Aqui faz a RDA normal pra poder rodar o Vif..colocar os nomes das suas vari?veis...
rdacomplete_env<-  rda(formula = macro18S_hell ~ depth + temperature + salinity + quadrat_lai + quadrat_microepiphyte_mg + quadrat_macroalgae_g + bed_area_m2, data = macro18S_env)
vif.cca(rdacomplete_env) ###test for multicollinearity = If any value is higher than 10, remove this variable and run again

# Forward.sel
sel_env_macro18S <-forward.sel(macro18S_hell, macro18S_env)
sel_env_macro18S
names(macro18S_env)
macro18S_env_sel <- macro18S_env %>% dplyr::select(c(depth, temperature, salinity, quadrat_macroalgae_g, bed_area_m2))

### macro18STERA USING PCNMs ###
pcnm_macro18S <- pcnm(distgeo_macro18S)
pcnm_macro18S_vectors <- pcnm_macro18S$vectors

sel_pcnm_macro18S<-forward.sel(macro18S_hell, pcnm_macro18S_vectors)
sel_pcnm_macro18S
### In this case, all PCNMs were selected
pcnm_macro18S_vectors <- as.data.frame(pcnm_macro18S_vectors)
pcnm_macro18S_selected <- pcnm_macro18S_vectors %>% dplyr::select(c(PCNM1, PCNM2, PCNM3, PCNM5, PCNM7, PCNM4))

###VARPART WITH YEAR
# include dummy variable year A = 2015, B = 2016
macro18S_right_order$year_label
macro18S_year_factor <-as.data.frame(macro18S_right_order$year_label)

# varpart
macro18S_varpart_separately_YEARS <- varpart(macro18S_hell, pcnm_macro18S_selected, macro18S_env_sel,macro18S_year_factor )
jpeg('R_Code_and_Analysis/varpart/macro18S_varpart_all_years.jpg')
plot(macro18S_varpart_separately_YEARS, digits = 1, Xnames = c('space', 'env', "year"), bg = c('navy', 'tomato', "orange"))
title(main = "Macro18S_all_years")
dev.off()

# Component [a], pure spatial
anova(rda(macro18S_hell, pcnm_macro18S_vectors, cbind(macro18S_env_sel, year_factor)))

# Component [b], pure abiotic
anova(rda(macro18S_hell, macro18S_env_sel, cbind(pcnm_macro18S_vectors, year_factor)))

# Component [c], pure year
anova(rda(macro18S_hell, year_factor, cbind(pcnm_macro18S_vectors, macro18S_env_sel)))

