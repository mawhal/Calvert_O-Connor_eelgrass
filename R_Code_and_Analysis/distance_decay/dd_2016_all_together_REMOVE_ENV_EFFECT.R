### Distance Decay of Similarity 2016 ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 10th, 2020 ###

library(vegan)
library(stats)
library(ggplot2)
library(readr)
library(reshape)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(emmeans)

#########################################
############ 16S prokaryotes GENUS ######
#########################################
### Read table metadata and abundances
microbes_16S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv", header=T)
attach(microbes_16S_genus)
names(microbes_16S_genus)

microbes_16S_genus_2016 <- filter(microbes_16S_genus, year == "2016") %>% 
  arrange(site)

samples_env_data_prokary_2016 <- c("ZosCFIoldB16","ZosCFIoldF16","ZosCSSoldC16","ZosCSSoldF16","ZosCSSoldE16",
"ZosCSSoldB16","ZosCSSoldA16","ZosCSSoldG16","ZosGSWoldF16","ZosGSWoldD16","ZosGSWoldA16","ZosGSWoldB16","ZosGSWoldG16","ZosGSWoldC16","ZosGSWoldH16","ZosPBSoldA16","ZosPBSoldB16","ZosPBSoldC16","ZosPBSoldD16","ZosPBSoldG16","ZosPBSoldF16","ZosPBSoldH16","ZosPBSoldE16","ZosPBPoldH16","ZosPBPoldB16","ZosPBPoldA16","ZosPBPoldE16","ZosPBPoldF16","ZosPBPoldC16","ZosPBPoldD16","ZosTQNoldG16","ZosTQNoldH16","ZosTQNoldF16","ZosTQNoldD16","ZosTQNoldB16","ZosTQNoldE16","ZosTQNoldC16","ZosTQNoldA16","ZosTQSoldG16","ZosTQSoldA16","ZosTQSoldH16","ZosTQSoldB16","ZosTQSoldF16","ZosTQSoldC16","ZosTQSoldD16")
microbes_16S_genus_2016 <- microbes_16S_genus_2016 %>% 
  dplyr::filter(SampleID %in% samples_env_data_prokary_2016 )%>% 
  arrange(site)

microbes_16S_genus_2016 <- microbes_16S_genus_2016 %>%
  unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
  select(SampleID, quadrat_year, everything())

prokary_2016_samples <- as.data.frame(microbes_16S_genus_2016$SampleID)
write.csv(prokary_2016_samples, "Data/prokaryotes/distance_decay_prokary_samples_2016_REMOVE_ENV.csv", row.names = F)

### Add environmental variables
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

prok_quad_envir <- inner_join(environmental_hakai, microbes_16S_genus_2016, by = "quadrat_year")

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
prokary_dist_2016<-read.csv("Data/prokaryotes/prokary_geogr_dist_matrix_2016_REMOVE_ENV.csv", header=TRUE) 

labels <- as.data.frame(prokary_dist_2016$SampleID)
names(labels)[1] <- "SampleID"

prokary_dist2_2016 <- dplyr::select(prokary_dist_2016, c(-(SampleID)))
nrow(prokary_dist2_2016)

prokary_matrix_2016 <- as.dist(prokary_dist2_2016, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
final_prokary_right_order <- left_join(labels, prok_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(final_prokary_right_order$SampleID)
check <- cbind(label_table, labels)

### get only abundance
names(final_prokary_right_order)[1:40]
prokary_abund_2016 <- final_prokary_right_order %>% dplyr::select(starts_with("ASV")) 
nrow(prokary_abund_2016)

### CREATE DISTANCE MATRIX FOR SPECIES
spe.dist.prokary_2016 <-vegdist(prokary_abund_2016,method="bray")
spe.dist.prokary_2016
spe.dist.prokary_2016.matrix <- as.matrix(spe.dist.prokary_2016)
spe.dist.matrix_test_df <- as.data.frame(spe.dist.prokary_2016.matrix)

### get environmental data
names(final_prokary_right_order)[1:40]
prokary_env_2016 <- final_prokary_right_order %>% dplyr::select(c("depth", "temperature", "salinity", "quadrat_lai", "quadrat_microepiphyte_mg", "quadrat_biomass_g", "quadrat_macroalgae_g", "bed_area_m2"))
prokary_env_2016 <- prokary_env_2016

### CREATE DISTANCE MATRIX FOR ENVIRONMENT
env.dist.prokary_2016 <-vegdist(log(prokary_env_2016+1),method="euclidean")
env.dist.prokary_2016
env.dist.prokary_2016.matrix <- as.matrix(env.dist.prokary_2016)
env.dist.matrix_test_df <- as.data.frame(env.dist.prokary_2016.matrix)

### get vector from species and environmental data
df.spe.prokary_2016 <- data.frame(dissimilarity=spe.dist.prokary_2016[upper.tri(spe.dist.prokary_2016, diag = FALSE)])
df.spe.prokary_2016 #get row number that contain values to remove NAs
df.spe.prokary.vector_2016 <- df.spe.prokary_2016[1:231,]

df.env.prokary_2016 <- data.frame(dissimilarity=env.dist.prokary_2016[upper.tri(env.dist.prokary_2016, diag = FALSE)])
df.env.prokary_2016 #get row number that contain values to remove NAs
df.env.prokary.vector_2016 <- df.env.prokary_2016[1:231,]


###ggplot only work with data frames, so we need to convert this data into data frame
library(reshape2)
df.spe.prokary.vector_2016 <- melt(as.matrix(spe.dist.prokary_2016), varnames = c("row", "col"))
names(df.spe.prokary.vector_2016)[names(df.spe.prokary.vector_2016) == 'value'] <- 'dissimilarity.prokary'
df.spe.prokary.vector_2016
dissimilarity_new_prokary_2016 <- df.spe.prokary.vector_2016$dissimilarity.prokary

df.spa.prokary.vector_2016 <- melt(as.matrix(prokary_matrix_2016), varnames = c("row", "col"))
names(df.spa.prokary.vector_2016)[names(df.spa.prokary.vector_2016) == 'value'] <- 'spatial_distance.prokary'
df.spa.prokary.vector_2016
spatial_prokary_2016 <- df.spa.prokary.vector_2016$spatial_distance.prokary

df.env.prokary.vector_2016 <- melt(as.matrix(env.dist.prokary_2016), varnames = c("row", "col"))
names(df.env.prokary.vector_2016)[names(df.env.prokary.vector_2016) == 'value'] <- 'environmental_distance.prokary'
df.env.prokary.vector_2016
environment_prokary_2016 <- df.env.prokary.vector_2016$env_distance.prokary

as.data.frame(df.spe.prokary.vector_2016)
as.data.frame(df.spa.prokary.vector_2016)
as.data.frame(df.env.prokary.vector_2016)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_prokary_2016 <- bind_cols(df.spe.prokary.vector_2016, df.spa.prokary.vector_2016, df.env.prokary.vector_2016) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_prokary_2016)

### Transforming into similarities = 1 - Bray
all_data_prokary_2016$similarity.prokary <- 1 - all_data_prokary_2016$dissimilarity.prokary
head(all_data_prokary_2016)


###############################################
############ 18S microeukaryotes GENUS ########
###############################################
### Read table metadata and abundances
microbes_18S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv", header=T)
attach(microbes_18S_genus)
names(microbes_18S_genus)

microbes_18S_genus_2016 <- filter(microbes_18S_genus, year == "2016") %>% 
  arrange(site)

samples_env_data_microeuk_2016 <- c("ZosCFIoldB16","ZosCFIoldG16","ZosCSSoldC16","ZosCSSoldF16","ZosCSSoldE16","ZosCSSoldA16","ZosCSSoldG16","ZosGSWoldF16","ZosGSWoldD16","ZosGSWoldA16","ZosGSWoldB16","ZosGSWoldG16","ZosGSWoldC16","ZosGSWoldE16","ZosGSWoldH16","ZosPBSoldA16","ZosPBSoldB16","ZosPBSoldC16","ZosPBSoldD16","ZosPBSoldG16","ZosPBSoldF16","ZosPBSoldH16","ZosPBSoldE16","ZosPBPoldH16","ZosPBPoldB16","ZosPBPoldG16","ZosPBPoldE16","ZosPBPoldC16","ZosPBPoldD16","ZosTQNoldG16","ZosTQNoldH16","ZosTQNoldF16","ZosTQNoldD16","ZosTQNoldB16","ZosTQNoldE16","ZosTQNoldC16","ZosTQNoldA16","ZosTQSoldG16","ZosTQSoldA16","ZosTQSoldH16","ZosTQSoldB16","ZosTQSoldF16","ZosTQSoldC16","ZosTQSoldD16")

microbes_18S_genus_2016 <- microbes_18S_genus_2016 %>% 
  dplyr::filter(SampleID %in% samples_env_data_microeuk_2016 )%>% 
  arrange(site)

microbes_18S_genus_2016_abund <- as.data.frame(microbes_18S_genus_2016 %>%
                                                 dplyr::select(-(1:9)) %>% 
                                                 rowSums())

# samples ZosPBPoldH16 and ZosPBPoldG16 have zero sum abundances, need to remove them
remove_zero_sum_samples <- c("ZosPBPoldH16", "ZosPBPoldG16")
microbes_18S_genus_2016 <- microbes_18S_genus_2016 %>% 
  filter(!SampleID %in% remove_zero_sum_samples)

#create a unique site_quadrat_id_year column
microbes_18S_genus_2016 <- microbes_18S_genus_2016 %>%
  unite(quadrat_year, site_quadrat_id, year, sep = "_" , remove = FALSE) %>%  
  select(SampleID, quadrat_year, everything())

microeuk_2016_samples <- as.data.frame(microbes_18S_genus_2016$SampleID)
write.csv(microeuk_2016_samples, "Data/micro_eukaryotes/distance_decay_microeuk_samples_2016_REMOVE_ENV.csv", row.names = F)

### Add environmental variables
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

microeuk_quad_envir <- inner_join(environmental_hakai, microbes_18S_genus_2016, by = "quadrat_year")

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
microeuk_dist_2016 <-read.csv("Data/micro_eukaryotes/microeuk_geogr_dist_matrix_2016_REMOVE_ENV.csv", header=TRUE) 
nrow(microeuk_dist_2016)
ncol(microeuk_dist_2016)
labels <- as.data.frame(microeuk_dist_2016$SampleID)
names(labels)[1] <- "SampleID"

microeuk_dist2_2016 <- dplyr::select(microeuk_dist_2016, c(-(SampleID)))
nrow(microeuk_dist2_2016)
ncol(microeuk_dist2_2016)
microeuk_matrix_2016 <- as.dist(microeuk_dist2_2016, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
final_microeuk_right_order <- left_join(labels, microeuk_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(final_microeuk_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
###get only abundance
names(final_microeuk_right_order)[1:40]
microeuk_abund_2016 <- final_microeuk_right_order %>% 
  dplyr::select(-(1:31)) 
names(microeuk_abund_2016)
nrow(microeuk_abund_2016)

### CREATE DISTANCE MATRIX FOR SPECIES
spe.dist.microeuk_2016 <-vegdist(microeuk_abund_2016,method="bray")
spe.dist.microeuk_2016
spe.dist.microeuk_2016.matrix <- as.matrix(spe.dist.microeuk_2016)
spe.dist.matrix_test_df <- as.data.frame(spe.dist.microeuk_2016.matrix)

### get environmental data
names(final_microeuk_right_order)[1:40]
microeuk_env_2016 <- final_microeuk_right_order %>% dplyr::select(c("depth", "temperature", "salinity", "quadrat_lai", "quadrat_microepiphyte_mg", "quadrat_biomass_g", "quadrat_macroalgae_g", "bed_area_m2"))
microeuk_env_2016 <- microeuk_env_2016

### CREATE DISTANCE MATRIX FOR ENVIRONMENT
env.dist.microeuk_2016 <-vegdist(log(microeuk_env_2016+1),method="euclidean")
env.dist.microeuk_2016
env.dist.microeuk_2016.matrix <- as.matrix(env.dist.microeuk_2016)
env.dist.matrix_test_df <- as.data.frame(env.dist.microeuk_2016.matrix)

df.spe.microeuk_2016 <- data.frame(dissimilarity=spe.dist.microeuk_2016[upper.tri(spe.dist.microeuk_2016, diag = FALSE)])
df.spe.microeuk_2016 #get row number that contain values to remove NAs
df.spe.microeuk.vector_2016 <- df.spe.microeuk_2016[1:210,]

df.env.microeuk_2016 <- data.frame(dissimilarity=env.dist.microeuk_2016[upper.tri(env.dist.microeuk_2016, diag = FALSE)])
df.env.microeuk_2016 #get row number that contain values to remove NAs
df.env.microeuk.vector_2016 <- df.env.microeuk_2016[1:210,]

###ggplot only work with data frames, so we need to convert this data into data frame
library(reshape2)
df.spe.microeuk.vector_2016 <- melt(as.matrix(spe.dist.microeuk_2016), varnames = c("row", "col"))
names(df.spe.microeuk.vector_2016)[names(df.spe.microeuk.vector_2016) == 'value'] <- 'dissimilarity.microeuk'
df.spe.microeuk.vector_2016
dissimilarity_new_microeuk_2016 <- df.spe.microeuk.vector_2016$dissimilarity.microeuk

df.spa.microeuk.vector_2016 <- melt(as.matrix(microeuk_matrix_2016), varnames = c("row", "col"))
names(df.spa.microeuk.vector_2016)[names(df.spa.microeuk.vector_2016) == 'value'] <- 'spatial_distance.microeuk'
df.spa.microeuk.vector_2016
spatial_microeuk_2016 <- df.spa.microeuk.vector_2016$spatial_distance.microeuk

df.env.microeuk.vector_2016 <- melt(as.matrix(env.dist.microeuk_2016), varnames = c("row", "col"))
names(df.env.microeuk.vector_2016)[names(df.env.microeuk.vector_2016) == 'value'] <- 'environmental_distance.microeuk'
df.env.microeuk.vector_2016
environment_microeuk_2016 <- df.env.microeuk.vector_2016$env_distance.microeuk

as.data.frame(df.spe.microeuk.vector_2016)
as.data.frame(df.spa.microeuk.vector_2016)
as.data.frame(df.env.microeuk.vector_2016)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_microeuk_2016 <- bind_cols(df.spe.microeuk.vector_2016, df.spa.microeuk.vector_2016, df.env.microeuk.vector_2016) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_microeuk_2016)

### Transforming into similarities = 1 - Bray
all_data_microeuk_2016$similarity.microeuk <- 1 - all_data_microeuk_2016$dissimilarity.microeuk
head(all_data_microeuk_2016)

#################################################
############ Inverts Macroeukaryotes ############
#################################################

inverts_finest <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")

inverts_finest_2016 <- filter(inverts_finest, year == "2016") %>% 
  arrange(site) 

inverts_finest_2016 <- inverts_finest_2016 %>% 
  dplyr::rename("SampleID" = "ID_year")

# *choked_inner_1_2016 = negative macroalgae
remove_no_env_inverts <- c("choked_inner_1_2016", "triquet_south_5_2016", "triquet_south_4_2017", "triquet_south_5_2017")
inverts_finest_2016 <- inverts_finest_2016 %>% 
  dplyr::filter(!SampleID %in% remove_no_env_inverts )

#create a unique site_quadrat_id_year column
inverts_finest_2016 <- inverts_finest_2016 %>%
  unite(quadrat_year, ID, year, sep = "_" , remove = FALSE) %>%  
  select(SampleID, quadrat_year, everything())

Sample_ID_inverts <- as.data.frame(inverts_finest_2016$SampleID)
write.csv(Sample_ID_inverts, "Data/macro_eukaryotes/distance_decay_macroeuk_samples_2016_REMOVE_ENV.csv", row.names = F)

### Add environmental variables
environmental_hakai <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv", header=T)
names(environmental_hakai)

inverts_quad_envir <- inner_join(environmental_hakai, inverts_finest_2016, by = "quadrat_year")

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
inverts_dist_2016 <-read.csv("Data/macro_eukaryotes/macroeuk_geogr_dist_matrix_2016_REMOVE_ENV.csv", header=TRUE) # had to put NA in all blank values

nrow(inverts_dist_2016)
ncol(inverts_dist_2016)
labels <- as.data.frame(inverts_dist_2016$SampleID)
names(labels)[1] <- "SampleID"

inverts_dist2_2016 <- dplyr::select(inverts_dist_2016, c(-(SampleID)))
nrow(inverts_dist2_2016)
ncol(inverts_dist2_2016)
inverts_matrix_2016 <- as.dist(inverts_dist2_2016, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
final_inverts_right_order <- left_join(labels, inverts_quad_envir, by="SampleID")
# check if it worked
label_table <- as.data.frame(final_inverts_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
###get only abundance
names(final_inverts_right_order)[1:40]
inverts_abund_2016 <- final_inverts_right_order %>% dplyr::select(-(1:29)) 
nrow(inverts_abund_2016)

### CREATE DISTANCE MATRIX FOR SPECIES
spe.dist.inverts_2016 <-vegdist(inverts_abund_2016,method="bray")
spe.dist.inverts_2016
spe.dist.inverts_2016.matrix <- as.matrix(spe.dist.inverts_2016)
spe.dist.matrix_test_df <- as.data.frame(spe.dist.inverts_2016.matrix)

### get environmental data
names(final_inverts_right_order)[1:40]
inverts_env_2016 <- final_inverts_right_order %>% dplyr::select(c("depth", "temperature", "salinity", "quadrat_lai", "quadrat_microepiphyte_mg", "quadrat_biomass_g", "quadrat_macroalgae_g", "bed_area_m2"))
inverts_env_2016 <- inverts_env_2016
d <- log(inverts_env_2016+1)
### CREATE DISTANCE MATRIX FOR ENVIRONMENT
env.dist.inverts_2016 <-vegdist(log(inverts_env_2016+1),method="euclidean")
env.dist.inverts_2016
env.dist.inverts_2016.matrix <- as.matrix(env.dist.inverts_2016)
env.dist.matrix_test_df <- as.data.frame(env.dist.inverts_2016.matrix)

df.spe.inverts_2016 <- data.frame(dissimilarity=spe.dist.inverts_2016[upper.tri(spe.dist.inverts_2016, diag = FALSE)])
df.spe.inverts_2016 #get row number that contain values to remove NAs
df.spe.inverts.vector_2016 <- df.spe.inverts_2016[1:120,]

df.env.inverts_2016 <- data.frame(dissimilarity=env.dist.inverts_2016[upper.tri(env.dist.inverts_2016, diag = FALSE)])
df.env.inverts_2016 #get row number that contain values to remove NAs
df.env.inverts.vector_2016 <- df.env.inverts_2016[1:120,]

###ggplot only work with data frames, so we need to convert this data into data frame
library(reshape2)
df.spe.inverts.vector_2016 <- melt(as.matrix(spe.dist.inverts_2016), varnames = c("row", "col"))
names(df.spe.inverts.vector_2016)[names(df.spe.inverts.vector_2016) == 'value'] <- 'dissimilarity.inverts'
df.spe.inverts.vector_2016
dissimilarity_new_inverts_2016 <- df.spe.inverts.vector_2016$dissimilarity.inverts

df.spa.inverts.vector_2016 <- melt(as.matrix(inverts_matrix_2016), varnames = c("row", "col"))
names(df.spa.inverts.vector_2016)[names(df.spa.inverts.vector_2016) == 'value'] <- 'spatial_distance.inverts'
df.spa.inverts.vector_2016
spatial_inverts_2016 <- df.spa.inverts.vector_2016$spatial_distance.inverts

df.env.inverts.vector_2016 <- melt(as.matrix(env.dist.inverts_2016), varnames = c("row", "col"))
names(df.env.inverts.vector_2016)[names(df.env.inverts.vector_2016) == 'value'] <- 'environmental_distance.inverts'
df.env.inverts.vector_2016
environmental_inverts_2016 <- df.env.inverts.vector_2016$environmental_distance.inverts

as.data.frame(df.spe.inverts.vector_2016)
as.data.frame(df.spa.inverts.vector_2016)
as.data.frame(df.env.inverts.vector_2016)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_inverts_2016 <- bind_cols(df.spe.inverts.vector_2016, df.spa.inverts.vector_2016, df.env.inverts.vector_2016) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_inverts_2016)

### Transforming into similarities = 1 - Bray
all_data_inverts_2016$similarity.inverts <- 1 - all_data_inverts_2016$dissimilarity.inverts
head(all_data_inverts_2016)

#####################################
######### Plot all together #########
#####################################

####### Prokaryotes ####### 
all_data_prokary_2016
nrow(all_data_prokary_2016)
head(all_data_prokary_2016)
# create labelling column for prokary
all_data_prokary_2016$host <- ("prokary_2016")
host_prokary_2016 <- as.data.frame(all_data_prokary_2016$host)
names(host_prokary_2016) <- "host"
head(host_prokary_2016)
prokary_ggplot_2016 <- all_data_prokary_2016 %>% dplyr::select(c(dissimilarity.prokary, spatial_distance.prokary, environmental_distance.prokary, host))
# rename(iris, Length = Sepal.Length)
prokary_ggplot2_2016 <- prokary_ggplot_2016 %>% dplyr::rename(dissimilarity = dissimilarity.prokary, spatial_distance = spatial_distance.prokary, environmental_distance = environmental_distance.prokary)
head(prokary_ggplot2_2016)

####### Micoeukaryotes ####### 
all_data_microeuk_2016
nrow(all_data_microeuk_2016)
head(all_data_microeuk_2016)
# create labelling column for microeuk
all_data_microeuk_2016$host <- ("microeuk_2016")
host_microeuk_2016 <- as.data.frame(all_data_microeuk_2016$host)
names(host_microeuk_2016) <- "host"
head(host_microeuk_2016)
microeuk_ggplot_2016 <- all_data_microeuk_2016 %>% dplyr::select(c(dissimilarity.microeuk, spatial_distance.microeuk, environmental_distance.microeuk, host))
# rename(iris, Length = Sepal.Length)
microeuk_ggplot2_2016 <- microeuk_ggplot_2016 %>% dplyr::rename(dissimilarity = dissimilarity.microeuk, spatial_distance = spatial_distance.microeuk, environmental_distance = environmental_distance.microeuk)
head(microeuk_ggplot2_2016)

####### Inverts ####### 
all_data_inverts_2016
nrow(all_data_inverts_2016)
head(all_data_inverts_2016)
# create labelling column for inverts
all_data_inverts_2016$host <- ("inverts_2016")
host_inverts_2016 <- as.data.frame(all_data_inverts_2016$host)
names(host_inverts_2016) <- "host"
head(host_inverts_2016)
inverts_ggplot_2016 <- all_data_inverts_2016 %>% dplyr::select(c(dissimilarity.inverts, spatial_distance.inverts, environmental_distance.inverts, host))
# rename(iris, Length = Sepal.Length)
inverts_ggplot2_2016 <- inverts_ggplot_2016 %>% dplyr::rename(dissimilarity = dissimilarity.inverts, spatial_distance = spatial_distance.inverts, environmental_distance = environmental_distance.inverts)
head(inverts_ggplot2_2016)

### combining all
all_ggplot <- bind_rows(prokary_ggplot2_2016, microeuk_ggplot2_2016, inverts_ggplot2_2016)
head(all_ggplot)
nrow(all_ggplot)
attach(all_ggplot)
class(host)
all_ggplot <- all_ggplot %>% mutate(host = factor(host, levels = c("prokary_2016", "microeuk_2016", "inverts_2016")))

### PLOT SIMILARITY X SPATIAL DISTANCE
all_ggplot$similarity <- 1 - all_ggplot$dissimilarity

mantel.prokary_geog <-mantel(spe.dist.prokary_2016, prokary_matrix_2016,permutations=9999, method = "pearson")
mantel.prokary_geog # Mantel statistic r: 0.3001 p<0.001

mantel.prokary_geog_remove_env <-mantel.partial(spe.dist.prokary_2016, prokary_matrix_2016, env.dist.prokary_2016,permutations=9999, method = "pearson")
mantel.prokary_geog_remove_env # Mantel statistic r: 0.2214 p=0.0011

mantel.microeuk_geog <-mantel(spe.dist.microeuk_2016, microeuk_matrix_2016,permutations=9999, method = "pearson")
mantel.microeuk_geog # Mantel statistic r: 0.29 p<0.001

mantel.microeuk_geog_remove_env <-mantel.partial(spe.dist.microeuk_2016, microeuk_matrix_2016, env.dist.microeuk_2016,permutations=9999, method = "pearson")
mantel.microeuk_geog_remove_env # Mantel statistic r: 0.2288 p=0.0014

mantel.inverts_geog <-mantel(spe.dist.inverts_2016, inverts_matrix_2016,permutations=9999, method = "pearson")
mantel.inverts_geog # Mantel statistic r: 0.4844 p<0.001

mantel.inverts_geog_remove_env <-mantel.partial(spe.dist.inverts_2016, inverts_matrix_2016, env.dist.inverts_2016,permutations=9999, method = "pearson")
mantel.inverts_geog_remove_env # Mantel statistic r: 0.4417 p<0.001

### DISTANCE IN KM 
# use mutate to create new column (spatial_distance_km) dividing spatial_distance by 1000 only if >0, else spatial_distance is kept
all_ggplot_km <- all_ggplot %>% mutate(spatial_distance_km = ifelse(spatial_distance > 0, spatial_distance/1000, spatial_distance))

# Fit regression line for similarity of prokary dist in km
all_data_prokary_2016_km <- all_data_prokary_2016 %>% mutate(spatial_distance.prokary_km = ifelse(spatial_distance.prokary > 0, spatial_distance.prokary/1000, spatial_distance.prokary))
require(stats)
reg_prokary_2016 <-lm(similarity.prokary ~ spatial_distance.prokary_km, data = all_data_prokary_2016_km)
reg_prokary_2016 #
summary(reg_prokary_2016)
# intercept: 0.4417084 +/- SE: 0.0051288
# slope: -0.0037059  +/- SE:0.0002211
reg_prokary_2016_REMOVE_ENV <-lm(similarity.prokary ~ spatial_distance.prokary_km + environmental_distance.prokary, data = all_data_prokary_2016_km)
reg_prokary_2016_REMOVE_ENV #
summary(reg_prokary_2016_REMOVE_ENV)
# intercept: 0.5707988 +/- SE: 0.0068034
# slope: -0.0022744  +/- SE:0.0002012

# Fit regression line for similarity of microeuk dist in km
all_data_microeuk_2016_km <- all_data_microeuk_2016 %>% mutate(spatial_distance.microeuk_km = ifelse(spatial_distance.microeuk > 0, spatial_distance.microeuk/1000, spatial_distance.microeuk))
require(stats)
reg_microeuk_2016 <-lm(similarity.microeuk ~ spatial_distance.microeuk_km, data = all_data_microeuk_2016_km)
reg_microeuk_2016 #
coeff_microeuk_2016 <- coefficients(reg_microeuk_2016) 
summary(reg_microeuk_2016)
# intercept: 0.4115357 +/- SE: 0.0107779
# slope: -0.0069718  +/- SE:0.0004546
reg_microeuk_2016_REMOVE_ENV <-lm(similarity.microeuk ~ spatial_distance.microeuk_km + environmental_distance.microeuk, data = all_data_microeuk_2016_km)
reg_microeuk_2016_REMOVE_ENV #
summary(reg_microeuk_2016_REMOVE_ENV)
# intercept: 0.6273193 +/- SE: 0.0148287
# slope: -0.0046735  +/- SE:0.0004296

# Fit regression line for similarity of inverts dist in km
all_data_inverts_2016_km <- all_data_inverts_2016 %>% mutate(spatial_distance.inverts_km = ifelse(spatial_distance.inverts > 0, spatial_distance.inverts/1000, spatial_distance.inverts))
require(stats)
reg_inverts_2016 <-lm(similarity.inverts ~ spatial_distance.inverts_km, data = all_data_inverts_2016_km)
reg_inverts_2016 #
coeff_inverts_2016 <- coefficients(reg_inverts_2016) 
summary(reg_inverts_2016)
# intercept: 0.5327516 +/- SE: 0.0088132
# slope: -0.0078764  +/- SE:0.0003923
reg_inverts_2016_REMOVE_ENV <-lm(similarity.inverts ~ spatial_distance.inverts_km + environmental_distance.inverts, data = all_data_inverts_2016_km)
reg_inverts_2016_REMOVE_ENV #
summary(reg_inverts_2016_REMOVE_ENV)
# intercept: 0.6441338 +/- SE: 0.0124449
# slope: -0.0064365  +/- SE:0.0003879

body_size_simi_km <- ggplot(all_ggplot_km, aes(x=spatial_distance_km, y=similarity)) +
  geom_point(size = 3,aes(color = host, alpha = host)) +
  scale_alpha_manual(values = c(0.8, 0.7, 0.6, 0.5), guide = FALSE) +
  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title + spacing top, right, bottom, left
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title = element_blank(), #remove title from legend
         legend.text = element_text(size=16), # size of legend text
         legend.text.align = 0,
         legend.position = c(0.82, 0.95), # set position to top right
         #legend.position = "right",
         legend.key.size = unit(0.8, "cm"), # increase distance between legend items
         legend.spacing = unit(5.0, 'cm')) # add spacing between legend text

body_size_simi_km <- body_size_simi_km + scale_colour_manual(values=c("#CA3542","#27647B", "#AEC0C9"))

body_size_simi_km <- body_size_simi_km + ggtitle("") +
  xlab("Spatial distance (km)") + ylab("Species composition similarity") +    
  # prokary
  geom_abline(intercept =  0.4417084, slope =  -0.0037059 , color="#CA3542",linetype="dashed", size = 1 ) +
  # prokary NO ENV
  geom_abline(intercept =  0.5707988, slope =  -0.0022744 , color="#CA3542",linetype="dotted", size = 1 ) +
  
  # microeuk
  geom_abline(intercept = 0.4115357, slope = -0.0069718, color="#27647B",linetype="dashed", size = 1 ) +
  # microeuk NO ENV
  geom_abline(intercept = 0.6273193, slope = -0.0046735, color="#27647B",linetype="dotted", size = 1 ) +
  
  # inverts
  geom_abline(intercept =  0.5327516, slope =   -0.0078764, color="#AEC0C9",linetype="dashed", size = 1) +
  # inverts NO ENV
  geom_abline(intercept =  0.6441338, slope =   -0.0064365, color="#AEC0C9",linetype="dotted", size = 1) 

body_size_simi_km
# ggsave("R_Code_and_Analysis/distance_decay/distance_decay_2016_REMOVE_ENV.tiff", plot = body_size_simi_km, width=310, height=250, units="mm",dpi=300, compression = "lzw", type = "cairo")

ggsave("R_Code_and_Analysis/distance_decay/distance_decay_2016_REMOVE_ENV.png", plot = body_size_simi_km, width=270, height=220, units="mm",dpi=300)


### COMPARE SLOPES 2016 ### 

# Fit the model, the covariate goes first
model_2016 <- lm(similarity ~ spatial_distance_km + environmental_distance * host, data = all_ggplot_km)
anova(model_2016)
summary(model_2016)
