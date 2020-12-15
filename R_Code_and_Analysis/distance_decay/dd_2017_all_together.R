### Distance Decay of Similarity 2017 ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 14th, 2020 ###

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

microbes_16S_genus_2017 <- filter(microbes_16S_genus, year == "2017") %>% 
  arrange(site)
prokary_2017_samples <- as.data.frame(microbes_16S_genus_2017$SampleID)
write.csv(prokary_2017_samples, "Data/prokaryotes/distance_decay_prokary_samples_2017.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
prokary_dist_2017 <-read.csv("Data/prokaryotes/prokary_geogr_dist_matrix_2017.csv", header=TRUE) 

labels <- as.data.frame(prokary_dist_2017$SampleID)
names(labels)[1] <- "SampleID"

prokary_dist2_2017 <- dplyr::select(prokary_dist_2017, c(-(SampleID)))
nrow(prokary_dist2_2017)

prokary_matrix_2017 <- as.dist(prokary_dist2_2017, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
final_prokary_right_order <- left_join(labels, microbes_16S_genus_2017, by="SampleID")
# check if it worked
label_table <- as.data.frame(final_prokary_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
###get only abundance
names(final_prokary_right_order)[1:20]
prokary_abund_2017 <- final_prokary_right_order %>% dplyr::select(-(1:12)) 
nrow(prokary_abund_2017)

### CREATE DISTANCE MATRIX FOR SPECIES
spe.dist.prokary_2017 <-vegdist(prokary_abund_2017,method="bray")
spe.dist.prokary_2017
spe.dist.prokary_2017.matrix <- as.matrix(spe.dist.prokary_2017)
spe.dist.matrix_test_df <- as.data.frame(spe.dist.prokary_2017.matrix)

df.spe.prokary_2017 <- data.frame(dissimilarity=spe.dist.prokary_2017[upper.tri(spe.dist.prokary_2017, diag = FALSE)])
df.spe.prokary_2017 #get row number that contain values to remove NAs
df.spe.prokary.vector_2017 <- df.spe.prokary_2017[1:91,]

###ggplot only work with data frames, so we need to convert this data into data frame
library(reshape2)
df.spe.prokary.vector_2017 <- melt(as.matrix(spe.dist.prokary_2017), varnames = c("row", "col"))
names(df.spe.prokary.vector_2017)[names(df.spe.prokary.vector_2017) == 'value'] <- 'dissimilarity.prokary'
df.spe.prokary.vector_2017
dissimilarity_new_prokary_2017 <- df.spe.prokary.vector_2017$dissimilarity.prokary

df.spa.prokary.vector_2017 <- melt(as.matrix(prokary_matrix_2017), varnames = c("row", "col"))
names(df.spa.prokary.vector_2017)[names(df.spa.prokary.vector_2017) == 'value'] <- 'spatial_distance.prokary'
df.spa.prokary.vector_2017
spatial_prokary_2017 <- df.spa.prokary.vector_2017$spatial_distance.prokary

as.data.frame(df.spe.prokary.vector_2017)
as.data.frame(df.spa.prokary.vector_2017)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_prokary_2017 <- bind_cols(df.spe.prokary.vector_2017, df.spa.prokary.vector_2017) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_prokary_2017)

### Transforming into similarities = 1 - Bray
all_data_prokary_2017$similarity.prokary <- 1 - all_data_prokary_2017$dissimilarity.prokary
head(all_data_prokary_2017)

###############################################
############ 18S microeukaryotes GENUS ########
###############################################
### Read table metadata and abundances
microbes_18S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv", header=T)
attach(microbes_18S_genus)
names(microbes_18S_genus)

microbes_18S_genus_2017 <- filter(microbes_18S_genus, year == "2017") %>%
  arrange(site)

microbes_18S_genus_2017_abund <- as.data.frame(microbes_18S_genus_2017 %>%
                                                 dplyr::select(-(1:9)) %>% 
                                                 rowSums())
# sample ZosCSSoldD17 have zero sum abundances, need to remove it
remove_zero_sum_samples <- c("ZosCSSoldD17")
microbes_18S_genus_2017 <- microbes_18S_genus_2017 %>% 
  filter(!SampleID %in% remove_zero_sum_samples)

microeuk_2017_samples <- as.data.frame(microbes_18S_genus_2017$SampleID)
write.csv(microeuk_2017_samples, "Data/micro_eukaryotes/distance_decay_microeuk_samples_2017.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
microeuk_dist_2017 <-read.csv("Data/micro_eukaryotes/microeuk_geogr_dist_matrix_2017.csv", header=TRUE) 
nrow(microeuk_dist_2017)
ncol(microeuk_dist_2017)
labels <- as.data.frame(microeuk_dist_2017$SampleID)
names(labels)[1] <- "SampleID"

microeuk_dist2_2017 <- dplyr::select(microeuk_dist_2017, c(-(SampleID)))
nrow(microeuk_dist2_2017)
ncol(microeuk_dist2_2017)
microeuk_matrix_2017 <- as.dist(microeuk_dist2_2017, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
final_microeuk_right_order <- left_join(labels, microbes_18S_genus_2017, by="SampleID")
# check if it worked
label_table <- as.data.frame(final_microeuk_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
###get only abundance
names(final_microeuk_right_order)[1:20]
microeuk_abund_2017 <- final_microeuk_right_order %>% dplyr::select(-(1:9))
nrow(microeuk_abund_2017)

### CREATE DISTANCE MATRIX FOR SPECIES
spe.dist.microeuk_2017 <-vegdist(microeuk_abund_2017,method="bray")
spe.dist.microeuk_2017
spe.dist.microeuk_2017.matrix <- as.matrix(spe.dist.microeuk_2017)
spe.dist.matrix_test_df <- as.data.frame(spe.dist.microeuk_2017.matrix)

df.spe.microeuk_2017 <- data.frame(dissimilarity=spe.dist.microeuk_2017[upper.tri(spe.dist.microeuk_2017, diag = FALSE)])
df.spe.microeuk_2017 #get row number that contain values to remove NAs
df.spe.microeuk.vector_2017 <- df.spe.microeuk_2017[1:6,]

###ggplot only work with data frames, so we need to convert this data into data frame
library(reshape2)
df.spe.microeuk.vector_2017 <- melt(as.matrix(spe.dist.microeuk_2017), varnames = c("row", "col"))
names(df.spe.microeuk.vector_2017)[names(df.spe.microeuk.vector_2017) == 'value'] <- 'dissimilarity.microeuk'
df.spe.microeuk.vector_2017
dissimilarity_new_microeuk_2017 <- df.spe.microeuk.vector_2017$dissimilarity.microeuk

df.spa.microeuk.vector_2017 <- melt(as.matrix(microeuk_matrix_2017), varnames = c("row", "col"))
names(df.spa.microeuk.vector_2017)[names(df.spa.microeuk.vector_2017) == 'value'] <- 'spatial_distance.microeuk'
df.spa.microeuk.vector_2017
spatial_microeuk_2017 <- df.spa.microeuk.vector_2017$spatial_distance.microeuk

as.data.frame(df.spe.microeuk.vector_2017)
as.data.frame(df.spa.microeuk.vector_2017)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_microeuk_2017 <- bind_cols(df.spe.microeuk.vector_2017, df.spa.microeuk.vector_2017) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_microeuk_2017)

### Transforming into similarities = 1 - Bray
all_data_microeuk_2017$similarity.microeuk <- 1 - all_data_microeuk_2017$dissimilarity.microeuk
head(all_data_microeuk_2017)

#################################################
############ Inverts Macroeukaryotes ############
#################################################

inverts_finest <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")

inverts_finest_2017 <- filter(inverts_finest, year == "2017") %>% 
  arrange(site) 

inverts_finest_2017 <- inverts_finest_2017 %>% 
  dplyr::rename("SampleID" = "ID_year")

inverts_finest_2017_abund <- as.data.frame(inverts_finest_2017 %>% 
  select(-c(1:7)) %>% 
  rowSums())

Sample_ID_inverts <- as.data.frame(inverts_finest_2017$SampleID)
write.csv(Sample_ID_inverts, "Data/macro_eukaryotes/distance_decay_macroeuk_samples_2017.csv", row.names = F)

### CREATE DISTANCE MATRIX FOR METERS CALCULATED IN GOOGLE EARTH
inverts_dist_2017 <-read.csv("Data/macro_eukaryotes/macroeuk_geogr_dist_matrix_2017.csv", header=TRUE) # had to put NA in all blank values
nrow(inverts_dist_2017)
ncol(inverts_dist_2017)
labels <- as.data.frame(inverts_dist_2017$SampleID)
names(labels)[1] <- "SampleID"

inverts_dist2_2017 <- dplyr::select(inverts_dist_2017, c(-(SampleID)))
nrow(inverts_dist2_2017)
ncol(inverts_dist2_2017)
inverts_matrix_2017 <- as.dist(inverts_dist2_2017, diag = TRUE)

########  IMPORTANT! TAXA TABLE HAVE TO BE IN THE SAME ORDER AS DISTANCE MATRIX!!! ###########
# putting FINAL TABLE in the same order as the distance matrix
final_inverts_right_order <- left_join(labels, inverts_finest_2017, by="SampleID")
# check if it worked
label_table <- as.data.frame(final_inverts_right_order$SampleID)
check <- cbind(label_table, labels)

# ### NOW CAN CREATE DIST MATRIX FOR ZOS
###get only abundance
names(final_inverts_right_order)[1:20]
inverts_abund_2017 <- final_inverts_right_order %>% dplyr::select(-(1:7)) 
nrow(inverts_abund_2017)

### CREATE DISTANCE MATRIX FOR SPECIES
spe.dist.inverts_2017 <-vegdist(inverts_abund_2017,method="bray")
spe.dist.inverts_2017
spe.dist.inverts_2017.matrix <- as.matrix(spe.dist.inverts_2017)
spe.dist.matrix_test_df <- as.data.frame(spe.dist.inverts_2017.matrix)

df.spe.inverts_2017 <- data.frame(dissimilarity=spe.dist.inverts_2017[upper.tri(spe.dist.inverts_2017, diag = FALSE)])
df.spe.inverts_2017 #get row number that contain values to remove NAs
df.spe.inverts.vector_2017 <- df.spe.inverts_2017[1:66,]

###ggplot only work with data frames, so we need to convert this data into data frame
library(reshape2)
df.spe.inverts.vector_2017 <- melt(as.matrix(spe.dist.inverts_2017), varnames = c("row", "col"))
names(df.spe.inverts.vector_2017)[names(df.spe.inverts.vector_2017) == 'value'] <- 'dissimilarity.inverts'
df.spe.inverts.vector_2017
dissimilarity_new_inverts_2017 <- df.spe.inverts.vector_2017$dissimilarity.inverts

df.spa.inverts.vector_2017 <- melt(as.matrix(inverts_matrix_2017), varnames = c("row", "col"))
names(df.spa.inverts.vector_2017)[names(df.spa.inverts.vector_2017) == 'value'] <- 'spatial_distance.inverts'
df.spa.inverts.vector_2017
spatial_inverts_2017 <- df.spa.inverts.vector_2017$spatial_distance.inverts

as.data.frame(df.spe.inverts.vector_2017)
as.data.frame(df.spa.inverts.vector_2017)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_inverts_2017 <- bind_cols(df.spe.inverts.vector_2017, df.spa.inverts.vector_2017) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_inverts_2017)

### Transforming into similarities = 1 - Bray
all_data_inverts_2017$similarity.inverts <- 1 - all_data_inverts_2017$dissimilarity.inverts
head(all_data_inverts_2017)

#####################################
######### Plot all together #########
#####################################

####### Prokaryotes ####### 
all_data_prokary_2017
nrow(all_data_prokary_2017)
head(all_data_prokary_2017)
# create labelling column for prokary
all_data_prokary_2017$host <- ("prokary_2017")
host_prokary_2017 <- as.data.frame(all_data_prokary_2017$host)
names(host_prokary_2017) <- "host"
head(host_prokary_2017)
prokary_ggplot_2017 <- all_data_prokary_2017 %>% dplyr::select(c(dissimilarity.prokary, spatial_distance.prokary, host))
# rename(iris, Length = Sepal.Length)
prokary_ggplot2_2017 <- prokary_ggplot_2017 %>% dplyr::rename(dissimilarity = dissimilarity.prokary, spatial_distance = spatial_distance.prokary)
head(prokary_ggplot2_2017)

####### Microeukaryotes ####### 
all_data_microeuk_2017
nrow(all_data_microeuk_2017)
head(all_data_microeuk_2017)
# create labelling column for microeuk
all_data_microeuk_2017$host <- ("microeuk_2017")
host_microeuk_2017 <- as.data.frame(all_data_microeuk_2017$host)
names(host_microeuk_2017) <- "host"
head(host_microeuk_2017)
microeuk_ggplot_2017 <- all_data_microeuk_2017 %>% dplyr::select(c(dissimilarity.microeuk, spatial_distance.microeuk, host))
# rename(iris, Length = Sepal.Length)
microeuk_ggplot2_2017 <- microeuk_ggplot_2017 %>% dplyr::rename(dissimilarity = dissimilarity.microeuk, spatial_distance = spatial_distance.microeuk)
head(microeuk_ggplot2_2017)

####### Inverts ####### 
all_data_inverts_2017
nrow(all_data_inverts_2017)
head(all_data_inverts_2017)
# create labelling column for inverts
all_data_inverts_2017$host <- ("inverts_2017")
host_inverts_2017 <- as.data.frame(all_data_inverts_2017$host)
names(host_inverts_2017) <- "host"
head(host_inverts_2017)
inverts_ggplot_2017 <- all_data_inverts_2017 %>% dplyr::select(c(dissimilarity.inverts, spatial_distance.inverts, host))
# rename(iris, Length = Sepal.Length)
inverts_ggplot2_2017 <- inverts_ggplot_2017 %>% dplyr::rename(dissimilarity = dissimilarity.inverts, spatial_distance = spatial_distance.inverts)
head(inverts_ggplot2_2017)

### combining all
all_ggplot <- bind_rows(prokary_ggplot2_2017, microeuk_ggplot2_2017, inverts_ggplot2_2017)
head(all_ggplot)
nrow(all_ggplot)
attach(all_ggplot)
class(host)
all_ggplot <- all_ggplot %>% mutate(host = factor(host, levels = c("prokary_2017", "microeuk_2017", "inverts_2017")))

### PLOT SIMILARITY X SPATIAL DISTANCE
all_ggplot$similarity <- 1 - all_ggplot$dissimilarity

# mantel.zostera_2015 <-mantel(spe.dist.zostera_2015, zos_matrix_2015,permutations=9999, method = "pearson") 
# mantel.zostera_2015 # Mantel statistic r: 0.4897  
# 
# mantel.zostera_2017 <-mantel(spe.dist.zostera_2017, zos_matrix_2017,permutations=9999, method = "pearson") 
# mantel.zostera_2017 # Mantel statistic r: 0.4728 
# 
# mantel.seawater_2015 <-mantel(spe.dist.seawater_2015, seawater_matrix_2015,permutations=9999, method = "pearson") 
# mantel.seawater_2015 # Mantel statistic r: 0.2914 
# 
# mantel.seawater_2017 <-mantel(spe.dist.seawater_2017, seawater_matrix_2017,permutations=9999, method = "pearson") 
# mantel.seawater_2017 # Mantel statistic r: 0.5645 

### DISTANCE IN KM 
# use mutate to create new column (spatial_distance_km) dividing spatial_distance by 1000 only if >0, else spatial_distance is kept
all_ggplot_km <- all_ggplot %>% mutate(spatial_distance_km = ifelse(spatial_distance > 0, spatial_distance/1000, spatial_distance))

# Fit regression line for similarity of prokary dist in km
all_data_prokary_2017_km <- all_data_prokary_2017 %>% mutate(spatial_distance.prokary_km = ifelse(spatial_distance.prokary > 0, spatial_distance.prokary/1000, spatial_distance.prokary))
require(stats)
reg_zos_2017 <-lm(similarity.prokary ~ spatial_distance.prokary_km, data = all_data_prokary_2017_km)
reg_zos_2017 #
coeff_zos_2017 <- coefficients(reg_zos_2017) 
summary(reg_zos_2017)
# intercept: 0.4485770 +/- SE: 0.0116323
# slope: -0.0088253  +/- SE:0.0007543

# Fit regression line for similarity of microeuk dist in km
all_data_microeuk_2017_km <- all_data_microeuk_2017 %>% mutate(spatial_distance.microeuk_km = ifelse(spatial_distance.microeuk > 0, spatial_distance.microeuk/1000, spatial_distance.microeuk))
require(stats)
reg_zos_2017 <-lm(similarity.microeuk ~ spatial_distance.microeuk_km, data = all_data_microeuk_2017_km)
reg_zos_2017 #
coeff_zos_2017 <- coefficients(reg_zos_2017) 
summary(reg_zos_2017)
# intercept: 0.504840 +/- SE: 0.043331
# slope: -0.022027  +/- SE:0.002805

# Fit regression line for similarity of inverts dist in km
all_data_inverts_2017_km <- all_data_inverts_2017 %>% mutate(spatial_distance.inverts_km = ifelse(spatial_distance.inverts > 0, spatial_distance.inverts/1000, spatial_distance.inverts))
require(stats)
reg_zos_2017 <-lm(similarity.inverts ~ spatial_distance.inverts_km, data = all_data_inverts_2017_km)
reg_zos_2017 #
coeff_zos_2017 <- coefficients(reg_zos_2017) 
summary(reg_zos_2017)
# intercept: 0.5784167 +/- SE: 0.0100680
# slope: -0.0115578  +/- SE:0.0006948

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
  geom_abline(intercept =  0.438849, slope =  -0.008338 , color="#CA3542",linetype="dashed", size = 1 ) +
  # microeuk
  geom_abline(intercept = 0.496475, slope = -0.021645, color="#27647B",linetype="dashed", size = 1 ) +
  # inverts
  geom_abline(intercept =  0.5882494, slope =   -0.0123179, color="#AEC0C9",linetype="dashed", size = 1) 

ggsave("R_Code_and_Analysis/distance_decay/distance_decay_2017.tiff", plot = body_size_simi_km, width=310, height=250, units="mm",dpi=300, compression = "lzw", type = "cairo")

