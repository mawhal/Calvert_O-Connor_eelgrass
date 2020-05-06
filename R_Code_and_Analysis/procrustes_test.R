### PROCRUSTES microbes and grazers ###
### Author: Bianca Trevizan Segovia ###
### Date created: March 02, 2020 ###

library(tidyverse)
library(vegan)
#library(usedist)

###################################
############### 2016 ##############
###################################

#load grazers metadata
metadata_grazers_2016 <- read_csv("R_Code_and_Analysis/output_data/2016_grazer_metadata.csv")

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
metadata_grazers_2016$labels <- vegan::make.cepnames(metadata_grazers_2016$sample)

#load grazers bray curtis dissimilarity matrix
df_grazers_2016 <- read.csv("R_Code_and_Analysis/mantel_microbes_grazers/2016_grazer_braycurtis.csv")

df_grazers_2016$sample <- colnames(df_grazers_2016)

#load 16S microbial distance matrix GENUS
df_16S_2016_genus <- read.csv("R_Code_and_Analysis/mantel_microbes_grazers/genus_16S_2016_braycurtis.csv")
df_16S_2016_genus <- df_16S_2016_genus %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix GENUS
df_18S_2016_genus <- read.csv("R_Code_and_Analysis/mantel_microbes_grazers/genus_18S_2016_braycurtis.csv")
df_18S_2016_genus <- df_18S_2016_genus %>% 
  dplyr::rename("sample" = "X")


##############################################
######### PROCRUSTES GRAZERS VS 16S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2016 <- as.vector(df_grazers_2016$sample)
labels_16S_genus_2016 <- as.vector(df_16S_2016_genus$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2016, labels_16S_genus_2016)

df_grazers_2016_mantel <- df_grazers_2016 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2016 <- dist(df_grazers_2016_mantel, diag = TRUE, upper = FALSE)

df_16S_genus_2016_mantel <- df_16S_2016_genus %>% 
  dplyr::filter(sample %in% common_samples)
dist_16S_genus_2016 <- dist(df_16S_genus_2016_mantel , diag = TRUE, upper = FALSE)


mds.invert <- monoMDS(dist_grazers_2016)
mds.16S <- monoMDS(dist_16S_genus_2016)
proc_invert_16S <- procrustes(mds.invert, mds.16S)
proc_invert_16S 
summary(proc_invert_16S)
plot(proc_invert_16S)
plot(proc_invert_16S, kind=2)
residuals(proc_invert_16S)
protest(mds.invert, mds.16S)
