### Mantel test microbes and grazers ###
### Author: Bianca Trevizan Segovia ###
### Date created: January 13, 2020 ###

library(dplyr)
library(vegan)
library(usedist)

###################################
############### 2017 ##############
###################################

#load grazers metadata
metadata_grazers_2017 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/2017_grazer_metadata.csv")

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
metadata_grazers_2017$labels <- vegan::make.cepnames(metadata_grazers_2017$sample)

#load grazers bray curtis dissimilarity matrix
df_grazers_2017 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/2017_grazer_braycurtis.csv")

#load 16S microbial distance matrix GENUS
df_16S_2017_genus <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/genus_16S_2017_braycurtis.csv")
df_16S_2017_genus <- df_16S_2017_genus %>% 
  dplyr::rename("sample" = "X")

#load 16S microbial distance matrix FAMILY
df_16S_2017_family <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/family_16S_2017_braycurtis.csv")
df_16S_2017_family <- df_16S_2017_family %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix GENUS
df_18S_2017_genus <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/genus_18S_2017_braycurtis.csv")
df_18S_2017_genus <- df_18S_2017_genus %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix FAMILY
df_18S_2017_family <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/family_18S_2017_braycurtis.csv")
df_18S_2017_family <- df_18S_2017_family %>% 
  dplyr::rename("sample" = "X")

##############################################
######### MANTEL GRAZERS VS 16S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2017 <- as.vector(df_grazers_2017$sample)
labels_16S_genus_2017 <- as.vector(df_16S_2017_genus$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2017, labels_16S_genus_2017)

df_grazers_2017_mantel <- df_grazers_2017 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2017 <- dist(df_grazers_2017_mantel, diag = TRUE, upper = FALSE)

df_16S_genus_2017_mantel <- df_16S_2017_genus %>% 
  dplyr::filter(sample %in% common_samples)
dist_16S_genus_2017 <- dist(df_16S_genus_2017_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 16S genus in 2017

mantel_grazes_16S_genus <-  mantel(dist_grazers_2017, dist_16S_genus_2017, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_16S <- capture.output(mantel_grazes_16S_genus)
write.table(as.data.frame(out_mantel_grazes_16S ), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_16S_genus_2017.csv", quote=F, row.names=F, col.names=T, sep="n")

#############################################
######### MANTEL GRAZERS VS 18S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2017 <- as.vector(df_grazers_2017$sample)
labels_18S_genus_2017 <- as.vector(df_18S_2017_genus$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2017, labels_18S_genus_2017) ### ******* only 5 samples were common here

df_grazers_2017_mantel <- df_grazers_2017 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2017 <- dist(df_grazers_2017_mantel, diag = TRUE, upper = FALSE)

df_18S_genus_2017_mantel <- df_18S_2017_genus %>% 
  dplyr::filter(sample %in% common_samples)
dist_18S_genus_2017 <- dist(df_18S_genus_2017_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 18S genus in 2017

mantel_grazes_18S_genus <-  mantel(dist_grazers_2017, dist_18S_genus_2017, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_18S <- capture.output(mantel_grazes_18S_genus)
write.table(as.data.frame(out_mantel_grazes_18S), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_18S_genus_2017.csv", quote=F, row.names=F, col.names=T, sep="n")

##############################################
######### MANTEL GRAZERS VS 16S FAMILY ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2017 <- as.vector(df_grazers_2017$sample)
labels_16S_family_2017 <- as.vector(df_16S_2017_family$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2017, labels_16S_family_2017)

df_grazers_2017_mantel <- df_grazers_2017 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2017 <- dist(df_grazers_2017_mantel, diag = TRUE, upper = FALSE)

df_16S_family_2017_mantel <- df_16S_2017_family %>% 
  dplyr::filter(sample %in% common_samples)
dist_16S_family_2017 <- dist(df_16S_family_2017_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 16S family in 2017

mantel_grazes_16S_family <-  mantel(dist_grazers_2017, dist_16S_family_2017, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_16S <- capture.output(mantel_grazes_16S_family)
write.table(as.data.frame(out_mantel_grazes_16S ), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_16S_family_2017.csv", quote=F, row.names=F, col.names=T, sep="n")

#############################################
######### MANTEL GRAZERS VS 18S FAMILY ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2017 <- as.vector(df_grazers_2017$sample)
labels_18S_family_2017 <- as.vector(df_18S_2017_family$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2017, labels_18S_family_2017) ### ******* only 5 samples were common here

df_grazers_2017_mantel <- df_grazers_2017 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2017 <- dist(df_grazers_2017_mantel, diag = TRUE, upper = FALSE)

df_18S_family_2017_mantel <- df_18S_2017_family %>% 
  dplyr::filter(sample %in% common_samples)
dist_18S_family_2017 <- dist(df_18S_family_2017_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 18S family in 2017

mantel_grazes_18S_family <-  mantel(dist_grazers_2017, dist_18S_family_2017, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_18S <- capture.output(mantel_grazes_18S_family)
write.table(as.data.frame(out_mantel_grazes_18S), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_18S_family_2017.csv", quote=F, row.names=F, col.names=T, sep="n")
