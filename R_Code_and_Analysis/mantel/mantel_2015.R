### Mantel test microbes and grazers ###
### Author: Bianca Trevizan Segovia ###
### Date created: February 24, 2020 ###

library(dplyr)
library(vegan)
library(usedist)

###################################
############### 2015 ##############
###################################

#load grazers metadata
metadata_grazers_2015 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/2015_grazer_metadata.csv")

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
metadata_grazers_2015$labels <- vegan::make.cepnames(metadata_grazers_2015$sample)

#load grazers bray curtis dissimilarity matrix
df_grazers_2015 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/2015_grazer_braycurtis.csv")

#load 16S microbial distance matrix GENUS
df_16S_2015_genus <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/genus_16S_2015_braycurtis.csv")
df_16S_2015_genus <- df_16S_2015_genus %>% 
  dplyr::rename("sample" = "X")

#load 16S microbial distance matrix FAMILY
df_16S_2015_family <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/family_16S_2015_braycurtis.csv")
df_16S_2015_family <- df_16S_2015_family %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix GENUS
df_18S_2015_genus <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/genus_18S_2015_braycurtis.csv")
df_18S_2015_genus <- df_18S_2015_genus %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix FAMILY
df_18S_2015_family <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/family_18S_2015_braycurtis.csv")
df_18S_2015_family <- df_18S_2015_family %>% 
  dplyr::rename("sample" = "X")

##############################################
######### MANTEL GRAZERS VS 16S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2015 <- as.vector(df_grazers_2015$sample)
labels_16S_genus_2015 <- as.vector(df_16S_2015_genus$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2015, labels_16S_genus_2015)

df_grazers_2015_mantel <- df_grazers_2015 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2015 <- dist(df_grazers_2015_mantel, diag = TRUE, upper = FALSE)

df_16S_genus_2015_mantel <- df_16S_2015_genus %>% 
  dplyr::filter(sample %in% common_samples)
dist_16S_genus_2015 <- dist(df_16S_genus_2015_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 16S genus in 2015

mantel_grazes_16S_genus <-  mantel(dist_grazers_2015, dist_16S_genus_2015, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_16S <- capture.output(mantel_grazes_16S_genus)
write.table(as.data.frame(out_mantel_grazes_16S ), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_16S_genus_2015.csv", quote=F, row.names=F, col.names=T, sep="n")

#############################################
######### MANTEL GRAZERS VS 18S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2015 <- as.vector(df_grazers_2015$sample)
labels_18S_genus_2015 <- as.vector(df_18S_2015_genus$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2015, labels_18S_genus_2015) ### ******* only 5 samples were common here

df_grazers_2015_mantel <- df_grazers_2015 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2015 <- dist(df_grazers_2015_mantel, diag = TRUE, upper = FALSE)

df_18S_genus_2015_mantel <- df_18S_2015_genus %>% 
  dplyr::filter(sample %in% common_samples)
dist_18S_genus_2015 <- dist(df_18S_genus_2015_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 18S genus in 2015

mantel_grazes_18S_genus <-  mantel(dist_grazers_2015, dist_18S_genus_2015, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_18S <- capture.output(mantel_grazes_18S_genus)
write.table(as.data.frame(out_mantel_grazes_18S), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_18S_genus_2015.csv", quote=F, row.names=F, col.names=T, sep="n")

##############################################
######### MANTEL GRAZERS VS 16S FAMILY ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2015 <- as.vector(df_grazers_2015$sample)
labels_16S_family_2015 <- as.vector(df_16S_2015_family$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2015, labels_16S_family_2015)

df_grazers_2015_mantel <- df_grazers_2015 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2015 <- dist(df_grazers_2015_mantel, diag = TRUE, upper = FALSE)

df_16S_family_2015_mantel <- df_16S_2015_family %>% 
  dplyr::filter(sample %in% common_samples)
dist_16S_family_2015 <- dist(df_16S_family_2015_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 16S family in 2015

mantel_grazes_16S_family <-  mantel(dist_grazers_2015, dist_16S_family_2015, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_16S <- capture.output(mantel_grazes_16S_family)
write.table(as.data.frame(out_mantel_grazes_16S ), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_16S_family_2015.csv", quote=F, row.names=F, col.names=T, sep="n")

#############################################
######### MANTEL GRAZERS VS 18S FAMILY ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2015 <- as.vector(df_grazers_2015$sample)
labels_18S_family_2015 <- as.vector(df_18S_2015_family$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2015, labels_18S_family_2015) ### ******* only 5 samples were common here

df_grazers_2015_mantel <- df_grazers_2015 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2015 <- dist(df_grazers_2015_mantel, diag = TRUE, upper = FALSE)

df_18S_family_2015_mantel <- df_18S_2015_family %>% 
  dplyr::filter(sample %in% common_samples)
dist_18S_family_2015 <- dist(df_18S_family_2015_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 18S family in 2015

mantel_grazes_18S_family <-  mantel(dist_grazers_2015, dist_18S_family_2015, method = "spearman", permutations = 9999, na.rm = TRUE)
out_mantel_grazes_18S <- capture.output(mantel_grazes_18S_family)
write.table(as.data.frame(out_mantel_grazes_18S), file = "~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/out_mantel_grazes_18S_family_2015.csv", quote=F, row.names=F, col.names=T, sep="n")