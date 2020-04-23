### Mantel test microbes and grazers ###
### Author: Bianca Trevizan Segovia ###
### Date created: January 13, 2020 ###

library(dplyr)
library(vegan)
library(usedist)

###################################
############### 2016 ##############
###################################

#load grazers metadata
metadata_grazers_2016 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/2016_grazer_metadata.csv")

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
metadata_grazers_2016$labels <- vegan::make.cepnames(metadata_grazers_2016$sample)

#load grazers bray curtis dissimilarity matrix
df_grazers_2016 <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/2016_grazer_braycurtis.csv")

df_grazers_2016$sample <- colnames(df_grazers_2016)

### move labels from last column to first
df_grazers_2016 <- df_grazers_2016 %>%
  dplyr::select(sample, everything())

#load 16S microbial distance matrix GENUS
df_16S_2016_genus <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/genus_16S_2016_braycurtis.csv")
df_16S_2016_genus <- df_16S_2016_genus %>% 
  dplyr::rename("sample" = "X")

#load 16S microbial distance matrix FAMILY
df_16S_2016_family <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/family_16S_2016_braycurtis.csv")
df_16S_2016_family <- df_16S_2016_family %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix GENUS
df_18S_2016_genus <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/genus_18S_2016_braycurtis.csv")
df_18S_2016_genus <- df_18S_2016_genus %>% 
  dplyr::rename("sample" = "X")

#load 18S microbial distance matrix FAMILY
df_18S_2016_family <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/mantel_microbes_grazers/family_18S_2016_braycurtis.csv")
df_18S_2016_family <- df_18S_2016_family %>% 
  dplyr::rename("sample" = "X")

##############################################
######### MANTEL GRAZERS VS 16S GENUS ######### 
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

# MANTEL test on grazers vs 16S genus in 2016

mantel_grazers_16S_genus <-  mantel(dist_grazers_2016, dist_16S_genus_2016, method = "spearman", permutations = 9999, na.rm = TRUE)

#############################################
######### MANTEL GRAZERS VS 18S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2016 <- as.vector(df_grazers_2016$sample)
labels_18S_genus_2016 <- as.vector(df_18S_2016_genus$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2016, labels_18S_genus_2016) ### ******* only 5 samples were common here

df_grazers_2016_mantel <- df_grazers_2016 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2016 <- dist(df_grazers_2016_mantel, diag = TRUE, upper = FALSE)

df_18S_genus_2016_mantel <- df_18S_2016_genus %>% 
  dplyr::filter(sample %in% common_samples)
dist_18S_genus_2016 <- dist(df_18S_genus_2016_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 18S genus in 2016

mantel_grazers_18S_genus <-  mantel(dist_grazers_2016, dist_18S_genus_2016, method = "spearman", permutations = 9999, na.rm = TRUE)

##############################################
######### MANTEL GRAZERS VS 16S FAMILY ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2016 <- as.vector(df_grazers_2016$sample)
labels_16S_family_2016 <- as.vector(df_16S_2016_family$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2016, labels_16S_family_2016)

df_grazers_2016_mantel <- df_grazers_2016 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2016 <- dist(df_grazers_2016_mantel, diag = TRUE, upper = FALSE)

df_16S_family_2016_mantel <- df_16S_2016_family %>% 
  dplyr::filter(sample %in% common_samples)
dist_16S_family_2016 <- dist(df_16S_family_2016_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 16S family in 2016

mantel_grazers_16S_family <-  mantel(dist_grazers_2016, dist_16S_family_2016, method = "spearman", permutations = 9999, na.rm = TRUE)

#############################################
######### MANTEL GRAZERS VS 18S FAMILY ######### 
##############################################

# which samples are common to both matrices?
labels_grazers_2016 <- as.vector(df_grazers_2016$sample)
labels_18S_family_2016 <- as.vector(df_18S_2016_family$sample)

# filter only samples that are common to both 
common_samples <- intersect(labels_grazers_2016, labels_18S_family_2016) ### ******* only 5 samples were common here

df_grazers_2016_mantel <- df_grazers_2016 %>% 
  dplyr::filter(sample %in% common_samples)
dist_grazers_2016 <- dist(df_grazers_2016_mantel, diag = TRUE, upper = FALSE)

df_18S_family_2016_mantel <- df_18S_2016_family %>% 
  dplyr::filter(sample %in% common_samples)
dist_18S_family_2016 <- dist(df_18S_family_2016_mantel , diag = TRUE, upper = FALSE)

# MANTEL test on grazers vs 18S family in 2016

mantel_grazers_18S_family <-  mantel(dist_grazers_2016, dist_18S_family_2016, method = "spearman", permutations = 9999, na.rm = TRUE)


#### WRITE OUTPUTS ALTOGETHER ###
mantel_extract <- function(x) {
  data.frame(stat = x$statistic, p= x$signif)
}

mantel_extract(mantel_grazers_16S_genus)
