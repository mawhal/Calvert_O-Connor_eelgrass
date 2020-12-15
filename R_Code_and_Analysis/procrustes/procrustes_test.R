### PROCRUSTES microbes and grazers ###
### Author: Bianca Trevizan Segovia ###
### Date created: March 02, 2020 ###

library(tidyverse)
library(vegan)
#library(usedist)

# pick a year
year <- 2018

#load grazers metadata
metadata_macro <- read_csv(paste0("R_Code_and_Analysis/mantel/",year,"_macroeuk_metadata.csv") )

#The row and column names in the distance matrix can be created from metadata$sample using vegan::make.cepnames()
metadata_macro$labels <- vegan::make.cepnames(metadata_macro$sample)

#load grazers bray curtis dissimilarity matrix
df_macro <- read_csv(paste0("R_Code_and_Analysis/mantel/",year,"_macroeuk_braycurtis.csv") )

df_macro$sample <- colnames(df_macro)

#load 16S microbial distance matrix GENUS
df_16S_meta  <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_metadata.csv") )
df_16S_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_16S_",year,"_braycurtis.csv") )
df_16S_genus <- df_16S_genus %>% 
  # dplyr::rename("sample" = "X1")
  mutate( sample = df_16S_meta$labels )

#load 18S microbial distance matrix GENUS
df_18S_meta  <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_metadata.csv") )
df_18S_genus <- read_csv(paste0("R_Code_and_Analysis/mantel/genus_18S_",year,"_braycurtis.csv") )
df_18S_genus <- df_18S_genus %>% 
  # dplyr::rename("sample" = "X1")
  mutate( sample = df_18S_meta$labels )


##############################################
######### PROCRUSTES MACROEukaryotes VS 16S GENUS ######### 
##############################################

# which samples are common to both matrices?
labels_macro <- as.vector(df_macro$sample)
# labels_macro <- as.vector(metadata_macro$labels)
labels_16S_genus <- as.vector(df_16S_genus$sample)
labels_18S_genus <- as.vector(df_18S_genus$sample)

# filter only samples that are common to both 
common_samples1 <- intersect(labels_macro, labels_16S_genus)

df_macro_mantel <- df_macro %>% 
  dplyr::filter(sample %in% common_samples1)
dist_macro <- dist(df_macro_mantel, diag = TRUE, upper = FALSE)

df_16S_genus_mantel <- df_16S_genus %>% 
  dplyr::filter(sample %in% common_samples1)
dist_16S_genus <- dist(df_16S_genus_mantel , diag = TRUE, upper = FALSE)


mds.target <- monoMDS(dist_macro)
mds.rotated <- monoMDS(dist_16S_genus)
proc_macro_16S <- procrustes(mds.target, mds.rotated)
proc_macro_16S 
summary(proc_macro_16S)
plot(proc_macro_16S)
plot(proc_macro_16S, kind=2)
residuals(proc_macro_16S)
p1 <- protest(mds.target, mds.rotated)
p1$signif
p1$ss

##############################################
######### PROCRUSTES MACROEukaryotes VS 18S GENUS ######### 
##############################################

# which samples are common to both matrices?

# filter only samples that are common to both 
common_samples2 <- intersect(labels_macro, labels_18S_genus)

df_macro_mantel <- df_macro %>% 
  dplyr::filter(sample %in% common_samples2)
dist_macro <- dist(df_macro_mantel, diag = TRUE, upper = FALSE)

df_18S_genus_mantel <- df_18S_genus %>% 
  dplyr::filter(sample %in% common_samples2)
dist_18S_genus <- dist(df_18S_genus_mantel , diag = TRUE, upper = FALSE)


mds.target <- monoMDS(dist_macro)
mds.rotated <- monoMDS(dist_18S_genus)
proc_macro_16S <- procrustes(mds.target, mds.rotated)
proc_macro_16S 
summary(proc_macro_16S)
plot(proc_macro_16S)
plot(proc_macro_16S, kind=2)
residuals(proc_macro_16S)
p2 <- protest(mds.target, mds.rotated)
p2$signif
p2$ss

##############################################
######### PROCRUSTES 18S GENUS VS 16S GENUS ######### 
##############################################

# which samples are common to both matrices?

# filter only samples that are common to both 
common_samples3 <- intersect(labels_18S_genus, labels_16S_genus)

df_16S_genus_mantel <- df_16S_genus %>% 
  dplyr::filter(sample %in% common_samples3)
dist_16S_genus <- dist(df_16S_genus_mantel , diag = TRUE, upper = FALSE)
df_18S_genus_mantel <- df_18S_genus %>% 
  dplyr::filter(sample %in% common_samples3)
dist_18S_genus <- dist(df_18S_genus_mantel , diag = TRUE, upper = FALSE)


mds.target <- monoMDS(dist_18S_genus)
mds.rotated <- monoMDS(dist_16S_genus)
proc_macro_16S <- procrustes(mds.target, mds.rotated)
proc_macro_16S 
summary(proc_macro_16S)
plot(proc_macro_16S)
plot(proc_macro_16S, kind=2)
residuals(proc_macro_16S)
p3 <- protest(mds.target, mds.rotated, permutations = how(nperm = 9999) )
p3$signif
p3$ss
# sqrt( 1-p3$ss )


# wrap results together
res <- data.frame(year=year, pair=c("macro-16S","macro-18S","18S-16S"),
                  samples=c(length(common_samples1),length(common_samples2),length(common_samples3)),
                  matrix( c(p1$ss,p1$signif,p2$ss,p2$signif,p3$ss,p3$signif), ncol=2, byrow=T ))
names(res) <- c('year','pair','samples','ss',"signif")

res$corr <- sqrt(1-res$ss)
res

# write to disk
write_csv( res, paste0("R_Code_and_Analysis/procrustes/procrustes_",year,".csv") )

           