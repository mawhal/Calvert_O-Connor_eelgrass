### Mantel test microbes and grazers ###
### Author: Bianca Trevizan Segovia ###
### Date created: January 13, 2020 ###
### modified by Matt Whalen #---> ###
### (starting 24 Februrary 2020) ###

# load libraries
library(tidyverse)
library(vegan)
# library(usedist)


#---> calculate distance matrix for core microbes
# function to get community matrix for each year
make_comm <- function(z,cols=7) {
  z=z[!duplicated(z$site_quadrat_id),]
  m=as.matrix(as.data.frame(select(z,-c(1:cols))))
  row.names(m) = z$site_quadrat_id
  return(m)
}
make_dist <- function(comm,diag=TRUE,upper=TRUE){
  df.dist=vegdist(comm,upper=upper)
  df.dist=as.matrix(df.dist, labels=TRUE)
  colnames(df.dist) <- rownames(df.dist) <- row.names(comm)
  df.dist <- df.dist[ order(colnames(df.dist)), order(colnames(df.dist)) ]
  return(as.dist(df.dist,upper=upper,diag=diag))
}

m16_core_family = read_csv("../../Data/data_parfrey/16S/family_16S_core_microbes.csv")
m16cf = m16_core_family %>% 
  group_split( year )
m16cf_comm <- lapply( m16cf, make_comm )
m16cf_dist <- lapply( m16cf_comm, make_dist )
m16_core_genus  = read_csv("../../Data/data_parfrey/16S/genus_16S_core_microbes.csv")
m16cg = m16_core_genus %>% 
  group_split( year )
m16cg_comm <- lapply( m16cg, make_comm )
m16cg_dist <- lapply( m16cg_comm, make_dist )
m18_core_family = read_csv("../../Data/data_parfrey/18S/family_18S_core_microbes.csv")
m18cf = m18_core_family %>% 
  group_split( year )
m18cf_comm <- lapply( m18cf, make_comm, cols=6 )
m18cf_dist <- lapply( m18cf_comm, make_dist )
m18_core_genus  = read_csv("../../Data/data_parfrey/18S/genus_18S_core_microbes.csv")
m18cg = m18_core_genus %>% 
  group_split( year )
m18cg_comm <- lapply( m18cg, make_comm, cols=6 )
m18cg_dist <- lapply( m18cg_comm, make_dist )
years = unique(m16_core_family$year)
#---> 

#### GRAZER DATA
###################################
############### 2015 ##############
###################################
# load grazers metadata
metadata_grazers_2015 <- read_csv("../output_data/2015_grazer_metadata.csv")
# change letters to numbers
split2015 <- strsplit( metadata_grazers_2015$sample, split = "_" )
metadata_grazers_2015$letts <- tolower(unlist(lapply( split2015, function(z) z[length(z)] )))
metadata_grazers_2015$site <- unlist(lapply( split2015, function(z) paste(z[-length(z)],collapse="_") ))
# read table for 2015 data
qc2015 <- read_csv( "../../Data/data_oconnor/grazers/2015_quadrat_sample_match_old_names.csv" )
# remove NA rows for sites
qc2015 <- filter( qc2015,!is.na(site) )
# change names in this table to match names decided upon by the whole group
qc2015$site[ qc2015$site == "Goose N" ] <- "goose_north"
qc2015$site[ qc2015$site == "McMullins N" ] <- "mcmullins_north"
qc2015$site[ qc2015$site == "McMullins S" ] <- "mcmullins_south"
qc2015$site[ qc2015$site == "Goose W" ] <- "goose_south_west"
qc2015$site[ qc2015$site == "Goose E" ] <- "goose_south_east"
qc2015$site[ qc2015$site == "Triquet N" ] <- "triquet_north"
qc2015$site[ qc2015$site == "Triquet S" ] <- "triquet_south"
qc2015$site[ qc2015$site == "Sandspit" ] <- "choked_sandspit"
qc2015$site[ qc2015$site == "Choked Pass, S. Pigu" ] <- "choked_inner"

metadata_grazers_2015 <- left_join(metadata_grazers_2015,qc2015,by=c("site","letts"="sample"))
metadata_grazers_2015$sample <- with(metadata_grazers_2015, paste(site,quad_number,sep="_") )
# load grazers bray curtis dissimilarity matrix
df_grazers_2015 <- read_csv("../output_data/2015_grazer_braycurtis.csv")
colnames(df_grazers_2015) <- rownames(df_grazers_2015) <- as.vector(metadata_grazers_2015$sample)
df_grazers_2015 <- df_grazers_2015[ order(colnames(df_grazers_2015)),order(colnames(df_grazers_2015)) ] 
# df_grazers_2015 <- as.dist(df_grazers_2015,upper=TRUE,diag=TRUE)
###################################
############### 2016 ##############
###################################
# load grazers metadata
metadata_grazers_2016 <- read_csv("../output_data/2016_grazer_metadata.csv")
# load grazers bray curtis dissimilarity matrix
df_grazers_2016 <- read_csv("../output_data/2016_grazer_braycurtis.csv")
colnames(df_grazers_2016) <- rownames(df_grazers_2016) <- as.vector(metadata_grazers_2016$sample)
df_grazers_2016 <- df_grazers_2016[ order(colnames(df_grazers_2016)),order(colnames(df_grazers_2016)) ] 
# df_grazers_2016 <- as.dist(df_grazers_2016,upper=TRUE,diag=TRUE)
###################################
############### 2017 ##############
###################################
# load grazers metadata
metadata_grazers_2017 <- read_csv("../output_data/2017_grazer_metadata.csv")
# load grazers bray curtis dissimilarity matrix
df_grazers_2017 <- read_csv("../output_data/2017_grazer_braycurtis.csv")
colnames(df_grazers_2017) <- rownames(df_grazers_2017) <- as.vector(metadata_grazers_2017$sample)
df_grazers_2017 <- df_grazers_2017[ order(colnames(df_grazers_2017)),order(colnames(df_grazers_2017)) ] 
# df_grazers_2017 <- as.dist(df_grazers_2017,upper=TRUE,diag=TRUE)




###################################
############### 2015 ##############
###################################

### Match samples to 2015 grazers
d16f15 <- m16cf_dist[[which(years==2015)]]
d16g15 <- m16cg_dist[[which(years==2015)]]
d18f15 <- m18cf_dist[[which(years==2015)]]
d18g15 <- m18cg_dist[[which(years==2015)]]

# which samples are common to both matrices?
labels_grazers_2015 <- as.vector(metadata_grazers_2015$sample)
labels_16S_2015 <- colnames(as.matrix(d16f15))
labels_18S_2015 <- colnames(as.matrix(d18f15))
common_samples16s <- intersect(labels_grazers_2015, labels_16S_2015)
common_samples18s <- intersect(labels_grazers_2015, labels_18S_2015)
common_samples1618 <- intersect(labels_16S_2015, labels_18S_2015)

# select the appropriate rows and columns in the dist objects
common_cols16s <- which( names(df_grazers_2015) %in% common_samples16s )
common_cols18s <- which( names(df_grazers_2015) %in% common_samples18s )
graz15_dist16s <- as.dist( df_grazers_2015[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
graz15_dist18s <- as.dist( df_grazers_2015[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
common_cols16s <- which( rownames(as.matrix(d16f15)) %in% common_samples16s )
common_cols18s <- which( rownames(as.matrix(d18f15)) %in% common_samples18s )
common_cols1618 <- which( rownames(as.matrix(d16f15)) %in% common_samples1618 )
common_cols1816 <- which( rownames(as.matrix(d18f15)) %in% common_samples1618 )
f15_dist16s <- as.dist( as.matrix(d16f15)[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
g15_dist16s <- as.dist( as.matrix(d16g15)[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
f15_dist18s <- as.dist( as.matrix(d18f15)[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
g15_dist18s <- as.dist( as.matrix(d18f15)[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
f15_dist1618 <- as.dist( as.matrix(d16f15)[common_cols1618,common_cols1618], upper=TRUE, diag=TRUE )
f15_dist1816 <- as.dist( as.matrix(d18f15)[common_cols1816,common_cols1816], upper=TRUE, diag=TRUE )
g15_dist1618 <- as.dist( as.matrix(d16g15)[common_cols1618,common_cols1618], upper=TRUE, diag=TRUE )
g15_dist1816 <- as.dist( as.matrix(d18g15)[common_cols1816,common_cols1816], upper=TRUE, diag=TRUE )




###################################
############### 2016 ##############
###################################

# #load 16S microbial distance matrix GENUS
# df_16S_2016_genus <- read.csv("../mantel_microbes_grazers/genus_16S_2016_braycurtis.csv")
# df_16S_2016_genus <- df_16S_2016_genus %>% 
#   rename("sample" = "X")
# 
# #load 16S microbial distance matrix FAMILY
# df_16S_2016_family <- read.csv("../mantel_microbes_grazers/family_16S_2016_braycurtis.csv")
# df_16S_2016_family <- df_16S_2016_family %>% 
#   rename("sample" = "X")
# 
# #load 18S microbial distance matrix GENUS
# df_18S_2016_genus <- read.csv("../mantel_microbes_grazers/genus_18S_2016_braycurtis.csv")
# df_18S_2016_genus <- df_18S_2016_genus %>% 
#   rename("sample" = "X")
# 
# #load 18S microbial distance matrix FAMILY
# df_18S_2016_family <- read.csv("../mantel_microbes_grazers/family_18S_2016_braycurtis.csv")
# df_18S_2016_family <- df_18S_2016_family %>% 
#   rename("sample" = "X")

#--->
### Match samples to 2016 grazers
d16f16 <- m16cf_dist[[which(years==2016)]]
d16g16 <- m16cg_dist[[which(years==2016)]]
d18f16 <- m18cf_dist[[which(years==2016)]]
d18g16 <- m18cg_dist[[which(years==2016)]]

# which samples are common to both matrices?
labels_grazers_2016 <- as.vector(metadata_grazers_2016$sample)
labels_16S_2016 <- colnames(as.matrix(d16f16))
labels_18S_2016 <- colnames(as.matrix(d18f16))
common_samples16s <- intersect(labels_grazers_2016, labels_16S_2016)
common_samples18s <- intersect(labels_grazers_2016, labels_18S_2016)
common_samples1618 <- intersect(labels_16S_2016, labels_18S_2016)

# select the appropriate rows and columns in the dist objects
common_cols16s <- which( names(df_grazers_2016) %in% common_samples16s )
common_cols18s <- which( names(df_grazers_2016) %in% common_samples18s )
graz16_dist16s <- as.dist( df_grazers_2016[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
graz16_dist18s <- as.dist( df_grazers_2016[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
common_cols16s <- which( rownames(as.matrix(d16f16)) %in% common_samples16s )
common_cols18s <- which( rownames(as.matrix(d18f16)) %in% common_samples18s )
common_cols1618 <- which( rownames(as.matrix(d16f16)) %in% common_samples1618 )
common_cols1816 <- which( rownames(as.matrix(d18f16)) %in% common_samples1618 )
f16_dist16s <- as.dist( as.matrix(d16f16)[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
g16_dist16s <- as.dist( as.matrix(d16g16)[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
f16_dist18s <- as.dist( as.matrix(d18f16)[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
g16_dist18s <- as.dist( as.matrix(d18f16)[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
f16_dist1618 <- as.dist( as.matrix(d16f16)[common_cols1618,common_cols1618], upper=TRUE, diag=TRUE )
f16_dist1816 <- as.dist( as.matrix(d18f16)[common_cols1816,common_cols1816], upper=TRUE, diag=TRUE )
g16_dist1618 <- as.dist( as.matrix(d16g16)[common_cols1618,common_cols1618], upper=TRUE, diag=TRUE )
g16_dist1816 <- as.dist( as.matrix(d18g16)[common_cols1816,common_cols1816], upper=TRUE, diag=TRUE )





###################################
############### 2017 ##############
###################################

### Match samples to 2017 grazers
d16f17 <- m16cf_dist[[which(years==2017)]]
d16g17 <- m16cg_dist[[which(years==2017)]]
d18f17 <- m18cf_dist[[which(years==2017)]]
d18g17 <- m18cg_dist[[which(years==2017)]]

# which samples are common to both matrices?
labels_grazers_2017 <- as.vector(metadata_grazers_2017$sample)
labels_16S_2017 <- colnames(as.matrix(d16f17))
labels_18S_2017 <- colnames(as.matrix(d18f17))
common_samples16s <- intersect(labels_grazers_2017, labels_16S_2017)
# for now get rid of choked_inner_2 for grazers
common_samples16s <- common_samples16s[ -1]
common_samples18s <- intersect(labels_grazers_2017, labels_18S_2017)
common_samples1618 <- intersect(labels_16S_2017, labels_18S_2017)


# select the appropriate rows and columns in the dist objects
common_cols16s <- which( names(df_grazers_2017) %in% common_samples16s )
common_cols18s <- which( names(df_grazers_2017) %in% common_samples18s )
graz17_dist16s <- as.dist( df_grazers_2017[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
graz17_dist18s <- as.dist( df_grazers_2017[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
common_cols16s <- which( rownames(as.matrix(d16f17)) %in% common_samples16s )
common_cols18s <- which( rownames(as.matrix(d18f17)) %in% common_samples18s )
common_cols1618 <- which( rownames(as.matrix(d16f17)) %in% common_samples1618 )
common_cols1816 <- which( rownames(as.matrix(d18f17)) %in% common_samples1618 )
f17_dist16s <- as.dist( as.matrix(d16f17)[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
g17_dist16s <- as.dist( as.matrix(d16g17)[common_cols16s,common_cols16s], upper=TRUE, diag=TRUE )
f17_dist18s <- as.dist( as.matrix(d18f17)[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
g17_dist18s <- as.dist( as.matrix(d18f17)[common_cols18s,common_cols18s], upper=TRUE, diag=TRUE )
f17_dist1618 <- as.dist( as.matrix(d16f17)[common_cols1618,common_cols1618], upper=TRUE, diag=TRUE )
f17_dist1816 <- as.dist( as.matrix(d18f17)[common_cols1816,common_cols1816], upper=TRUE, diag=TRUE )
g17_dist1618 <- as.dist( as.matrix(d16g17)[common_cols1618,common_cols1618], upper=TRUE, diag=TRUE )
g17_dist1816 <- as.dist( as.matrix(d18g17)[common_cols1816,common_cols1816], upper=TRUE, diag=TRUE )






################################################
###### MANTEL TESTS ####################
#############################
####################
#############
#######
#


### 2015
## family level
# grazer-16s
m_graz_16s_fam_2015 <- mantel( graz15_dist16s,f15_dist16s, method = "spearman", permutations = 9999, na.rm = FALSE )
# grazer-18s
m_graz_18s_fam_2015 <- mantel( graz15_dist18s,f15_dist18s, method = "spearman", permutations = 9999, na.rm = FALSE )
# 16s-18s
m_16s_18s_fam_2015 <- mantel( f15_dist1618,f15_dist1816, method = "spearman", permutations = 9999, na.rm = FALSE )
## genus level
# grazer-16s
m_graz_16s_gen_2015 <- mantel( graz15_dist16s,g15_dist16s, method = "spearman", permutations = 9999, na.rm = FALSE )
# grazer-18s
m_graz_18s_gen_2015 <- mantel( graz15_dist18s,g15_dist18s, method = "spearman", permutations = 9999, na.rm = FALSE )
# 16s-18s
m_16s_18s_gen_2015 <- mantel( g15_dist1618,g15_dist1816, method = "spearman", permutations = 9999, na.rm = FALSE )

### 2016
## family level
# grazer-16s
m_graz_16s_fam_2016 <- mantel( graz16_dist16s,f16_dist16s, method = "spearman", permutations = 9999, na.rm = FALSE )
# grazer-18s
m_graz_18s_fam_2016 <- mantel( graz16_dist18s,f16_dist18s, method = "spearman", permutations = 9999, na.rm = FALSE )
# 16s-18s
m_16s_18s_fam_2016 <- mantel( f16_dist1618,f16_dist1816, method = "spearman", permutations = 9999, na.rm = FALSE )
## genus level
# grazer-16s
m_graz_16s_gen_2016 <- mantel( graz16_dist16s,g16_dist16s, method = "spearman", permutations = 9999, na.rm = FALSE )
# grazer-18s
m_graz_18s_gen_2016 <- mantel( graz16_dist18s,g16_dist18s, method = "spearman", permutations = 9999, na.rm = FALSE )
# 16s-18s
m_16s_18s_gen_2016 <- mantel( g16_dist1618,g16_dist1816, method = "spearman", permutations = 9999, na.rm = FALSE )

### 2017
## family level
# grazer-16s
m_graz_16s_fam_2017 <- mantel( graz17_dist16s,f17_dist16s, method = "spearman", permutations = 9999, na.rm = FALSE )
# grazer-18s
m_graz_18s_fam_2017 <- mantel( graz17_dist18s,f17_dist18s, method = "spearman", permutations = 9999, na.rm = FALSE )
# 16s-18s
m_16s_18s_fam_2017 <- mantel( f17_dist1618,f17_dist1816, method = "spearman", permutations = 9999, na.rm = FALSE )
## genus level
# grazer-16s
m_graz_16s_gen_2017 <- mantel( graz17_dist16s,g17_dist16s, method = "spearman", permutations = 9999, na.rm = FALSE )
# grazer-18s
m_graz_18s_gen_2017 <- mantel( graz17_dist18s,g17_dist18s, method = "spearman", permutations = 9999, na.rm = FALSE )
# 16s-18s
m_16s_18s_gen_2017 <- mantel( g17_dist1618,g17_dist1816, method = "spearman", permutations = 9999, na.rm = FALSE )


## combine all tests into a list
mantel_genus  <- list(m_graz_16s_gen_2015,m_graz_18s_gen_2015,m_16s_18s_gen_2015,
                      m_graz_16s_gen_2016,m_graz_18s_gen_2016,m_16s_18s_gen_2016,
                      m_graz_16s_gen_2017,m_graz_18s_gen_2017,m_16s_18s_gen_2017)
mantel_family <- list(m_graz_16s_fam_2015,m_graz_18s_fam_2015,m_16s_18s_fam_2015,
                      m_graz_16s_fam_2016,m_graz_18s_fam_2016,m_16s_18s_fam_2016,
                      m_graz_16s_fam_2017,m_graz_18s_fam_2017,m_16s_18s_fam_2017)


#### WRITE OUTPUTS ALTOGETHER ###
mantel_extract <- function(x) {
  data.frame(stat = x$statistic, p= x$signif)
}

mantel_family <- do.call( rbind, lapply( mantel_family, mantel_extract ))
mantel_family <- data.frame( year=rep(2015:2017,each=3),lhs=rep(c("grazer","grazer","16S"),3),rhs=rep(c("16S","18S","18S"),3), mantel_family )

mantel_genus <- do.call( rbind, lapply( mantel_genus, mantel_extract ))
mantel_genus <- data.frame( year=rep(2015:2017,each=3),lhs=rep(c("grazer","grazer","16S"),3),rhs=rep(c("16S","18S","18S"),3), mantel_genus )



write_csv( mantel_family, "mantel_tests_family.csv" )
write_csv( mantel_genus, "mantel_tests_genus.csv" )
