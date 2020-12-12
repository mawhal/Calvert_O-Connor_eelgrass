### Inverts pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ### *Katy Davis provided the base code for creating the phyloseq object
### Date created: October 22nd, 2020 ###

######################
#### Finest level ####
######################

#### Create phyloseq object ####
library(phyloseq)
library(tidyverse)

# Read in file
inverts <- read.csv(file="Data/R_Code_for_Data_Prep/master_data/MASTER_grazers_finest_to_phyloseq.csv")

# add unique sample labels
inverts <- inverts %>% 
  dplyr::mutate(ID_year = paste(ID, year, sep = "_"))

inverts <- inverts %>% dplyr::select(ID_year, everything())

# Separate to metadata and species table
# metadata
invert.meta <- inverts %>% dplyr::select(ID_year,ID, year, group, region, site, sample) 

# taxa table
invert.otus <- inverts %>% dplyr::select(ID_year, Actiniidae:ncol(inverts))

# Make sure species data frame is numeric matrix and transpose so taxa are rows
invert.otus_t <- setNames(as.data.frame(t(invert.otus[,-1])), invert.otus[,1])
invert.otus_t_matrix <- data.matrix(invert.otus_t)

# Desired format for species table
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30

# Make sure metadata has ID as rownames
invert.meta.ID <- data.frame(invert.meta[,-1], row.names=invert.meta[,1])

# Here you can create a taxonomy table if you want (could do in excel?) in the matrix format
# taxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
taxonomy_finest <- read.csv( "Data/macro_eukaryotes/taxonomy_table_inverts.csv", header=TRUE)
taxonomy_finest_f <- taxonomy_finest %>% 
  column_to_rownames(var = "finest")
tax.mat <- as.matrix(taxonomy_finest_f)

##       Domain Phylum Class Order Family Genus Species
## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o" 
# You don't have to have those specific taxonomic ranks, you could have any number of columns and can be called anything

#### Save each object in "phyloseq" format to be combined in a phyloseq object ####
OTU <- otu_table(invert.otus_t_matrix, taxa_are_rows = T)
TAX <- tax_table(tax.mat)
META <- sample_data(invert.meta.ID)

# invert.phyloseq <- phyloseq(OTU, TAX, META)
invert.phyloseq <- phyloseq(OTU, TAX,META)
sample_names(OTU)
sample_names(META)

saveRDS(invert.phyloseq, "Data/macro_eukaryotes/invert_phyloseq.rds")

#########################################################
### normalizing data using coverage based rarefaction ###
#########################################################
library(devtools)
#install_github("vmikk/metagMisc")
#install_github('JohnsonHsieh/iNEXT')
library(metagMisc)
library(iNEXT)
# phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

# Samples standardized by size will have different degrees of completeness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).

# phyloseq object to be used
invert_phyloseq <- readRDS("Data/macro_eukaryotes/invert_phyloseq.rds")
View(as.data.frame(otu_table(invert_phyloseq)))
# check if taxa are rows in the phyloseq object
taxa_are_rows(invert_phyloseq)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(invert_phyloseq)),
  correct_singletons = F) ### IMPORTANT! DO NOT CORRECT FOR SINGLETONS!!!

# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

inverts_otu_table <- as.data.frame(otu_table(invert_phyloseq))
View(inverts_otu_table)
#Due to the stochasticity introduced in random subsampling results could be slightly different.
#So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.
# 
# # run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)

invert_COVERAGE_RAREF_1000 <- phyloseq_coverage_raref(physeq=invert_phyloseq, coverage = 0.8, iter = 1000, replace = F, correct_singletons = F, drop_lowcoverage = T)   ### IMPORTANT! DO NOT CORRECT FOR SINGLETONS!!!

#Warning message:
#  In phyloseq_coverage_raref(physeq = invert_phyloseq, coverage = 0.8,  :
#                              Samples with coverage lower than the selected threshold #were discarded (n = 6).

saveRDS(invert_COVERAGE_RAREF_1000, "Data/macro_eukaryotes/invert_COVERAGE_RAREF_1000.rds")

invert_COVERAGE_RAREF_1000 <- readRDS("Data/macro_eukaryotes/invert_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_1000 <- invert_COVERAGE_RAREF_1000[c(1:1000)]

# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
otu_tables_1000 <- lapply(subset_phylo_objects_1000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# # average all matrices to get the mean abundance across all iterations
average_otu_tables_1000 <- Reduce("+",otu_tables_1000)/length(otu_tables_1000)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round <- average_otu_tables_1000 %>% dplyr::mutate_at(c(1:ncol(average_otu_tables_1000)), funs(round(., 0)))
View(as.data.frame(average_otu_tables_1000_round))
# add SampleID column back
average_otu_tables_1000_round$ID_year <- rownames(average_otu_tables_1000)
average_otu_tables_1000_round <- average_otu_tables_1000_round %>%
  dplyr::select(ID_year, everything())

write.csv(average_otu_tables_1000_round, "Data/macro_eukaryotes/inverts_average_otu_tables_1000.csv", quote=F, row.names=F )

inverts_1000 <- read.csv("Data/macro_eukaryotes/inverts_average_otu_tables_1000.csv")
View(invert.meta)
inverts_1000_metadata <- inner_join(invert.meta, inverts_1000, by = "ID_year")

write.csv(inverts_1000_metadata, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv", quote=F, row.names=F)

######################
#### Family level ####
######################

#### Create phyloseq object ####
library(phyloseq)
library(tidyverse)

# Read in file
inverts_family <- read.csv(file="Data/R_Code_for_Data_Prep/master_data/MASTER_grazers_family_to_phyloseq.csv")

# add unique sample labels
inverts_family <- inverts_family %>% 
  dplyr::mutate(ID_year = paste(ID, year, sep = "_"))

inverts_family <- inverts_family %>% dplyr::select(ID_year, everything())

# Separate to metadata and species table
# metadata
invert.meta_family <- inverts_family %>% dplyr::select(ID_year,ID, year, group, region, site, sample) 

# taxa table
invert.otus_family <- inverts_family %>% dplyr::select(ID_year, Actiniidae:ncol(inverts_family))

# Make sure species data frame is numeric matrix and transpose so taxa are rows
invert.otus_t_family <- setNames(as.data.frame(t(invert.otus_family[,-1])), invert.otus_family[,1])
invert.otus_t_matrix_family <- data.matrix(invert.otus_t_family)

# Desired format for species table
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30

# Make sure metadata has ID as rownames
invert.meta.ID_family <- data.frame(invert.meta_family[,-1], row.names=invert.meta_family[,1])

# Here you can create a taxonomy table if you want (could do in excel?) in the matrix format
# tax.mat <- matrix(taxa.data)

##       Domain Phylum Class Order Family Genus Species
## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o" 
# You don't have to have those specific taxonomic ranks, you could have any number of columns and can be called anything

#### Save each object in "phyloseq" format to be combined in a phyloseq object ####
OTU_family <- otu_table(invert.otus_t_matrix_family, taxa_are_rows = T)
#TAX <- tax_table(tax.mat)
META_family <- sample_data(invert.meta.ID_family)

# invert.phyloseq <- phyloseq(OTU, TAX, META)
invert.phyloseq_family <- phyloseq(OTU_family, META_family)
sample_names(OTU_family)
sample_names(META_family)

saveRDS(invert.phyloseq_family, "Data/macro_eukaryotes/invert_phyloseq_family.rds")

#########################################################
### normalizing data using coverage based rarefaction ###
#########################################################
library(devtools)
#install_github("vmikk/metagMisc")
#install_github('JohnsonHsieh/iNEXT')
library(metagMisc)
library(iNEXT)
# phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

# Samples standardized by size will have different degrees of completeness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).

# phyloseq object to be used
invert_phyloseq_family <- readRDS("Data/macro_eukaryotes/invert_phyloseq_family.rds")
#View(as.data.frame(otu_table(invert_phyloseq_family)))
# check if taxa are rows in the phyloseq object
taxa_are_rows(invert_phyloseq_family)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(invert_phyloseq_family)),
  correct_singletons = F) ### IMPORTANT! DO NOT CORRECT FOR SINGLETONS!!!

# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

inverts_otu_table_family <- as.data.frame(otu_table(invert_phyloseq_family))
#View(inverts_otu_table)
#Due to the stochasticity introduced in random subsampling results could be slightly different.
#So you have to average diversity estimates or sample dissimilarities across multiple rarefactions.
# 
# # run coverage-based rarefaction (Chao & Jost, 2012) correcting for singletons (Chiu & Chao 2016)

invert_COVERAGE_RAREF_1000_family <- phyloseq_coverage_raref(physeq=invert_phyloseq_family, coverage = 0.8, iter = 1000, replace = F, correct_singletons = F, drop_lowcoverage = T)   ### IMPORTANT! DO NOT CORRECT FOR SINGLETONS!!!

#Warning message:
#  In phyloseq_coverage_raref(physeq = invert_phyloseq, coverage = 0.8,  :
#                              Samples with coverage lower than the selected threshold #were discared (n = 6).

saveRDS(invert_COVERAGE_RAREF_1000_family, "Data/macro_eukaryotes/invert_family_COVERAGE_RAREF_1000.rds")

invert_COVERAGE_RAREF_1000_family <- readRDS("Data/macro_eukaryotes/invert_family_COVERAGE_RAREF_1000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_1000_family <- invert_COVERAGE_RAREF_1000_family[c(1:1000)]

# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
otu_tables_1000_family <- lapply(subset_phylo_objects_1000_family, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# # average all matrices to get the mean abundance across all iterations
average_otu_tables_1000_family <- Reduce("+",otu_tables_1000_family)/length(otu_tables_1000_family)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_1000_round_family <- average_otu_tables_1000_family %>% dplyr::mutate_at(c(1:ncol(average_otu_tables_1000_family)), funs(round(., 0)))
#View(as.data.frame(average_otu_tables_1000_round))
# add SampleID column back
average_otu_tables_1000_round_family$ID_year <- rownames(average_otu_tables_1000_family)
average_otu_tables_1000_round_family <- average_otu_tables_1000_round_family %>%
  dplyr::select(ID_year, everything())

write.csv(average_otu_tables_1000_round_family, "Data/macro_eukaryotes/inverts_family_average_otu_tables_1000.csv", quote=F, row.names=F )

inverts_1000_family <- read.csv("Data/macro_eukaryotes/inverts_family_average_otu_tables_1000.csv")
View(invert.meta_family)
inverts_1000_metadata_family <- left_join(invert.meta_family, inverts_1000_family, by = "ID_year")

write.csv(inverts_1000_metadata_family, "/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_family_1000_COVERAGE_RAREF.csv", quote=F, row.names=F)