### Microeukaryotes pipeline phyloseq ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 02, 2020 ###

library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(ape)
library(data.table)
library(vegan)

#### Importing files ####
all_years_18S_unfiltered <- readRDS("/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/micro_eukaryotes/18S_allyears_unfiltered.RDS")

#### QUALITY FILTERING TAXA DATA ####

# 1. Remove unassigned taxa and zostera OTUs #these are from the host plant, we want only microbial community
all_years_18S_filtered <- all_years_18S_unfiltered %>%
  subset_taxa(Rank1 != "Unassigned"| is.na(Rank1))  %>% 
  subset_taxa(Rank6 != "Zostera"| is.na(Rank6)) 

# 2. Remove seaweeds and plants
all_years_18S_filtered  <- all_years_18S_filtered %>%
  subset_taxa(Rank3 != "Rhodophyceae"| is.na(Rank3))  %>% 
  subset_taxa(Rank5 != "Ulvophyceae"| is.na(Rank5))  %>% 
  subset_taxa(Rank5 != "Ulvales"| is.na(Rank5))  %>% 
  subset_taxa(Rank5 != "Embryophyta"| is.na(Rank5))  %>% 
  subset_taxa(Rank5 != "Phaeophyceae"| is.na(Rank5))  %>%
  subset_taxa(Rank5 != "Desmarestiales"| is.na(Rank5))  %>%
  subset_taxa(Rank4 != "Charophyta"| is.na(Rank4)) %>% 
  subset_taxa(Rank5 != "Ectocarpales"| is.na(Rank5))

# 3. Remove large-bodied Metazoans
all_years_18S_filtered  <- all_years_18S_filtered  %>%
  subset_taxa(Rank4 == "Metazoa_(Animalia)")  # | is.na(Rank4) removed because it was getting other stuff together
all_years_18S_filtered  <- all_years_18S_filtered  %>%
  subset_taxa(Rank6 != "Mammalia"| is.na(Rank6))

#View(as.data.frame(tax_table(all_years_18S_filtered )))

# # FILTERING per sample
# # 5. Remove ASVs with less than ~ 2-5 reads in a given sample - PER SAMPLE FILTERING
# otu <- as.data.frame(otu_table(all_years_18S_filtered)) #get ASV table
# otu_table(all_years_18S_filtered)[otu <= 3] <- 0 # for entries where the raw abundance of an ASV in a sample is less than 3 reads, set the raw read count to 0 - free of noise

# FILTERING overall
# 6. Remove OTUs with less than N total reads. (N = 250 in example) whole dataset
all_years_18S_filtered <- prune_taxa(taxa_sums(all_years_18S_filtered) >= 250, all_years_18S_filtered) 

# 7. Remove samples with less than N reads. (N = 1000 in example) whole dataset
all_years_18S_filtered <- prune_samples(sample_sums(all_years_18S_filtered) >= 1000, all_years_18S_filtered)

# 8. look at minimum, mean, and maximum sample counts, if desired
# smin <- 
#   min(sample_sums(all_years_18S_filtered))
# meanreads <- 
#   mean(sample_sums(all_years_18S_filtered))
# smax <- 
#   max(sample_sums(all_years_18S_filtered))
# totalreads <- 
#   sum(sample_sums(all_years_18S_filtered))

all_years_18S_filtered_meso <- all_years_18S_filtered %>% subset_samples(survey_type == "meso_quadrat") 

all_years_18S_filtered_meso_Zos <- all_years_18S_filtered_meso %>% subset_samples(sample_type =="leaf_old") 

# get rid of taxa that sums to zero
minTotRelAbun = 0
x = taxa_sums(all_years_18S_filtered_meso_Zos)
keepTaxa = (x / sum(x)) > minTotRelAbun
macro18S_final = prune_taxa(keepTaxa, all_years_18S_filtered_meso_Zos)

###################################################################
### normalizing data using coverage based rarefaction ASV LEVEL ###
###################################################################
library(devtools)
###install_github("vmikk/metagMisc")
###install_github('JohnsonHsieh/iNEXT')
library(metagMisc)
library(iNEXT)
### phyloseq_coverage_raref(physeq, coverage = NULL, iter = 1, replace = F, correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...)

### Samples standardized by size will have different degrees of completeness. When we compare samples with the same coverage, we are making sure that samples are equally complete and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).
### i.e. a seagrass sample with 10,000 total reads will have a different coverage than a seawater sample with 10,000 reads if seawater samples have many more species

### phyloseq object to be used
macro18S_final_for_CB <- macro18S_final

### the following steps are necessary before running the function because it needs the otu table in the phyloseq object to have taxa as rows

# check if taxa are rows in the phyloseq object
taxa_are_rows(macro18S_final_for_CB)

# transpose so taxa are rows
otu_table(macro18S_final_for_CB) <- t(otu_table(macro18S_final_for_CB))

taxa_are_rows(macro18S_final_for_CB)

# prepare otu table for function
x <- metagMisc::prepare_inext(
  as.data.frame(otu_table(macro18S_final_for_CB)),
  correct_singletons = T)

# check if read counts are correct (samples should show "numeric" in the second column)
SC <- plyr::llply(.data = x, .fun = function(z){ try( iNEXT:::Chat.Ind(z, sum(z)) ) })
plyr::ldply(.data = SC, .fun = class)

macro18S_COVERAGE_RAREF_5000 <- phyloseq_coverage_raref(physeq=macro18S_final_for_CB, coverage = 0.9, iter = 5000, replace = F, correct_singletons = F, drop_lowcoverage = T)

saveRDS(macro18S_COVERAGE_RAREF_5000, "/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/macro_18S/macro_18S_COVERAGE_RAREF_5000.rds")
macro18S_COVERAGE_RAREF_5000 <- readRDS("/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/macro_18S/macro_18S_COVERAGE_RAREF_5000.rds")

### Average otu tables from all iterations to get a final robust table ###
subset_phylo_objects_5000 <- macro18S_COVERAGE_RAREF_5000[c(1:5000)]
# first, extract otu tables from phyloseq objects
# this is how you do it for a single phyloseq object:
# y <- as.data.frame(t(phyloseq::otu_table(all_16S_COVERAGE_RAREF)))
# now do it for the list of phyloseq objects
otu_tables_5000 <- lapply(subset_phylo_objects_5000, function(z) as.data.frame(t(phyloseq::otu_table(z))))

# average all matrices to get the mean abundance across all iterations

average_otu_tables_5000 <- Reduce("+",otu_tables_5000)/length(otu_tables_5000)
# IMPORTANT! NEED TO ROUND IT OTHERWISE iNEXT WILL NOT WORK!
average_otu_tables_5000_round <- average_otu_tables_5000 %>% mutate_at(vars(starts_with("ASV")), funs(round(., 0)))
# add SampleID column back
average_otu_tables_5000_round$SampleID <- rownames(average_otu_tables_5000) 
average_otu_tables_5000_round <- average_otu_tables_5000_round %>% 
  select(SampleID, everything())

write.csv(average_otu_tables_5000_round, "/Users/parfreylab/Desktop/lab_member_files/bia/Calvert_O-Connor_eelgrass/Data/macro_18S/macro_18S_average_otu_tables_5000.csv", quote=F, row.names=F )
