### Species richness estimates iNEXT (interpolation extrapolation) ###
### Tutorial by https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
### Author: Bianca Trevizan Segovia ###
### Date created: July 20th, 2020 ###

# ## install iNEXT package from CRAN
# install.packages("iNEXT")
# 
# ## install the latest version from github
# install.packages('devtools')
# library(devtools)
# install_github('JohnsonHsieh/iNEXT')

## import packages
library(iNEXT)
library(ggplot2)
library(dplyr)
library(tidyverse)

#########################################
############ 16S prokaryotes ############
#########################################

###########################################
### WITH COVERAGE-BASED NORMALIZED DATA ###
###########################################

#######################################################
############# TWO HUNDRED ITERATIONS  #################
#######################################################

### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_200_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_ASV)[1:17]

data_all_16S <- microbes_16S_ASV %>%
  dplyr::select(., -(c(1:16))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2015 <- microbes_16S_ASV %>%
  filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:16))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2016 <- microbes_16S_ASV %>%
  filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:16))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2017 <- microbes_16S_ASV %>%
  filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:16))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2018 <- microbes_16S_ASV %>%
  filter(year == "2018") %>% # filter by your treatment
  dplyr::select(., -(c(1:16))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

#create data input for rarefied richness estimates: 
input_2015 <- data_2015[,2] 
input_2016 <- data_2016[,2]
input_2017 <- data_2017[,2]
input_2018 <- data_2018[,2]


input_all <- cbind(input_2015, input_2016, input_2017, input_2018)

chao_richness_16S <- ChaoRichness(input_all, datatype="abundance")
chao_richness_16S

chao_richness_16S_all_years <- data_all_16S %>% 
  select(N)  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_16S_all_years
#Observed Estimator Est_s.e. 95% Lower 95% Upper
#N     1522  1867.954   46.481  1788.174  1971.646

estimate_richness_16S <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) # q=0 estimates richness

coverage_based_16S <- ggiNEXT(estimate_richness_16S, type=3)
coverage_based_16S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_200_CB_16S.png", plot = coverage_based_16S, width=250, height=200, units="mm",dpi=300)

sample_size_based_16S <- ggiNEXT(estimate_richness_16S, type=1)
sample_size_based_16S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_200_SSB_16S.png", plot = sample_size_based_16S, width=250, height=200, units="mm",dpi=300)


#########################################
########## 18S microeukaryotes ##########
#########################################

###########################################
### WITH COVERAGE-BASED NORMALIZED DATA ###
###########################################

##################################################
############# TWO HUNDRED ITERATIONS #############
##################################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_200_COVERAGE_RAREF.csv", header=T)

names(microbes_18S_ASV)[1:12]

data_all_18S <- microbes_18S_ASV %>%
  dplyr::select(., -(c(1:10))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2015 <- microbes_18S_ASV %>%
  filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:10))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2016 <- microbes_18S_ASV %>%
  filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:10))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2017 <- microbes_18S_ASV %>%
  filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:10))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2018 <- microbes_18S_ASV %>%
  filter(year == "2018") %>% # filter by your treatment
  dplyr::select(., -(c(1:10))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

#create data input for rarefied richness estimates: 
input_2015 <- data_2015[,2] 
input_2016 <- data_2016[,2]
input_2017 <- data_2017[,2]
input_2018 <- data_2018[,2]


input_all <- cbind(input_2015, input_2016, input_2017, input_2018)

chao_richness_18S <- ChaoRichness(input_all, datatype="abundance")
chao_richness_18S

chao_richness_18S_all_years <- data_all_18S %>% 
  select(N)  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_18S_all_years
#Observed Estimator Est_s.e. 95% Lower 95% Upper
#N  231   353.685    34.37   302.585   441.264

estimate_richness_18S <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) # q=0 estimates richness

coverage_based_18S <- ggiNEXT(estimate_richness_18S, type=3)
coverage_based_18S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_200_CB_18S.png", plot = coverage_based_18S, width=250, height=200, units="mm",dpi=300)

sample_size_based_18S <- ggiNEXT(estimate_richness_18S, type=1)
sample_size_based_18S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_200_SSB_18S.png", plot = sample_size_based_18S, width=250, height=200, units="mm",dpi=300)


#########################################
############ Macroeukaryotes ############
#########################################

#### Finest level ####

### Read table metadata and abundances
m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan
# define regions as first word of site
m$region <- unlist(lapply( strsplit(m$site,split = "_"), function(z) z[1]))
# unique sample ID to differential samples from different sites
m$ID <- with(m, paste(site,sample,sep = "_"))
# change year to character
m$group <- paste0( "year",m$year )

# filter taxa and sites # taxon 2 = finest
mfilt_finest <- m %>%
  dplyr::select( year, group, region, site, sample, ID, taxon = taxon2, remove, size ) %>% 
  dplyr::filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum_finest <- mfilt_finest %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean_finest <- m.sum_finest %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance)) 

# make a community dataset
m.meta_finest <- m.sum_finest %>% 
  spread( taxon, abundance, fill=0 ) 

names(m.meta_finest)[1:16]

m.meta_finest <- m.meta_finest %>% 
  ungroup( sample, ID,year, group, region, site )

data_all_inverts <- as.data.frame(m.meta_finest %>%
                             dplyr::select(., -(c(1:6))) %>% # remove metadata
                             summarise_all(., funs(sum)) %>% # sum abundances across all replicates
                             tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2014 <- as.data.frame(m.meta_finest %>%
  dplyr::filter(year == "2014") %>% # filter by your treatment
  dplyr::select(., -(c(1:6))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2015 <- as.data.frame(m.meta_finest %>%
  dplyr::filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:6))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2016 <- as.data.frame(m.meta_finest %>%
  dplyr::filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:6))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2017 <- as.data.frame(m.meta_finest %>%
  dplyr::filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:6))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

#create data input for rarefied richness estimates: 
input_2014 <- data_2014[,2]
input_2015 <- data_2015[,2] 
input_2016 <- data_2016[,2]
input_2017 <- data_2017[,2]

input_all <- cbind(input_2014,input_2015, input_2016, input_2017) 

chao_richness_inverts <- ChaoRichness(input_all, datatype="abundance")
chao_richness_inverts

chao_richness_inverts_all_years <- data_all_inverts %>% 
  dplyr::select(N)  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_inverts_all_years
#  Observed Estimator Est_s.e. 95% Lower 95% Upper
#N      105       109    3.742   105.846    123.92

estimate_richness_inverts <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) # q=0 estimates richness

coverage_based_inverts <- ggiNEXT(estimate_richness_inverts, type=3)
coverage_based_inverts
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_inverts.png", plot = coverage_based_inverts, width=250, height=200, units="mm",dpi=300)

sample_size_based_inverts <- ggiNEXT(estimate_richness_inverts, type=1)
sample_size_based_inverts
ggsave("R_Code_and_Analysis/gammadiversity/sample_size_based_raref_inverts.png", plot = sample_size_based_inverts, width=250, height=200, units="mm",dpi=300)
