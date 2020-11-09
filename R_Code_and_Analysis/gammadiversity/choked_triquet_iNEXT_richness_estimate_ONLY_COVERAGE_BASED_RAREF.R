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
############# ONE THOUSAND ITERATIONS  #################
#######################################################

### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_ASV)[1:17]

microbes_16S_ASV <- microbes_16S_ASV %>% 
  filter(region == "choked" | region == "triquet")

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
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# input_2015      649   844.986   37.257   784.471   932.533
# input_2016      368   562.935   47.105   490.212   678.933
# input_2017      577  1020.642   78.226   891.840  1202.138
# input_2018      494   615.310   25.490   574.718   676.316

chao_richness_16S_all_years <- data_all_16S %>% 
  select(N)  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_16S_all_years
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# N     1129  1540.453   57.648  1442.068  1669.757

estimate_richness_16S <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) # q=0 estimates richness

coverage_based_16S <- ggiNEXT(estimate_richness_16S, type=3)
coverage_based_16S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/choked_triquet_cov_based_raref_1000_CB_16S.png", plot = coverage_based_16S, width=250, height=200, units="mm",dpi=300)

sample_size_based_16S <- ggiNEXT(estimate_richness_16S, type=1)
sample_size_based_16S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/choked_triquet_cov_based_raref_1000_SSB_16S.png", plot = sample_size_based_16S, width=250, height=200, units="mm",dpi=300)

#########################################
########## 18S microeukaryotes ##########
#########################################

###########################################
### WITH COVERAGE-BASED NORMALIZED DATA ###
###########################################

##################################################
############# ONE THOUSAND ITERATIONS #############
##################################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)

microbes_18S_ASV <- microbes_18S_ASV %>% 
  filter(region == "choked" | region == "triquet")

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
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# input_2015       49   123.741   41.840    75.862   256.963
# input_2016       47    83.023   23.849    58.051   164.421
# input_2017       31    62.797   23.181    39.845   145.314
# input_2018       62   110.877   26.245    80.225   193.081

chao_richness_18S_all_years <- data_all_18S %>% 
  select(N)  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_18S_all_years
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# N      119   190.219   25.769   154.812   260.633

estimate_richness_18S <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) # q=0 estimates richness

coverage_based_18S <- ggiNEXT(estimate_richness_18S, type=3)
coverage_based_18S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/choked_triquet_cov_based_raref_1000_CB_18S.png", plot = coverage_based_18S, width=250, height=200, units="mm",dpi=300)

sample_size_based_18S <- ggiNEXT(estimate_richness_18S, type=1)
sample_size_based_18S
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/choked_triquet_cov_based_raref_1000_SSB_18S.png", plot = sample_size_based_18S, width=250, height=200, units="mm",dpi=300)


#########################################
############ Macroeukaryotes ############
#########################################

#### Finest level ####

### Read table metadata and abundances
inverts <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_1000_COVERAGE_RAREF.csv")
names(inverts)

inverts <- inverts %>% 
  filter(region == "choked" | region == "triquet")

data_all_inverts <- as.data.frame(inverts %>%
                             dplyr::select(., -(c(1:7))) %>% # remove metadata
                             summarise_all(., funs(sum)) %>% # sum abundances across all replicates
                             tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2014 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2014") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2015 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2016 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise_all(., funs(sum)) %>% # sum abundances across all replicates
  tidyr::gather(., "species", "N")) %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2017 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
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
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# input_2014       15    32.919   23.517    17.520   142.400
# input_2015       16    17.117    1.759    16.125    25.967
# input_2016       26    50.436   31.028    29.580   192.778
# input_2017       28    29.330    1.842    28.175    38.112

chao_richness_inverts_all_years <- data_all_inverts %>% 
  dplyr::select(N)  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_inverts_all_years
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# N       46    59.489    12.45    48.896   108.829

estimate_richness_inverts <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50) # q=0 estimates richness

coverage_based_inverts <- ggiNEXT(estimate_richness_inverts, type=3)
coverage_based_inverts
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/choked_triquet_cov_based_raref_1000_CB_inverts.png", plot = coverage_based_inverts, width=250, height=200, units="mm",dpi=300)

sample_size_based_inverts <- ggiNEXT(estimate_richness_inverts, type=1)
sample_size_based_inverts
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/choked_triquet_cov_based_raref_1000_SSB_inverts.png", plot = sample_size_based_inverts, width=250, height=200, units="mm",dpi=300)
