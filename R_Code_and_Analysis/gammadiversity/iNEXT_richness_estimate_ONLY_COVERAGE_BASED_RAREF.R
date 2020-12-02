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
library(ggpubr)

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

data_all_16S <- microbes_16S_ASV %>%
  dplyr::select(-(c(1:15))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2015 <- microbes_16S_ASV %>%
  filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:15))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2016 <- microbes_16S_ASV %>%
  filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:15))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2017 <- microbes_16S_ASV %>%
  filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:15))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2018 <- microbes_16S_ASV %>%
  filter(year == "2018") %>% # filter by your treatment
  dplyr::select(., -(c(1:15))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

#create data input for rarefied richness estimates: 
input_2015 <- data_2015[,2] 
colnames(input_2015)[1] <- "input_2015"
input_2016 <- data_2016[,2]
colnames(input_2016)[1] <- "input_2016"
input_2017 <- data_2017[,2]
colnames(input_2017)[1] <- "input_2017"
input_2018 <- data_2018[,2]
colnames(input_2018)[1] <- "input_2018"

input_all <- as.data.frame(bind_cols(input_2015, input_2016, input_2017, input_2018))

chao_richness_16S <- ChaoRichness(input_all, datatype="abundance")
chao_richness_16S
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# input_2015      918  1085.268   31.250  1034.345  1158.479
# input_2016      559   799.622   47.239   723.363   911.262
# input_2017      814  1304.774   73.842  1180.038  1472.017
# input_2018      575   730.854   31.153   680.740   804.720

chao_richness_16S_all_years <- as.data.frame(data_all_16S %>% 
  dplyr::select(N))  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_16S_all_years
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# N     1501   1879.15   50.851  1791.879  1992.605

# endpoint NULL = double the reference sample size.
estimate_richness_16S <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=500) # q=0 estimates richness

coverage_based_16S <- ggiNEXT(estimate_richness_16S, type=3)
coverage_based_16S <- coverage_based_16S + ggtitle("Prokaryotes")
coverage_based_16S <- coverage_based_16S + ggtitle("Prokaryotes")
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_1000_CB_16S.png", plot = coverage_based_16S, width=250, height=200, units="mm",dpi=300)

sample_size_based_16S <- ggiNEXT(estimate_richness_16S, type=1)
sample_size_based_16S <- sample_size_based_16S + ggtitle("Prokaryotes")
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_1000_SSB_16S.png", plot = sample_size_based_16S, width=250, height=200, units="mm",dpi=300)

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
names(microbes_18S_ASV)[1:12]

data_all_18S <- microbes_18S_ASV %>%
  dplyr::select(., -(c(1:9))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2015 <- microbes_18S_ASV %>%
  filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:9))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2016 <- microbes_18S_ASV %>%
  filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:9))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2017 <- microbes_18S_ASV %>%
  filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:9))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

data_2018 <- microbes_18S_ASV %>%
  filter(year == "2018") %>% # filter by your treatment
  dplyr::select(., -(c(1:9))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = starts_with("ASV"), names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N)) # descending order

#create data input for rarefied richness estimates: 
input_2015 <- data_2015[,2] 
colnames(input_2015)[1] <- "input_2015"
input_2016 <- data_2016[,2]
colnames(input_2016)[1] <- "input_2016"
input_2017 <- data_2017[,2]
colnames(input_2017)[1] <- "input_2017"
input_2018 <- data_2018[,2]
colnames(input_2018)[1] <- "input_2018"

input_all <- as.data.frame(bind_cols(input_2015, input_2016, input_2017, input_2018))

chao_richness_18S <- ChaoRichness(input_all, datatype="abundance")
chao_richness_18S
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# input_2015       76   140.725   30.808   102.694   232.938
# input_2016      138   199.887   23.338   168.284   264.469
# input_2017       59    84.905   13.445    68.946   126.469
# input_2018       94   185.915   40.136   134.531   302.441

chao_richness_18S_all_years <- as.data.frame(data_all_18S %>% 
  select(N))  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_18S_all_years
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# N      225     317.6   26.758   278.157   386.311

# endpoint NULL = double the reference sample size.
estimate_richness_18S <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=500) # q=0 estimates richness

coverage_based_18S <- ggiNEXT(estimate_richness_18S, type=3)
coverage_based_18S <- coverage_based_18S + ggtitle("Microeukaryotes")
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_1000_CB_18S.png", plot = coverage_based_18S, width=250, height=200, units="mm",dpi=300)

sample_size_based_18S <- ggiNEXT(estimate_richness_18S, type=1)
sample_size_based_18S <- sample_size_based_18S + ggtitle("Microeukaryotes")
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_1000_SSB_18S.png", plot = sample_size_based_18S, width=250, height=200, units="mm",dpi=300)

#########################################
############ Macroeukaryotes ############
#########################################

#### Finest level ####

### Read table metadata and abundances
inverts <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")
names(inverts)

data_all_inverts <- as.data.frame(inverts %>%
    dplyr::select(., -(c(1:7))) %>% # remove metadata
    summarise(across(everything(), sum)) %>% # sum abundances across all replicates
    pivot_longer(cols = Actiniidae:Thorlaksonius.subcarinatus, names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
    arrange(., desc(N))) # descending order

data_2014 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2014") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = Actiniidae:Thorlaksonius.subcarinatus, names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N))) # descending order

data_2015 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2015") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = Actiniidae:Thorlaksonius.subcarinatus, names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N))) # descending order

data_2016 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2016") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = Actiniidae:Thorlaksonius.subcarinatus, names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N))) # descending order

data_2017 <- as.data.frame(inverts %>%
  dplyr::filter(year == "2017") %>% # filter by your treatment
  dplyr::select(., -(c(1:7))) %>% # remove metadata
  summarise(across(everything(), sum)) %>% # sum abundances across all replicates
  pivot_longer(cols = Actiniidae:Thorlaksonius.subcarinatus, names_to = "species", values_to ="N") %>% # transform into 2 columns 1:species 2:N(summed abundances)
  arrange(., desc(N))) # descending order

#create data input for rarefied richness estimates: 
input_2014 <- as.data.frame(data_2014[,2])
colnames(input_2014)[1] <- "input_2014"
input_2015 <- as.data.frame(data_2015[,2]) 
colnames(input_2015)[1] <- "input_2015"
input_2016 <- as.data.frame(data_2016[,2])
colnames(input_2016)[1] <- "input_2016"
input_2017 <- as.data.frame(data_2017[,2])
colnames(input_2017)[1] <- "input_2017"

input_all <- as.data.frame(bind_cols(input_2014,input_2015, input_2016, input_2017))

chao_richness_inverts <- ChaoRichness(input_all, datatype="abundance")
chao_richness_inverts
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# input_2014       23    67.882   30.100    36.593   171.190
# input_2015       23    26.591    3.843    23.648    42.898
# input_2016       35    58.957   20.164    40.708   135.543
# input_2017       30    31.996    2.641    30.278    44.336

chao_richness_inverts_all_years <- as.data.frame(data_all_inverts %>% 
  dplyr::select(N))  %>% 
  ChaoRichness(datatype="abundance")
chao_richness_inverts_all_years
# Observed Estimator Est_s.e. 95% Lower 95% Upper
# N       58    77.589   14.352    63.421   128.791

# endpoint NULL = double the reference sample size.
estimate_richness_inverts <- iNEXT(input_all, q= c(0), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=500) # q=0 estimates richness

coverage_based_inverts <- ggiNEXT(estimate_richness_inverts, type=3)
coverage_based_inverts <- coverage_based_inverts + ggtitle("Macroeukaryotes")
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_1000_CB_inverts.png", plot = coverage_based_inverts, width=250, height=200, units="mm",dpi=300)

sample_size_based_inverts <- ggiNEXT(estimate_richness_inverts, type=1)
sample_size_based_inverts <- sample_size_based_inverts + ggtitle("Macroeukaryotes")
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/cov_based_raref_1000_SSB_inverts.png", plot = sample_size_based_inverts, width=250, height=200, units="mm",dpi=300)

############################
####### Final figure #######
############################

gamma_diversity_SSB <- ggarrange(sample_size_based_16S,
                             sample_size_based_18S,
                             sample_size_based_inverts,
                        labels = c("A", "B", "C"), ncol = 3, nrow = 1)
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/gamma_all_CB_normalized_SSB_curves.png", 
       plot = gamma_diversity_SSB , width=600, height=200, units="mm",dpi=300)

gamma_diversity_CB <- ggarrange(coverage_based_16S,
                                coverage_based_18S,
                             coverage_based_inverts,
                             labels = c("D", "E", "F"), ncol = 3, nrow = 1)
ggsave("R_Code_and_Analysis/gammadiversity/coverage_based_raref_tests/gamma_all_CB_normalized_CB_curves.png", 
       plot = gamma_diversity_CB , width=600, height=200, units="mm",dpi=300)
