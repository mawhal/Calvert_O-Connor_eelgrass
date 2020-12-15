### Alpha diversity using coverage based rarefaction (Chao & Jost, 2012) corrected for singletons (Chiu & Chao 2016) ###
### Author: Bianca Trevizan Segovia ###
### Date created: September 1st, 2020 ###

##### ALPHA DIVERSITY #####
library(vegan)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#########################################
############ 16S prokaryotes ############
#########################################

#########################################
### WITH COVERAGE-BASED RAREFIED DATA ###
#########################################
# 
# ### Read table metadata and abundances
# microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# names(microbes_16S_ASV)[1:17]
# 
# ### Creating an object to store abundances only
# abundances_16S <- microbes_16S_ASV %>% 
#   dplyr::select(-(1:15))
# 
# # Calculate alpha diversity metrics
# shannon <- diversity(abundances_16S, index = "shannon")
# Richness <- specnumber(abundances_16S) 
# pielou <- shannon/log(Richness)
# chao1 <- estimateR(abundances_16S)[2,] ### Chao1 (Estimated "real" Richness)
# 
# ## creating data frame with alpha metrics and metadata
# alpha_16S <- data.frame(shannon, pielou,chao1, microbes_16S_ASV$region, microbes_16S_ASV$year)
# names(alpha_16S)
# ### renaming columns new name = old name
# alpha_16S_metrics <- alpha_16S %>%
#   dplyr::rename(region =microbes_16S_ASV.region, year = microbes_16S_ASV.year) 
# 
# # remove pielou for a more compact graph/result
# alpha_16S_metrics <- alpha_16S_metrics %>% 
#   select(-c(pielou, shannon))
# 
# ###CREATE BOX_PLOT
# 
# # melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
# dat.16S = reshape2::melt(alpha_16S_metrics, id.var=c("year", "region"))
# dat.16S$year = factor(dat.16S$year, levels=c("2015","2016", "2017", "2018"))
# dat.16S$region <- factor(dat.16S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# 
# #### year ####
# p_16S_year <- ggplot(dat.16S, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
#   facet_wrap(. ~ variable, scale="free") +   
#   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
#   scale_fill_brewer(palette = "Dark2")
# 
# p_16S_year <- p_16S_year + labs(y = "alpha diversity measures prokaryotes") # y-axis label
# 
# p_16S_year <- p_16S_year + ggtitle("Prokaryotes") + #add title (optional)
#   theme_bw() +
#   theme (legend.position="bottom",
#          legend.title = element_blank(),
#          legend.spacing.x = unit(0.2, 'cm'),
#          legend.text= element_text(size = 16),
#          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#          axis.title.x = element_blank(),
#          axis.text.x = element_blank(), #font size of numbers in axis
#          axis.text.y = element_text(size = 18),
#          panel.grid.major = element_blank(), #remove major grid
#          panel.grid.minor = element_blank(), #remove minor grid
#          axis.line = element_line(colour = "black"), #draw line in the axis
#          strip.text.x = element_text(size = 18),   
#          panel.border = element_blank()) #remove lines outside the graph
# 
# p_16S_year
# 
# ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_16S_1000_COVER_BAS_RAR_years.png", plot = p_16S_year, width=250, height=200, units="mm",dpi=300)
# 
# #### region ####
# p_16S_region <- ggplot(dat.16S, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
#   facet_wrap(. ~ variable, scale="free") +   
#   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
#   scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 
# 
# p_16S_region <- p_16S_region + labs(y = "alpha diversity measures prokaryotes") # y-axis label
# 
# p_16S_region <- p_16S_region + ggtitle("") + #add title (optional)
#   theme_bw() +
#   theme (legend.position="bottom",
#          legend.title = element_blank(),
#          legend.spacing.x = unit(0.2, 'cm'),
#          legend.text= element_text(size = 16),
#          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#          axis.title.x = element_blank(),
#          axis.text.x = element_blank(), #font size of numbers in axis
#          axis.text.y = element_text(size = 18),
#          panel.grid.major = element_blank(), #remove major grid
#          panel.grid.minor = element_blank(), #remove minor grid
#          axis.line = element_line(colour = "black"), #draw line in the axis
#          strip.text.x = element_text(size = 18),   
#          panel.border = element_blank()) #remove lines outside the graph
# 
# p_16S_region
# 
# ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_16S_1000_COVER_BAS_RAR_regions.png", plot = p_16S_region, width=250, height=200, units="mm",dpi=300)

#########################################
############ 16S prokaryotes ############
#########################################

#########################################
### WITH COVERAGE-BASED RAREFIED DATA ###
#########################################

### Read table metadata and abundances
microbes_16S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_genus)[1:17]

### Creating an object to store abundances only
abundances_16S <- microbes_16S_genus %>% 
  dplyr::select(-(1:12))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_16S, index = "shannon")
Richness <- specnumber(abundances_16S) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_16S)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_16S <- data.frame(shannon, pielou,chao1, microbes_16S_genus$region, microbes_16S_genus$year)
names(alpha_16S)
### renaming columns new name = old name
alpha_16S_metrics <- alpha_16S %>%
  dplyr::rename(region =microbes_16S_genus.region, year = microbes_16S_genus.year) 

# remove pielou for a more compact graph/result
alpha_16S_metrics <- alpha_16S_metrics %>% 
  select(-c(pielou, shannon))

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.16S = reshape2::melt(alpha_16S_metrics, id.var=c("year", "region"))
dat.16S$year = factor(dat.16S$year, levels=c("2015","2016", "2017", "2018"))
dat.16S$region <- factor(dat.16S$region, levels=c("pruth", "choked", "triquet","goose","mcmullins"))

#### region ####
p_16S_region <- ggplot(dat.16S, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 

p_16S_region <- p_16S_region + labs(y = "alpha diversity measures prokaryotes") # y-axis label

p_16S_region <- p_16S_region + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 16),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),   
         panel.border = element_blank()) #remove lines outside the graph

p_16S_region

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_16S_1000_GENUS_COVER_BAS_RAR_regions.png", plot = p_16S_region, width=250, height=200, units="mm",dpi=300)

### year ####
p_16S_year <- ggplot(dat.16S, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_brewer(palette = "Dark2")

p_16S_year <- p_16S_year + labs(y = "alpha diversity measures prokaryotes") # y-axis label

p_16S_year <- p_16S_year + ggtitle("Prokaryotes") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 16),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),
         panel.border = element_blank()) #remove lines outside the graph

p_16S_year

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_16S_1000_GENUS_COVER_BAS_RAR_years.png", plot = p_16S_year, width=250, height=200, units="mm",dpi=300)

#########################################
############ 18S prokaryotes ############
# #########################################
# 
# #########################################
# ### WITH COVERAGE-BASED RAREFIED DATA ###
# #########################################
# 
# ### Read table metadata and abundances
# microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
# 
# names(microbes_18S_ASV)[1:17]
# 
# ### Creating an object to store abundances only
# abundances_18S <- microbes_18S_ASV %>% 
#   dplyr::select(-(1:9))
# 
# # Calculate alpha diversity metrics
# shannon <- diversity(abundances_18S, index = "shannon")
# Richness <- specnumber(abundances_18S) 
# pielou <- shannon/log(Richness)
# chao1 <- estimateR(abundances_18S)[2,] ### Chao1 (Estimated "real" Richness)
# 
# ## creating data frame with alpha metrics and metadata
# alpha_18S <- data.frame(shannon, pielou,chao1, microbes_18S_ASV$region, microbes_18S_ASV$year)
# names(alpha_18S)
# ### renaming columns new name = old name
# alpha_18S_metrics <- alpha_18S %>%
#   dplyr::rename(region =microbes_18S_ASV.region, year = microbes_18S_ASV.year) 
# 
# # remove pielou for a more compact graph/result
# alpha_18S_metrics <- alpha_18S_metrics %>% 
#   select(-c(pielou, shannon))
# 
# ###CREATE BOX_PLOT
# 
# # melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
# dat.18S = reshape2::melt(alpha_18S_metrics, id.var=c("year", "region"))
# dat.18S$year = factor(dat.18S$year, levels=c("2015","2016", "2017", "2018"))
# dat.18S$region <- factor(dat.18S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# 
# 
# #### year ####
# p_18S_year <- ggplot(dat.18S, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
#   facet_wrap(. ~ variable, scale="free") +   
#   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
#   scale_fill_brewer(palette = "Dark2")
# 
# p_18S_year <- p_18S_year + labs(y = "alpha diversity measures microeukaryotes") # y-axis label
# 
# p_18S_year <- p_18S_year + ggtitle("") + #add title (optional)
#   theme_bw() +
#   theme (legend.position="bottom",
#          legend.title = element_blank(),
#          legend.spacing.x = unit(0.2, 'cm'),
#          legend.text= element_text(size = 16),
#          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#          axis.title.x = element_blank(),
#          axis.text.x = element_blank(), #font size of numbers in axis
#          axis.text.y = element_text(size = 18),
#          panel.grid.major = element_blank(), #remove major grid
#          panel.grid.minor = element_blank(), #remove minor grid
#          axis.line = element_line(colour = "black"), #draw line in the axis
#          strip.text.x = element_text(size = 18),   
#          panel.border = element_blank()) #remove lines outside the graph
# 
# p_18S_year
# 
# ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_18S_1000_COVER_BAS_RAR_years.png", plot = p_18S_year, width=250, height=200, units="mm",dpi=300)
# 
# 
# #### region ####
# p_18S_region <- ggplot(dat.18S, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
#   facet_wrap(. ~ variable, scale="free") +   
#   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
#   scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 
# 
# p_18S_region <- p_18S_region + labs(y = "alpha diversity measures microeukaryotes") # y-axis label
# 
# p_18S_region <- p_18S_region + ggtitle("") + #add title (optional)
#   theme_bw() +
#   theme (legend.position="bottom",
#          legend.title = element_blank(),
#          legend.spacing.x = unit(0.2, 'cm'),
#          legend.text= element_text(size = 16),
#          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#          axis.title.x = element_blank(),
#          axis.text.x = element_blank(), #font size of numbers in axis
#          axis.text.y = element_text(size = 18),
#          panel.grid.major = element_blank(), #remove major grid
#          panel.grid.minor = element_blank(), #remove minor grid
#          axis.line = element_line(colour = "black"), #draw line in the axis
#          strip.text.x = element_text(size = 18),   
#          panel.border = element_blank()) #remove lines outside the graph
# 
# p_18S_region
# 
# ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_18S_1000_COVER_BAS_RAR_regions.png", plot = p_18S_region, width=250, height=200, units="mm",dpi=300)


#########################################
############ 18S prokaryotes ############
#########################################

#########################################
### WITH COVERAGE-BASED RAREFIED DATA ###
#########################################

### Read table metadata and abundances
microbes_18S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv", header=T)

names(microbes_18S_genus)[1:17]

### Creating an object to store abundances only
abundances_18S <- microbes_18S_genus %>% 
  dplyr::select(-(1:9))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_18S, index = "shannon")
Richness <- specnumber(abundances_18S) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_18S)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_18S <- data.frame(shannon, pielou,chao1, microbes_18S_genus$region, microbes_18S_genus$year)
names(alpha_18S)
### renaming columns new name = old name
alpha_18S_metrics <- alpha_18S %>%
  dplyr::rename(region =microbes_18S_genus.region, year = microbes_18S_genus.year) 

# remove pielou for a more compact graph/result
alpha_18S_metrics <- alpha_18S_metrics %>% 
  select(-c(pielou, shannon))

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.18S = reshape2::melt(alpha_18S_metrics, id.var=c("year", "region"))
dat.18S$year = factor(dat.18S$year, levels=c("2015","2016", "2017", "2018"))
dat.18S$region <- factor(dat.18S$region, levels=c("pruth", "choked", "triquet","goose","mcmullins"))

#### region ####
p_18S_region <- ggplot(dat.18S, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 

p_18S_region <- p_18S_region + labs(y = "alpha diversity measures microeukaryotes") # y-axis label

p_18S_region <- p_18S_region + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 16),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),   
         panel.border = element_blank()) #remove lines outside the graph

p_18S_region

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_18S_1000_GENUS_COVER_BAS_RAR_regions.png", plot = p_18S_region, width=250, height=200, units="mm",dpi=300)

#### year ####
p_18S_year <- ggplot(dat.18S, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_brewer(palette = "Dark2")

p_18S_year <- p_18S_year + labs(y = "alpha diversity measures microeukaryotes") # y-axis label

p_18S_year <- p_18S_year + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 16),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),
         panel.border = element_blank()) #remove lines outside the graph

p_18S_year

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_18S_1000_GENUS_COVER_BAS_RAR_years.png", plot = p_18S_year, width=250, height=200, units="mm",dpi=300)

#################################################
############ Inverts macroeukaryotes ############
#################################################

#### Finest level ####

### Read table metadata and abundances
inverts_finest <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")
names(inverts_finest)[1:12]
### Creating an object to store abundances only
abundances_inverts_finest <- inverts_finest %>% 
  dplyr::select(-(1:7))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_inverts_finest, index = "shannon")
Richness <- specnumber(abundances_inverts_finest) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_inverts_finest)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_inverts_finest <- data.frame(shannon, pielou,chao1, inverts_finest$region, inverts_finest$year)
names(alpha_inverts_finest)
### renaming columns new name = old name
alpha_inverts_metrics_finest <- alpha_inverts_finest %>%
  dplyr::rename(region=inverts_finest.region, year =inverts_finest.year) 

# remove pielou for a more compact graph/result
alpha_inverts_metrics_finest <- alpha_inverts_metrics_finest %>% 
  select(-c(pielou, shannon))

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.inverts_finest = reshape2::melt(alpha_inverts_metrics_finest, id.var=c("year", "region"))
dat.inverts_finest$year = factor(dat.inverts_finest$year, levels=c("2014","2015", "2016", "2017"))
dat.inverts_finest$region <- factor(dat.inverts_finest$region, levels=c("pruth", "choked", "triquet","goose","mcmullins"))

# set same years as same colours as microbes
#brewer.pal(n = 8, name = "Dark2") "#1B9E77" "#D95F02" "#7570B3"

#### year ####
p_inverts_finest_year <- ggplot(dat.inverts_finest, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_manual(values=c("#E6AB02", "#1B9E77", "#D95F02", "#7570B3"))

p_inverts_finest_year <- p_inverts_finest_year + labs(y = "alpha diversity measures macroeukaryotes finest") # y-axis label

p_inverts_finest_year <- p_inverts_finest_year + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 16),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),   
         panel.border = element_blank()) #remove lines outside the graph

p_inverts_finest_year

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_inverts_finest_1000_COVER_BAS_RAR_years.png", plot = p_inverts_finest_year, width=250, height=200, units="mm",dpi=300)

#### region ####
p_inverts_finest_region <- ggplot(dat.inverts_finest, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 

p_inverts_finest_region <- p_inverts_finest_region + labs(y = "alpha diversity measures macroeukaryotes finest") # y-axis label

p_inverts_finest_region <- p_inverts_finest_region + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 16),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),   
         panel.border = element_blank()) #remove lines outside the graph

p_inverts_finest_region

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_inverts_finest_1000_COVER_BAS_RAR_regions.png", plot = p_inverts_finest_region, width=250, height=200, units="mm",dpi=300)

#####################################
############ Macroeuk18S ############
#####################################

#########################################
### WITH COVERAGE-BASED RAREFIED DATA ###
#########################################

# ### Read table metadata and abundances
# macroeuk18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv", header=TRUE)
# names(macroeuk18S_ASV)[1:17]
# 
# ### Creating an object to store abundances only
# abundances_macroeuk18S <- macroeuk18S_ASV %>% 
#   dplyr::select(-(1:9))
# 
# # Calculate alpha diversity metrics
# shannon <- diversity(abundances_macroeuk18S, index = "shannon")
# Richness <- specnumber(abundances_macroeuk18S) 
# pielou <- shannon/log(Richness)
# chao1 <- estimateR(abundances_macroeuk18S)[2,] ### Chao1 (Estimated "real" Richness)
# 
# ## creating data frame with alpha metrics and metadata
# alpha_macroeuk18S <- data.frame(shannon, pielou,chao1, macroeuk18S_ASV$region, macroeuk18S_ASV$year)
# names(alpha_macroeuk18S)
# ### renaming columns new name = old name
# alpha_macroeuk18S_metrics <- alpha_macroeuk18S %>%
#   dplyr::rename(region = macroeuk18S_ASV.region, year = macroeuk18S_ASV.year) 
# 
# # remove pielou for a more compact graph/result
# alpha_macroeuk18S_metrics <- alpha_macroeuk18S_metrics %>% 
#   select(-c(pielou, shannon))
# 
# ###CREATE BOX_PLOT
# 
# # melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
# dat.macroeuk18S = reshape2::melt(alpha_macroeuk18S_metrics, id.var=c("year", "region"))
# dat.macroeuk18S$year = factor(dat.macroeuk18S$year, levels=c("2015","2016", "2017", "2018"))
# dat.macroeuk18S$region <- factor(dat.macroeuk18S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# 
# 
# #### year ####
# p_macroeuk18S_year <- ggplot(dat.macroeuk18S, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
#   facet_wrap(. ~ variable, scale="free") +   
#   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
#   scale_fill_brewer(palette = "Dark2")
# 
# p_macroeuk18S_year <- p_macroeuk18S_year + labs(y = "alpha diversity measures microeukaryotes") # y-axis label
# 
# p_macroeuk18S_year <- p_macroeuk18S_year + ggtitle("") + #add title (optional)
#   theme_bw() +
#   theme (legend.position="bottom",
#          legend.title = element_blank(),
#          legend.spacing.x = unit(0.2, 'cm'),
#          legend.text= element_text(size = 16),
#          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#          axis.title.x = element_blank(),
#          axis.text.x = element_blank(), #font size of numbers in axis
#          axis.text.y = element_text(size = 18),
#          panel.grid.major = element_blank(), #remove major grid
#          panel.grid.minor = element_blank(), #remove minor grid
#          axis.line = element_line(colour = "black"), #draw line in the axis
#          strip.text.x = element_text(size = 18),   
#          panel.border = element_blank()) #remove lines outside the graph
# 
# p_macroeuk18S_year
# 
# ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_macroeuk18S_1000_COVER_BAS_RAR_years.png", plot = p_macroeuk18S_year, width=250, height=200, units="mm",dpi=300)
# 
# 
# #### region ####
# p_macroeuk18S_region <- ggplot(dat.macroeuk18S, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
#   facet_wrap(. ~ variable, scale="free") +   
#   geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
#   scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 
# 
# p_macroeuk18S_region <- p_macroeuk18S_region + labs(y = "alpha diversity measures microeukaryotes") # y-axis label
# 
# p_macroeuk18S_region <- p_macroeuk18S_region + ggtitle("") + #add title (optional)
#   theme_bw() +
#   theme (legend.position="bottom",
#          legend.title = element_blank(),
#          legend.spacing.x = unit(0.2, 'cm'),
#          legend.text= element_text(size = 16),
#          axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#          axis.title.x = element_blank(),
#          axis.text.x = element_blank(), #font size of numbers in axis
#          axis.text.y = element_text(size = 18),
#          panel.grid.major = element_blank(), #remove major grid
#          panel.grid.minor = element_blank(), #remove minor grid
#          axis.line = element_line(colour = "black"), #draw line in the axis
#          strip.text.x = element_text(size = 18),   
#          panel.border = element_blank()) #remove lines outside the graph
# 
# p_macroeuk18S_region
# 
# ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_macroeuk18S_1000_COVER_BAS_RAR_regions.png", plot = p_macroeuk18S_region, width=250, height=200, units="mm",dpi=300)


############################
####### Final figure #######
############################

# arrange the three plots in a single row
p_16S_year_no_leg <- p_16S_year + theme(legend.position="bottom", 
                                        plot.margin=unit(c(t=0.3,r=1,b=0.5,l=0.3),"cm"), 
                                        plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + 
  labs(y = "Alpha diversity measures") + ggtitle("Prokaryotes")
p_18S_year_no_leg <- p_18S_year + theme(legend.position="bottom", plot.margin=unit(c(t=0.3,r=1,b=0.5,l=1),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "") + ggtitle("Microeukaryotes")
p_inverts_finest_year_no_leg <- p_inverts_finest_year + theme(legend.position="bottom", plot.margin=unit(c(t=0.3,r=1,b=0.5,l=1),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "") + ggtitle("Macroeukaryotes")
# p_macroeuk18S_year_no_leg <- p_macroeuk18S_year + theme(legend.position="bottom", plot.margin=unit(c(t=0.3,r=1,b=0.5,l=1),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "") + ggtitle("Macroeuk18S")

plots_year <- cowplot::plot_grid( p_16S_year_no_leg, p_18S_year_no_leg, p_inverts_finest_year_no_leg, ncol=3,  labels = c("A", "B", "C"), label_x =.05, hjust = 1, label_size=20)
plots_year 
# extract the legend from one of the plots
# legend <- get_legend(
#   p_16S_year + 
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom"))
# # add the legend to the row we made earlier. Give it one-third of 
# # the width of one plot (via rel_widths).
year <- plot_grid(plots_year, ncol = 1, rel_heights = c(1, .1))
year

# arrange the three plots in a single row
p_16S_region_no_leg <- p_16S_region + theme(legend.position="bottom", plot.margin=unit(c(t=0,r=1,b=0,l=0.3),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "Alpha diversity measures") + ggtitle("Prokaryotes")
p_18S_region_no_leg <-p_18S_region + theme(legend.position="bottom", plot.margin=unit(c(t=0,r=1,b=0,l=0.3),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "") + ggtitle("Microeukaryotes")
p_inverts_finest_region_no_leg <- p_inverts_finest_region + theme(legend.position="bottom", plot.margin=unit(c(t=0,r=1,b=0,l=0.3),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "") + ggtitle("Macroeukaryotes")
# p_macroeuk18S_region_no_leg <- p_macroeuk18S_region + theme(legend.position="bottom", plot.margin=unit(c(t=0,r=1,b=0,l=0.3),"cm"), plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0))) + labs(y = "") + ggtitle("Macroeuk18S")

plots_region <- cowplot::plot_grid( p_16S_region_no_leg, p_18S_region_no_leg, p_inverts_finest_region_no_leg, ncol=3,  labels = c("D", "E", "F"), label_x =.05, hjust = 1, label_size=20)
plots_region
# # extract the legend from one of the plots
# legend <- get_legend(
#   p_16S_region + 
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom"))
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
region <- plot_grid(plots_region, ncol = 1, rel_heights = c(1, .1))
region

year_region <- plot_grid(year, region, nrow=2)
year_region
ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_all_GENUS_microbes_FINEST_inverst_COVER_BAS_RAR.png", plot = year_region, width = 470, height = 270, units="mm", dpi=300)
