### Alpha diversity (Chao1 richness, Shannon and Pielou's eveness)  ###
### Author: Bianca Trevizan Segovia ###
### Date created: May 19th, 2020 ###

library(vegan)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readr)
library(tidyverse)
library(RColorBrewer)

#########################################
############ 16S prokaryotes ############
#########################################

### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", header=T)
names(microbes_16S_ASV)[1:16]

### Creating an object to store abundances only
abundances_16S <- microbes_16S_ASV %>% 
  dplyr::select(-(1:15))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_16S, index = "shannon")
Richness <- specnumber(abundances_16S) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_16S)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_16S <- data.frame(shannon, pielou,chao1, microbes_16S_ASV$region, microbes_16S_ASV$year)
names(alpha_16S)
### renaming columns new name = old name
alpha_16S_metrics <- alpha_16S %>%
  dplyr::rename(region =microbes_16S_ASV.region, year = microbes_16S_ASV.year) 

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.16S = reshape2::melt(alpha_16S_metrics, id.var=c("year", "region"))
dat.16S$year = factor(dat.16S$year, levels=c("2015","2016", "2017", "2018"))
dat.16S$region <- factor(dat.16S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

#### year ####
p_16S_year <- ggplot(dat.16S, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_brewer(palette = "Dark2")

p_16S_year <- p_16S_year + labs(y = "alpha diversity measures prokaryotes") # y-axis label

p_16S_year <- p_16S_year + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 20),
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

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_16S_years.png", plot = p_16S_year, width=250, height=200, units="mm",dpi=300)

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
         legend.text= element_text(size = 20),
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

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_16S_regions.png", plot = p_16S_region, width=250, height=200, units="mm",dpi=300)

#############################################
############ 18S microeukaryotes ############
#############################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", header=T)
names(microbes_18S_ASV)[1:18]

### Creating an object to store abundances only
abundances_18S <- microbes_18S_ASV %>% 
  dplyr::select(-(1:15))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_18S, index = "shannon")
Richness <- specnumber(abundances_18S) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_18S)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_18S <- data.frame(shannon, pielou,chao1, microbes_18S_ASV$region, microbes_18S_ASV$year)
names(alpha_18S)
### renaming columns new name = old name
alpha_18S_metrics <- alpha_18S %>%
  dplyr::rename(region =microbes_18S_ASV.region, year =microbes_18S_ASV.year) 

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.18S = reshape2::melt(alpha_18S_metrics, id.var=c("year", "region"))
dat.18S$year = factor(dat.18S$year, levels=c("2015","2016", "2017", "2018"))
dat.18S$region <- factor(dat.18S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
levels(dat.18S$year)
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
         legend.text= element_text(size = 20),
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

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_18S_years.png", plot = p_18S_year, width=250, height=200, units="mm",dpi=300)

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
         legend.text= element_text(size = 20),
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

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_18S_regions.png", plot = p_18S_region, width=250, height=200, units="mm",dpi=300)

#################################################
############ Inverts macroeukaryotes ############
#################################################

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
  select( year, group, region, site, sample, ID, taxon = taxon2, remove, size ) %>% 
  filter( is.na(remove), !is.na(taxon))

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
  ungroup( year, group, region, site, taxon )

### Creating an object to store abundances only
abundances_inverts_finest <- m.meta_finest %>% 
  dplyr::select(-(1:6))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_inverts_finest, index = "shannon")
Richness <- specnumber(abundances_inverts_finest) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_inverts_finest)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_inverts_finest <- data.frame(shannon, pielou,chao1, m.meta_finest$region, m.meta_finest$year)
names(alpha_inverts_finest)
### renaming columns new name = old name
alpha_inverts_metrics_finest <- alpha_inverts_finest %>%
  dplyr::rename(region= m.meta_finest.region, year = m.meta_finest.year) 

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.inverts_finest = reshape2::melt(alpha_inverts_metrics_finest, id.var=c("year", "region"))
dat.inverts_finest$year = factor(dat.inverts_finest$year, levels=c("2014","2015", "2016", "2017"))
dat.inverts_finest$region <- factor(dat.inverts_finest$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

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
         legend.text= element_text(size = 20),
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

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_inverts_finest_years.png", plot = p_inverts_finest_year, width=250, height=200, units="mm",dpi=300)

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
         legend.text= element_text(size = 20),
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

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_inverts_finest_regions.png", plot = p_inverts_finest_region, width=250, height=200, units="mm",dpi=300)

#### Family level ####

# filter taxa and sites # taxon 4 = family
mfilt_family <- m %>%
  select( year, group, region, site, sample, ID, taxon = taxon4, remove, size ) %>% 
  filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum_family <- mfilt_family %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean_family <- m.sum_family %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance)) 

# make a community dataset
m.meta_family <- m.sum_family %>% 
  spread( taxon, abundance, fill=0 ) 

names(m.meta_family)[1:16]

m.meta_family <- m.meta_family %>% 
  ungroup( year, group, region, site, taxon )

### Creating an object to store abundances only
abundances_inverts_family <- m.meta_family %>% 
  dplyr::select(-(1:6))

# Calculate alpha diversity metrics
shannon <- diversity(abundances_inverts_family, index = "shannon")
Richness <- specnumber(abundances_inverts_family) 
pielou <- shannon/log(Richness)
chao1 <- estimateR(abundances_inverts_family)[2,] ### Chao1 (Estimated "real" Richness)

## creating data frame with alpha metrics and metadata
alpha_inverts_family <- data.frame(shannon, pielou,chao1, m.meta_family$region, m.meta_family$year)
names(alpha_inverts_family)
### renaming columns new name = old name
alpha_inverts_metrics_family <- alpha_inverts_family %>%
  dplyr::rename(region= m.meta_family.region, year = m.meta_family.year) 

###CREATE BOX_PLOT

# melting alpha metrics into one column called variable (containing all three measure names) and another column called value containing all values for those measures 
dat.inverts_family = reshape2::melt(alpha_inverts_metrics_family, id.var=c("year", "region"))
dat.inverts_family$year = factor(dat.inverts_family$year, levels=c("2014","2015", "2016", "2017"))
dat.inverts_family$region <- factor(dat.inverts_family$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

# set same years as same colours as microbes
#brewer.pal(n = 8, name = "Dark2") "#1B9E77" "#D95F02" "#7570B3"

#### year ####
p_inverts_family_year <- ggplot(dat.inverts_family, aes(x = year, y = value, fill = year)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_manual(values=c("#E6AB02", "#1B9E77", "#D95F02", "#7570B3"))

p_inverts_family_year <- p_inverts_family_year + labs(y = "alpha diversity measures macroeukaryotes family") # y-axis label

p_inverts_family_year <- p_inverts_family_year + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 20),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),   
         panel.border = element_blank()) #remove lines outside the graph

p_inverts_family_year

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_inverts_family_years.png", plot = p_inverts_family_year, width=250, height=200, units="mm",dpi=300)

#### region ####
p_inverts_family_region <- ggplot(dat.inverts_family, aes(x = region, y = value, fill = region)) + #fill allows to set different colors
  facet_wrap(. ~ variable, scale="free") +   
  geom_boxplot(outlier.shape = NA, outlier.alpha = 0.4, notch = F) + geom_jitter(alpha = 0.3, width = 0.2) +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) 

p_inverts_family_region <- p_inverts_family_region + labs(y = "alpha diversity measures macroeukaryotes family") # y-axis label

p_inverts_family_region <- p_inverts_family_region + ggtitle("") + #add title (optional)
  theme_bw() +
  theme (legend.position="bottom",
         legend.title = element_blank(),
         legend.spacing.x = unit(0.2, 'cm'),
         legend.text= element_text(size = 20),
         axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
         axis.title.x = element_blank(),
         axis.text.x = element_blank(), #font size of numbers in axis
         axis.text.y = element_text(size = 18),
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         strip.text.x = element_text(size = 18),   
         panel.border = element_blank()) #remove lines outside the graph

p_inverts_family_region

ggsave("R_Code_and_Analysis/alphadiversity/alpha_diversity_inverts_family_regions.png", plot = p_inverts_family_region, width=250, height=200, units="mm",dpi=300)

# arrange the three plots in a single row
p_16S_year_no_leg <- p_16S_year + theme(legend.position="none", plot.margin=unit(c(t=0,r=1,b=0,l=1),"cm"))
p_18S_year_no_leg <-p_18S_year + theme(legend.position="none", plot.margin=unit(c(t=0,r=1,b=0,l=0),"cm"))
p_inverts_finest_year_no_leg <- p_inverts_finest_year + theme(legend.position="none", plot.margin=unit(c(t=0,r=1,b=0,l=0),"cm"))
plots_year <- cowplot::plot_grid( p_16S_year_no_leg, p_18S_year_no_leg, p_inverts_finest_year_no_leg, ncol=3)
# extract the legend from one of the plots
legend <- get_legend(
  p_16S_year + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
year <- plot_grid(plots_year, legend, ncol = 1, rel_heights = c(1, .1))
year

# arrange the three plots in a single row
p_16S_region_no_leg <- p_16S_region + theme(legend.position="none", plot.margin=unit(c(t=0,r=1,b=0,l=1),"cm"))
p_18S_region_no_leg <-p_18S_region + theme(legend.position="none", plot.margin=unit(c(t=0,r=1,b=0,l=0),"cm"))
p_inverts_finest_region_no_leg <- p_inverts_finest_region + theme(legend.position="none", plot.margin=unit(c(t=0,r=1,b=0,l=0),"cm"))
plots_region <- cowplot::plot_grid( p_16S_region_no_leg, p_18S_region_no_leg, p_inverts_finest_region_no_leg, ncol=3)
# extract the legend from one of the plots
legend <- get_legend(
  p_16S_region + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
region <- plot_grid(plots_region, legend, ncol = 1, rel_heights = c(1, .1))
region

year_region <- plot_grid(year, region, ncol=1)
ggsave(paste0("R_Code_and_Analysis/alphadiversity/alpha_diversity_all.png"), width = 18, height = 15  )
