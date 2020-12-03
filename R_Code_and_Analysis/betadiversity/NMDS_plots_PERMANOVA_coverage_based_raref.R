### NMDS plots ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 26, 2019 ###

library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library(tidyverse)

#########################################
############ 16S prokaryotes ############
#########################################

### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_ASV)[1:16]

### Creating an object to store abundances only
abundances_16S_NMDS <- microbes_16S_ASV %>% 
  dplyr::select(-(1:15))

### Get MDS stats
set.seed(2)
NMDS.16S.LOG <- metaMDS(log(abundances_16S_NMDS+1), distance = "bray", k=2)  
NMDS.16S.LOG 

stressplot(NMDS.16S.LOG)
plot(NMDS.16S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.16S.LOG$points[,1]
MDS2 = NMDS.16S.LOG$points[,2]
NMDS_16S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = microbes_16S_ASV$year, region = microbes_16S_ASV$region)
NMDS_16S

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_16S, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_16S
NMDS_16S$region

# re-order the factor levels before the plot
NMDS_16S$region <- factor(NMDS_16S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_16S$year <- factor(NMDS_16S$year, levels=c("2015", "2016", "2017","2018"))

nmds_prokaryotes <- ggplot(NMDS_16S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Prokaryotes") + 
  annotate("text", label = "stress = 0.20", x = 1.1, y = -1.5, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17,18))

nmds_prokaryotes <- nmds_prokaryotes +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size

nmds_prokaryotes
ggsave("R_Code_and_Analysis/betadiversity/NMDS_prokaryotes_1000_COVER_BAS_RAR.png", plot = nmds_prokaryotes, width=250, height=200, units="mm",dpi=300)

### PERMANOVA 16S ###
#### LOG-transformation
log_16S <- log1p(abundances_16S_NMDS)
####  Distance matrix * this is using Bray-Curtis  ###
bray_16S <- vegdist(log_16S ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_16S <-adonis2(bray_16S ~ region + year,
                                         data=microbes_16S_ASV, permutations=999, by = "margin")
permanova_16S
write.csv(permanova_16S, "R_Code_and_Analysis/betadiversity/permanova_16S.csv")

#### run PERMDISP
as.factor(microbes_16S_ASV$region)
permdisp_16S_region <- betadisper(bray_16S, region, type = c("median","centroid")) 
plot(permdisp_16S_region)
boxplot(permdisp_16S_region)
permutest(permdisp_16S_region, pairwise = TRUE)

as.factor(microbes_16S_ASV$year)
permdisp_16S_year <- betadisper(bray_16S, year, type = c("median","centroid")) 
plot(permdisp_16S_year)
boxplot(permdisp_16S_year)
permutest(permdisp_16S_year, pairwise = TRUE)

#############################################
############ 18S microeukaryotes ############
#############################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_18S_ASV)[1:16]

### Creating an object to store abundances only
abundances_18S_NMDS <- microbes_18S_ASV %>% 
  dplyr::select(-(1:9))

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(log(abundances_18S_NMDS+1), distance = "bray", k=2)  
NMDS.18S.LOG 

stressplot(NMDS.18S.LOG)
plot(NMDS.18S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.18S.LOG$points[,1]
MDS2 = NMDS.18S.LOG$points[,2]
NMDS_18S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = microbes_18S_ASV$year, region = microbes_18S_ASV$region)
NMDS_18S

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_18S, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_18S
NMDS_18S$region

# re-order the factor levels before the plot
NMDS_18S$region <- factor(NMDS_18S$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_18S$year <- factor(NMDS_18S$year, levels=c("2015", "2016", "2017","2018"))

nmds_microeukaryotes <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Microeukaryotes") + 
  annotate("text", label = "stress = 0.15", x = 1.3, y = -2.1, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17,18))

nmds_microeukaryotes <- nmds_microeukaryotes +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size

nmds_microeukaryotes
ggsave("R_Code_and_Analysis/betadiversity/NMDS_microeukaryotes_1000_COVER_BAS_RAR.png", plot = nmds_microeukaryotes, width=250, height=200, units="mm",dpi=300)

### PERMANOVA 18S ###
#### LOG-transformation
log_18S <- log1p(abundances_18S_NMDS)
####  Distance matrix * this is using Bray-Curtis  ###
bray_18S <- vegdist(log_18S ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_18S <-adonis2(bray_18S ~ region + year,
                        data=microbes_18S_ASV, permutations=999, by = "margin")
permanova_18S
write.csv(permanova_18S, "R_Code_and_Analysis/betadiversity/permanova_18S.csv")

#### run PERMDISP
as.factor(microbes_18S_ASV$region)
permdisp_18S_region <- betadisper(bray_18S, microbes_18S_ASV$region, type = c("median","centroid")) 
plot(permdisp_18S_region)
boxplot(permdisp_18S_region)
permutest(permdisp_18S_region, pairwise = TRUE)

as.factor(microbes_18S_ASV$year)
permdisp_18S_year <- betadisper(bray_18S, microbes_18S_ASV$year, type = c("median","centroid")) 
plot(permdisp_18S_year)
boxplot(permdisp_18S_year)
permutest(permdisp_18S_year, pairwise = TRUE)

#################################################
############ Inverts Macroeukaryotes ############
#################################################

#### Finest level ####

### Read table metadata and abundances
inverts_finest <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")
names(inverts_finest)[1:12]

### Creating an object to store abundances only
abundances_inverts_finest <- inverts_finest %>% 
  dplyr::select(-(1:7))

### Get MDS stats
set.seed(2)
NMDS.inverts.LOG <- metaMDS(log(abundances_inverts_finest+1), distance = "bray", k=2)  
NMDS.inverts.LOG 

stressplot(NMDS.inverts.LOG)
plot(NMDS.inverts.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.inverts.LOG$points[,1]
MDS2 = NMDS.inverts.LOG$points[,2]
NMDS_inverts = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = inverts_finest$year, region = inverts_finest$region)
NMDS_inverts

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_inverts, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_inverts
NMDS_inverts$region

# re-order the factor levels before the plot
NMDS_inverts$region <- factor(NMDS_inverts$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_inverts$year <- factor(NMDS_inverts$year, levels=c("2014", "2015", "2016","2017"))

nmds_macroeukaryotes <- ggplot(NMDS_inverts, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.19", x = 1.0, y = -1.5, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(0,19,8,17))

nmds_macroeukaryotes <- nmds_macroeukaryotes +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size

nmds_macroeukaryotes
ggsave("R_Code_and_Analysis/betadiversity/NMDS_macroeukaryotes_finest_1000_COVER_BAS_RAR.png", plot = nmds_macroeukaryotes, width=250, height=200, units="mm",dpi=300)

### PERMANOVA inverts ###
#### LOG-transformation
log_inverts <- log1p(abundances_inverts_finest)
####  Distance matrix * this is using Bray-Curtis  ###
bray_inverts <- vegdist(log_inverts ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_inverts <-adonis2(bray_inverts ~ region + year,
                                      data=inverts_finest, permutations=999, by = "margin")
permanova_inverts
write.csv(permanova_inverts, "R_Code_and_Analysis/betadiversity/permanova_inverts.csv")


#### Family level ####

inverts_family <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_family_1000_COVERAGE_RAREF.csv")
names(inverts_family)[1:12]
inverts_family <- inverts_family %>% 
  drop_na()  #remove samples with NAs (why did they drop out with CB rarefaction?)

### Creating an object to store abundances only
abundances_inverts_family <- inverts_family %>% 
  dplyr::select(-(1:7))

### Get MDS stats
set.seed(2)
NMDS.inverts.LOG <- metaMDS(log(abundances_inverts_family+1), distance = "bray", k=2)  
NMDS.inverts.LOG 

stressplot(NMDS.inverts.LOG)
plot(NMDS.inverts.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.inverts.LOG$points[,1]
MDS2 = NMDS.inverts.LOG$points[,2]
NMDS_inverts = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = inverts_family$year, region = inverts_family$region)
NMDS_inverts

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_inverts, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_inverts
NMDS_inverts$region

# re-order the factor levels before the plot
NMDS_inverts$region <- factor(NMDS_inverts$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_inverts$year <- factor(NMDS_inverts$year, levels=c("2014", "2015", "2016","2017"))

nmds_macroeukaryotes_family <- ggplot(NMDS_inverts, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.19", x = 1.0, y = -1.5, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(0,19,8,17))

nmds_macroeukaryotes_family <- nmds_macroeukaryotes_family +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size

nmds_macroeukaryotes_family
ggsave("R_Code_and_Analysis/betadiversity/NMDS_macroeukaryotes_family_1000_COVER_BAS_RAR.png", plot = nmds_macroeukaryotes_family, width=250, height=200, units="mm",dpi=300)


####### without 2014 ########
### macroeukaryotes (inverts) without 2014
inverts_finest_no_2014 <- inverts_finest %>% 
  filter(!year == 2014)

### Creating an object to store abundances only
abundances_inverts_NMDS_no_2014 <- inverts_finest_no_2014 %>% 
  dplyr::select(-(1:7))

set.seed(2)
NMDS.inverts.LOG_no_2014 <- metaMDS(log(abundances_inverts_NMDS_no_2014+1), distance = "bray", k=2)
NMDS.inverts.LOG_no_2014 

stressplot(NMDS.inverts.LOG_no_2014)
plot(NMDS.inverts.LOG_no_2014)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.inverts.LOG_no_2014$points[,1]
MDS2 = NMDS.inverts.LOG_no_2014$points[,2]
NMDS_inverts_no_2014 = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = inverts_finest_no_2014$year, region = inverts_finest_no_2014$region)
NMDS_inverts_no_2014

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_inverts_no_2014, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_inverts_no_2014
NMDS_inverts_no_2014$region

# re-order the factor levels before the plot
NMDS_inverts_no_2014$region <- factor(NMDS_inverts_no_2014$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_inverts_no_2014$year <- factor(NMDS_inverts_no_2014$year, levels=c("2015", "2016","2017"))

nmds_macroeukaryotes_no_2014 <- ggplot(NMDS_inverts_no_2014, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.15", x = 1.0, y = -1.5, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17))

nmds_macroeukaryotes_no_2014 <- nmds_macroeukaryotes_no_2014 +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size

nmds_macroeukaryotes_no_2014
ggsave("R_Code_and_Analysis/betadiversity/NMDS_macroeukaryotes_no_2014.png", plot = nmds_macroeukaryotes_no_2014, width=250, height=200, units="mm",dpi=300)

### PERMANOVA inverts ###
#### LOG-transformation
log_inverts_no_2014 <- log1p(abundances_inverts_NMDS_no_2014)
####  Distance matrix * this is using Bray-Curtis  ###
bray_inverts_no_2014 <- vegdist(log_inverts_no_2014 ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_inverts_no_2014 <-adonis2(bray_inverts_no_2014 ~ region + year,
                        data=m.meta_no_2014, permutations=999, by = "margin")
permanova_inverts_no_2014
write.csv(permanova_inverts_no_2014, "R_Code_and_Analysis/betadiversity/permanova_inverts_no_2014.csv")

### TEST Permanova Inverts only 2017 and 2018
### macroeukaryotes (inverts) without 2014
inverts_finest_2016_2017 <- inverts_finest %>% 
  filter(!year == 2014, !year == 2015)

### Creating an object to store abundances only
abundances_inverts_NMDS_2016_2017 <- inverts_finest_2016_2017 %>% 
  dplyr::select(-(1:7))

### PERMANOVA inverts ###
#### LOG-transformation
log_inverts_2016_2017 <- log1p(abundances_inverts_NMDS_2016_2017)
####  Distance matrix * this is using Bray-Curtis  ###
bray_inverts_2016_2017 <- vegdist(log_inverts_2016_2017 ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_inverts_2016_2017 <-adonis2(bray_inverts_2016_2017 ~ region + year,
                                    data=m.meta_2016_2017, permutations=999, by = "margin")
permanova_inverts_2016_2017
write.csv(permanova_inverts_2016_2017, "R_Code_and_Analysis/betadiversity/permanova_inverts_2016_2017.csv")
### GRAPH 2016 and 2017
set.seed(2)
NMDS.inverts.LOG_2016_2017 <- metaMDS(log(abundances_inverts_NMDS_2016_2017+1), distance = "bray", k=2)
NMDS.inverts.LOG_2016_2017 

stressplot(NMDS.inverts.LOG_2016_2017)
plot(NMDS.inverts.LOG_2016_2017)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.inverts.LOG_2016_2017$points[,1]
MDS2 = NMDS.inverts.LOG_2016_2017$points[,2]
NMDS_inverts_2016_2017 = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = inverts_finest_2016_2017$year, region = inverts_finest_2016_2017$region)
NMDS_inverts_2016_2017

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS_inverts_2016_2017, old=c("MDS1", "MDS2", "year", "region"), new=c("NMDS1","NMDS2", "year", "region"))
NMDS_inverts_2016_2017
NMDS_inverts_2016_2017$region

# re-order the factor levels before the plot
NMDS_inverts_2016_2017$region <- factor(NMDS_inverts_2016_2017$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

NMDS_inverts_2016_2017$year <- factor(NMDS_inverts_2016_2017$year, levels=c("2016","2017"))

nmds_macroeukaryotes_2016_2017 <- ggplot(NMDS_inverts_2016_2017, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.16", x = 0.5, y = -1.5, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(8,17))

nmds_macroeukaryotes_2016_2017 <- nmds_macroeukaryotes_2016_2017 +  theme_bw() + 
  theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 16), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 16), #font size of legend
         plot.title = element_text(hjust = 0.1, size = 20, face = "bold")) #center plot title and set font size

nmds_macroeukaryotes_2016_2017
ggsave("R_Code_and_Analysis/betadiversity/NMDS_macroeukaryotes_2016_2017.png", plot = nmds_macroeukaryotes_2016_2017, width=250, height=200, units="mm",dpi=300)


plots_top <-cowplot::plot_grid( nmds_prokaryotes, nmds_microeukaryotes, nmds_macroeukaryotes, ncol=3, labels = c("A", "B", "C"), label_x =.05, hjust = 1, label_size=20)

plots_bottom <- cowplot::plot_grid(NULL,nmds_macroeukaryotes_no_2014, nmds_macroeukaryotes_2016_2017,NULL, ncol=4, nrow=1,labels = c("", "D", "E", ""), label_x =.05, hjust = 1, label_size=20, rel_widths=c(0.1,0.21,0.21, 0.1))

cowplot::plot_grid(plots_top, plots_bottom, ncol=1)

ggsave(paste0("R_Code_and_Analysis/betadiversity/NMDS_all.png"), width = 25, height = 12, dpi=100)
