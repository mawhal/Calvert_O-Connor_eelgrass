### NMDS plots tests macroeuk 18S ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 26, 2019 ###

library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library(tidyverse)


#####################################
############ Macroeuk18S ############
#####################################

##### 1000 CORRECTED FOR SINGLETONS

master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv", header=TRUE)

### Creating an object to store abundances only
abundances_macroeuk18S_NMDS <- master_df_Macro18S %>%
  dplyr::select(-(1:9))

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(log(abundances_macroeuk18S_NMDS+1), distance = "bray", k=2)
NMDS.18S.LOG

stressplot(NMDS.18S.LOG)
plot(NMDS.18S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.18S.LOG$points[,1]
MDS2 = NMDS.18S.LOG$points[,2]
NMDS_18S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = master_df_Macro18S$year, region = master_df_Macro18S$region, SampleID=master_df_Macro18S$SampleID)
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

nmds_Macro18S <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macro18S") +
  #annotate("text", label = "stress = 0.15", x = 1.3, y = -2.1, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17,18))

nmds_Macro18S <- nmds_Macro18S +  theme_bw() +
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

nmds_Macro18S

ggsave("R_Code_and_Analysis/betadiversity/tests_macroeuk18S/macroeuk18S_1000_COR_SING.png", plot = nmds_Macro18S, width = 15, height = 12, dpi=100)

# nmds_Macro18S <- nmds_Macro18S + ggrepel::geom_text_repel(data = NMDS_18S, aes(x=NMDS1, y=NMDS2, label = SampleID), cex = 3, direction = "both", segment.size = 0.25)
# nmds_Macro18S

##### 1000 NOT CORRECTED FOR SINGLETONS

master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_NOT_COR_SING_COVERAGE_RAREF.csv", header=TRUE)

### Creating an object to store abundances only
abundances_macroeuk18S_NMDS <- master_df_Macro18S %>%
  dplyr::select(-(1:9))

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(log(abundances_macroeuk18S_NMDS+1), distance = "bray", k=2)
NMDS.18S.LOG

stressplot(NMDS.18S.LOG)
plot(NMDS.18S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.18S.LOG$points[,1]
MDS2 = NMDS.18S.LOG$points[,2]
NMDS_18S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = master_df_Macro18S$year, region = master_df_Macro18S$region, SampleID=master_df_Macro18S$SampleID)
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

nmds_Macro18S <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macro18S") +
  #annotate("text", label = "stress = 0.15", x = 1.3, y = -2.1, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17,18))

nmds_Macro18S <- nmds_Macro18S +  theme_bw() +
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

nmds_Macro18S

ggsave("R_Code_and_Analysis/betadiversity/tests_macroeuk18S/macroeuk18S_1000_NOT_COR_SING.png", plot = nmds_Macro18S, width = 15, height = 12, dpi=100)


##### 2000 NOT CORRECTED FOR SINGLETONS

master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_2000_NOT_COR_SING_COVERAGE_RAREF.csv", header=TRUE)

### Creating an object to store abundances only
abundances_macroeuk18S_NMDS <- master_df_Macro18S %>%
  dplyr::select(-(1:9))

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(log(abundances_macroeuk18S_NMDS+1), distance = "bray", k=2)
NMDS.18S.LOG

stressplot(NMDS.18S.LOG)
plot(NMDS.18S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.18S.LOG$points[,1]
MDS2 = NMDS.18S.LOG$points[,2]
NMDS_18S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = master_df_Macro18S$year, region = master_df_Macro18S$region, SampleID=master_df_Macro18S$SampleID)
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

nmds_Macro18S <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macro18S") +
  #annotate("text", label = "stress = 0.15", x = 1.3, y = -2.1, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17,18))

nmds_Macro18S <- nmds_Macro18S +  theme_bw() +
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

nmds_Macro18S

ggsave("R_Code_and_Analysis/betadiversity/tests_macroeuk18S/macroeuk18S_2000_NOT_COR_SING.png", plot = nmds_Macro18S, width = 15, height = 12, dpi=100)

##### 5000 NOT CORRECTED FOR SINGLETONS

master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_5000_NOT_COR_SING_COVERAGE_RAREF.csv", header=TRUE)

### Creating an object to store abundances only
abundances_macroeuk18S_NMDS <- master_df_Macro18S %>%
  dplyr::select(-(1:9))

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(log(abundances_macroeuk18S_NMDS+1), distance = "bray", k=2)
NMDS.18S.LOG

stressplot(NMDS.18S.LOG)
plot(NMDS.18S.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.18S.LOG$points[,1]
MDS2 = NMDS.18S.LOG$points[,2]
NMDS_18S = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = master_df_Macro18S$year, region = master_df_Macro18S$region, SampleID=master_df_Macro18S$SampleID)
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

nmds_Macro18S <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macro18S") +
  #annotate("text", label = "stress = 0.15", x = 1.3, y = -2.1, size = 4, colour = "black") +
  scale_colour_manual(values=c("slateblue1", "sienna1", "yellow3", "#2a9958", "hotpink2")) +
  scale_shape_manual(values=c(19,8,17,18))

nmds_Macro18S <- nmds_Macro18S +  theme_bw() +
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

nmds_Macro18S

ggsave("R_Code_and_Analysis/betadiversity/tests_macroeuk18S/macroeuk18S_5000_NOT_COR_SING.png", plot = nmds_Macro18S, width = 15, height = 12, dpi=100)
