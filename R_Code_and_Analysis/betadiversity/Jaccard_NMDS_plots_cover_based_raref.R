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

###SAME THING FOR PRESENCE-ABSENCE DATA = JACCARD
NMDS.16S_presabs <- abundances_16S_NMDS
NMDS.16S_presabs[NMDS.16S_presabs > 0] <- 1

### Get MDS stats
set.seed(2)
NMDS.16S.LOG <- metaMDS(NMDS.16S_presabs, distance = "jaccard", k=2)  
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
  ggtitle("Jaccard Prokaryotes") + 
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
ggsave("R_Code_and_Analysis/betadiversity/JACCARD_NMDS_prokaryotes_1000_COVER_BAS_RAR.png", plot = nmds_prokaryotes, width=250, height=200, units="mm",dpi=300)

#############################################
############ 18S microeukaryotes ############
#############################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_18S_ASV)[1:16]

### Creating an object to store abundances only
abundances_18S_NMDS <- microbes_18S_ASV %>% 
  dplyr::select(-(1:9))

###SAME THING FOR PRESENCE-ABSENCE DATA = JACCARD
NMDS.18S_presabs <- abundances_18S_NMDS
NMDS.18S_presabs[NMDS.18S_presabs > 0] <- 1

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(NMDS.18S_presabs, distance = "jaccard", k=2)  
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
  ggtitle("Jaccard Microeukaryotes") + 
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
ggsave("R_Code_and_Analysis/betadiversity/JACCARD_NMDS_microeukaryotes_1000_COVER_BAS_RAR.png", plot = nmds_microeukaryotes, width=250, height=200, units="mm",dpi=300)

#################################################
############ Inverts Macroeukaryotes ############
#################################################

#### Finest level ####

### Read table metadata and abundances
inverts_finest <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")
names(inverts_finest)[1:12]

### TEST Permanova Inverts only 2017 and 2018
### macroeukaryotes (inverts) without 2014
inverts_finest_2016_2017 <- inverts_finest %>% 
  filter(!year == 2014, !year == 2015)

### Creating an object to store abundances only
abundances_inverts_NMDS_2016_2017 <- inverts_finest_2016_2017 %>% 
  dplyr::select(-(1:7))

###SAME THING FOR PRESENCE-ABSENCE DATA = JACCARD
NMDS.inverts_presabs <- abundances_inverts_NMDS_2016_2017
NMDS.inverts_presabs[NMDS.inverts_presabs > 0] <- 1


### GRAPH 2016 and 2017
set.seed(2)
NMDS.inverts.LOG_2016_2017 <- metaMDS(NMDS.inverts_presabs, distance = "jaccard", k=2)
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
  ggtitle("Jaccard Macroeukaryotes") + 
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
ggsave("R_Code_and_Analysis/betadiversity/JACCARD_NMDS_macroeukaryotes_2016_2017.png", plot = nmds_macroeukaryotes_2016_2017, width=250, height=200, units="mm",dpi=300)


#####################################
############ Macroeuk18S ############
#####################################

### Read table metadata and abundances
master_df_Macro18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv", header=TRUE)

### Creating an object to store abundances only
abundances_macroeuk18S_NMDS <- master_df_Macro18S %>%
  dplyr::select(-(1:9))


###SAME THING FOR PRESENCE-ABSENCE DATA = JACCARD
NMDS.macroeuk18S_presabs <- abundances_macroeuk18S_NMDS
NMDS.macroeuk18S_presabs[NMDS.macroeuk18S_presabs > 0] <- 1

### Get MDS stats
set.seed(2)
NMDS.18S.LOG <- metaMDS(NMDS.macroeuk18S_presabs, distance = "jaccard", k=2)
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
  ggtitle("Jaccard Macroeuk18S") +
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

# 
# nmds_Macro18S <- nmds_Macro18S + ggrepel::geom_text_repel(data = NMDS_18S, aes(x=NMDS1, y=NMDS2, label = SampleID), cex = 3, direction = "both", segment.size = 0.25)
# nmds_Macro18S

ggsave("R_Code_and_Analysis/betadiversity/JACCARD_macroeuk18S_1000_COR_SING.png", plot = nmds_Macro18S, width=250, height=200, units="mm",dpi=100)

############################
####### Final figure #######
############################

plots_top <-cowplot::plot_grid( nmds_prokaryotes, nmds_microeukaryotes, ncol=2, labels = c("A", "B"), label_x =.05, hjust = 1, label_size=20)

plots_bottom <- cowplot::plot_grid( nmds_macroeukaryotes_2016_2017, nmds_Macro18S, ncol=2, labels = c("C", "D"), label_x =.05, hjust = 1, label_size=20)

plot_final <- cowplot::plot_grid(plots_top, plots_bottom, ncol=1)

ggsave("R_Code_and_Analysis/betadiversity/JACCARD_NMDS_all.png", plot = plot_final, width = 15, height = 12, dpi=100)

