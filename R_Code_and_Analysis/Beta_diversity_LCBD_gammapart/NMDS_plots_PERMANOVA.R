### NMDS plots ###
### Author: Bianca Trevizan Segovia ###
### Date created: November 26, 2019 ###

library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)

#########################################
############ 16S prokaryotes ############
#########################################

### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", header=T)
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
NMDS_16S$region <- factor(NMDS_16S$region, levels=c("choked", "pruth", "triquet","goose","mcmullin"))

NMDS_16S$year <- factor(NMDS_16S$year, levels=c("2015", "2016", "2017","2018"))

nmds_prokaryotes <- ggplot(NMDS_16S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Prokaryotes sequencing data") + 
  annotate("text", label = "stress = 0.18", x = 1.1, y = -1.2, size = 4, colour = "black") +
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
         plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) #center plot title and set font size

nmds_prokaryotes
ggsave("R_Code_and_Analysis/figs/NMDS_prokaryotes.png", plot = nmds_prokaryotes, width=250, height=200, units="mm",dpi=300)

### PERMANOVA 16S ###
#### LOG-transformation
log_16S <- log1p(abundances_16S_NMDS)
####  Distance matrix * this is using Bray-Curtis  ###
bray_16S <- vegdist(log_16S ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_16S <-adonis2(bray_16S ~ region + year,
                                         data=microbes_16S_ASV, permutations=999, by = "margin")

#### run PERMDISP
microbes_16S_ASV$region
attach(microbes_16S_ASV)
permdisp_16S <- betadisper(bray_16S, region, type = c("median","centroid")) 
plot(permdisp_16S)
boxplot(permdisp_16S,
        par(cex.lab=1.5))
permutest(permdisp_16S, pairwise = TRUE)

#############################################
############ 18S microeukaryotes ############
#############################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", header=T)
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
NMDS_18S$region <- factor(NMDS_18S$region, levels=c("choked", "pruth", "triquet","goose","mcmullin"))

NMDS_18S$year <- factor(NMDS_18S$year, levels=c("2015", "2016", "2017","2018"))

nmds_microeukaryotes <- ggplot(NMDS_18S, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Microeukaryotes sequencing data") + 
  annotate("text", label = "stress = 0.18", x = 1.3, y = -2.1, size = 4, colour = "black") +
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
         plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) #center plot title and set font size

nmds_microeukaryotes
ggsave("R_Code_and_Analysis/figs/NMDS_microeukaryotes.png", plot = nmds_microeukaryotes, width=250, height=200, units="mm",dpi=300)

### PERMANOVA 18S ###
#### LOG-transformation
log_18S <- log1p(abundances_18S_NMDS)
####  Distance matrix * this is using Bray-Curtis  ###
bray_18S <- vegdist(log_18S ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_18S <-adonis2(bray_18S ~ region + year,
                        data=microbes_18S_ASV, permutations=999, by = "margin")

#### run PERMDISP
microbes_18S_ASV$region
attach(microbes_18S_ASV)
permdisp_18S <- betadisper(bray_18S, region, type = c("median","centroid")) 
plot(permdisp_18S)
boxplot(permdisp_18S,
        par(cex.lab=1.5))
permutest(permdisp_18S, pairwise = TRUE)

#################################################
############ Inverts Macroeukaryotes ############
#################################################

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

# filter taxa and sites
mfilt <- m %>%
  select( year, group, region, site, sample, ID, taxon = taxon4, remove, size ) %>% 
  filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean <- m.sum %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance)) 

# make a community dataset
m.meta <- m.sum %>% 
  spread( taxon, abundance, fill=0 ) 

names(m.meta)[1:16]

m.meta <- m.meta %>% 
  ungroup( year, group, region, site, taxon )

### Creating an object to store abundances only
abundances_inverts_NMDS <- m.meta %>% 
  dplyr::select(-(1:6))

### Get MDS stats
set.seed(2)
NMDS.inverts.LOG <- metaMDS(log(abundances_inverts_NMDS+1), distance = "bray", k=2)  
NMDS.inverts.LOG 

stressplot(NMDS.inverts.LOG)
plot(NMDS.inverts.LOG)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.inverts.LOG$points[,1]
MDS2 = NMDS.inverts.LOG$points[,2]
NMDS_inverts = data.frame(MDS1 = MDS1, MDS2 = MDS2, year = m.meta$year, region = m.meta$region)
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
  annotate("text", label = "stress = 0.19", x = 1.3, y = -1.5, size = 4, colour = "black") +
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
         plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) #center plot title and set font size

nmds_macroeukaryotes
ggsave("R_Code_and_Analysis/figs/NMDS_macroeukaryotes.png", plot = nmds_macroeukaryotes, width=250, height=200, units="mm",dpi=300)

### macroeukaryotes (inverts) without 2014
NMDS_inverts_no_2014 <- NMDS_inverts %>% 
  filter(!year == 2014)

nmds_macroeukaryotes_no_2014 <- ggplot(NMDS_inverts_no_2014, aes(x=NMDS1, y=NMDS2, shape = year, colour=region)) +
  stat_ellipse(aes(colour =region, group = region), type = "t", linetype = 3, size = 1) +
  geom_point(size = 5, alpha = 0.8) +
  ggtitle("Macroeukaryotes") + 
  annotate("text", label = "stress = 0.19", x = 1.3, y = -1.5, size = 4, colour = "black") +
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
         plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) #center plot title and set font size

nmds_macroeukaryotes_no_2014
ggsave("R_Code_and_Analysis/figs/NMDS_macroeukaryotes_no_2014.png", plot = nmds_macroeukaryotes_no_2014, width=250, height=200, units="mm",dpi=300)

### PERMANOVA inverts ###
#### LOG-transformation
log_inverts <- log1p(abundances_inverts_NMDS)
####  Distance matrix * this is using Bray-Curtis  ###
bray_inverts <- vegdist(log_inverts ,m="bray")

#### run PERMANOVA with region and year - marginal
permanova_inverts <-adonis2(bray_inverts ~ region + year,
                        data=m.meta, permutations=999, by = "margin")

#### run PERMDISP
m.meta$region
m.meta$year
attach(m.meta)
permdisp_inverts <- betadisper(bray_inverts, region, type = c("median","centroid")) 
plot(permdisp_inverts)
boxplot(permdisp_inverts,
        par(cex.lab=1.5))
permutest(permdisp_inverts, pairwise = TRUE)

permdisp_inverts_year <- betadisper(bray_inverts, year, type = c("median","centroid")) 
plot(permdisp_inverts_year)
boxplot(permdisp_inverts_year,
        par(cex.lab=1.5))
permutest(permdisp_inverts_year, pairwise = TRUE)
