# started Nov 27 2024
# aim of this code is to compile as assess the CTD data from the Hakai data portal
# Zach M. very kindly provided a table of potentially relevant sites near where the seagrass was collected
# BUT they differ in the timing and frequency of sampling---let some may be more useful than others

# 1. read in the data and compile for each site:

# 2. Make a map---what sites are closes to our sampling sites?

# 3. Which sites have the best (most frequent) sampling frequency for the range of dates we are interested in?

rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) 
library(ggrepel)
require(stringr)

if(length(grep("deirdre", getwd())>0)) {  setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {   setwd("~/Documents/git/temp")
}

# read in an merge the different site data files
# 1. Data
filesData <- list.files(path = "m2m/input/ctdData", pattern =".csv" )
filesData

ctdDat <- vector()

for(i in 1: length(filesData)){
temp <- read.csv(paste("m2m/input/ctdData/", filesData[i], sep = ""))
temp$site <- filesData[i]
temp$site <- gsub("Data.csv","", temp$site)

ctdDat <- rbind(ctdDat, temp)
}

# 2. Drops
filesDrops <- list.files(path = "m2m/input/ctdDrops", pattern =".csv" )
filesDrops

ctdDrop <- vector()

for(i in 1: length(filesDrops)){
  temp <- read.csv(paste("m2m/input/ctdDrops/", filesDrops[i], sep = ""))
  temp$site <- filesDrops[i]
  temp$site <- gsub("Drops.csv","", temp$site)
  
  ctdDrop <- rbind(ctdDrop, temp)
}


# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
ctd<- merge(ctdDrop, ctdDat, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude", "site"))

temp <- str_split_fixed(ctd$Measurement.time, " ", 2)
ctd$date <- as.character(temp[,1])

temp2 <- data.frame(str_split_fixed(ctd$date, "-", 3))
names(temp2) <- c("year", "month", "day")

# temp2 <- data.frame(t(sapply(ctd$date, function(x) substring(x, first=c(1,5,7), last=c(4,6,8))))); head(temp2)
# names(temp2) <- c("year", "month", "day")
 ctd <- cbind(ctd, temp2)

ctd$Station[which(ctd$Station == "SEA6" )] <- "SEA06"
ctd$Station[which(ctd$Station == "SEA7" )] <- "SEA07"
ctd$Station[which(ctd$Station == "SEA8" )] <- "SEA08"

ctd$yrMonth <- paste(ctd$year, ctd$month, sep = "_")

names(ctd)
ctdSub <- ctd[,c("site","Station", "year", "month", "day", "yrMonth", "date","Station.Longitude", "Station.Latitude", "Drop.depth..m.", "Water.Depth..m.","Depth..m.",  
                   "Dissolved.O2..mL.L.", "Salinity..PSU.", "Temperature..deg.C.")]


ggplot(ctdSub) +
  geom_point(aes(x = year, y = Dissolved.O2..mL.L., col = Station)) +
  facet_wrap(vars(site))

sites <- unique(ctdSub$site)

for(i in 1:length(sites)){
  temp <- subset(ctdSub, site == sites[i])
  climPlot <- ggplot(temp) +
    geom_point(aes(x = year, y = Temperature..deg.C., col = Station)) +
    xlab("Year") + labs(title = sites[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  pdf(file = paste("m2m/figures/climVarTemp/", sites[i], ".pdf", sep = ""), height = 3, width = 3)
  climPlot
  dev.off()
  
}

subby <- unique(ctdSub[,c("Station", "site", "year")])
aggregate(subby["year"], subby[c("Station","site")], FUN = length)

#   Station      site year
# 1   KFPS04    choked    4 2015-2018 <<- 
# 2   KFPS08    choked    4 2015-2018 # less data
# 3    SEA06   gooseSE    3 2015, 2016, 2018 <<-
# 4   MACRO9   gooseSW    2 
# 5    SEA07   gooseSW    4 2015-2018 <<-
# 6    SEA08 mcMullinN    3 2015-2017 <<-
# 7      KC1     pruth    4 2015-2018 # most data <<-
# 8     KC13     pruth    4 2015-2018 # least data
# 9      KC4     pruth    4 2015-2018 : 
# 10  KFPC06   triquet    3
# 11  MACRO1   triquet    4 2015-2018 <<-
# 12   SEA11   triquet    3 2015-2017

datSites <- c("choked_sandspit", 
              "choked_inner", 
              "goose_south_west", 
              "goose_south_east", 
              "goose_north", # lat long needed
              "mcmullins_south", 
              "mcmullins_north", 
              "triquet_north", 
              "triquet_south", 
              "pruth_pocket", 
              "pruth_bay")

# seagrass <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_merged_explanatory_20200214.csv")
seagrass <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_seagrass_metrics_20200214.csv")
#mcMullins---2015 only, goose 2015 and 2016

ggplot(seagrass) +
  geom_point(aes(x = year, y = quadrat_shoot_density, col = site)) +
  labs(y = "Shoot density per quadrat") +
  geom_smooth(method = lm, aes(x = year, y = quadrat_shoot_density, col = site)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(seagrass) +
  geom_point(aes(x = year, y = quadrat_macroalgae_g, col = site)) +
  labs(y = "Macroalgae per quadrat") +
  geom_smooth(method = lm, aes(x = year, y = quadrat_macroalgae_g, col = site)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(seagrass) +
  geom_point(aes(x = year, y = quadrat_lai, col = site)) +
  labs(y = "LAI per quadrat") +
  geom_smooth(method = lm, aes(x = year, y = quadrat_lai, col = site)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(seagrass) +
  geom_point(aes(x = year, y = quadrat_biomass_g, col = site)) +
  labs(y = "Biomass per quadrat") +
  geom_smooth(method = lm, aes(x = year, y = quadrat_biomass_g, col = site)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(seagrass) +
  geom_point(aes(x = year, y = quadrat_microepiphyte_mg, col = site)) +
  labs(y = "Microphiphyte per quadrat") +
 # geom_smooth(method = lm, aes(x = year, y = quadrat_microepiphyte_mg, col = site)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

abiotic <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_abiotic_20200214.csv")
sitesRM <- c("goose kelp", "maye kelp", "triquet kelp", "west beach kelp", "wolf")

abiotic <- abiotic[!abiotic$site %in% sitesRM, ]

ggplot(abiotic) +
  geom_point(aes(x = lat, y = salinity, col = site)) +
  labs(y = "Salinity") +
   geom_smooth(method = lm, aes(x = lat, y = salinity)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(abiotic) +
  geom_point(aes(x = site, y = depth, col = site)) +
  labs(y = "Depth") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("m2m/figures/siteDepth.pdf", width =5, height = 4)
ggplot(abiotic) +
  geom_point(aes(x = site, y = site.depth, col = site)) +
  labs(y = "Site Depth") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# site depth ranges from 2.2-9.5 depth from 0 to 7

# What does the data look like? 
# Are any of the environmental factors highly correlated?
# Run a PCA?