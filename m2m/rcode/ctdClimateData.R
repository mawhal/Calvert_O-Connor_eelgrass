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
temp$site <- gsub(".csv","", temp$site)

ctdDat <- rbind(ctdDat, temp)
}

# 2. Drops
filesDrops <- list.files(path = "m2m/input/ctdDrops", pattern =".csv" )
filesDrops

ctdDrop <- vector()

for(i in 1: length(filesDrops)){
  temp <- read.csv(paste("m2m/input/ctdDrops/", filesDrops[i], sep = ""))
  temp$site <- filesDrops[i]
  temp$site <- gsub(".csv","", temp$site)
  
  ctdDrop <- rbind(ctdDrop, temp)
}

chokedDrops <- read.csv("m2m/input/chokedDrops.csv")
chokedData <- read.csv("m2m/input/chokedData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
ctd<- merge(ctdDrop, ctdDat, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))

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

# ctdChoked <- subset(ctd, site == "choked")
# 
# ggplot(ctdChoked, aes(x = yrMonth, y = Salinity..PSU.)) +
#   geom_point(aes(col = Station)) 
# 
# 
# ctdCoord <- unique(ctd[, c("Station", "Station.Longitude", "Station.Latitude")])

ctdSub <- ctd[,c("site","Station", "year", "month", "day", "yrMonth", "date","Station.Longitude", "Station.Latitude", 
                   "Dissolved.O2..mL.L.", "Salinity..PSU.", "Temperature..deg.C.")]

ctdCoord <- unique(ctdSub[, c("Station", "Station.Longitude", "Station.Latitude")])


trt.dataset <- ctdSub%>%
  group_by(site, Station, year) %>%
  summarise(no_rows = length(unique(date)), .groups = 'drop')

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
seagrass <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_abiotic_20200214.csv")
siteCoord <- unique(seagrass[, c("site","lat", "long")])
siteCoord <- subset(siteCoord, !is.na(lat))
siteCoord <- siteCoord[siteCoord$site %in% datSites, ]

state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"))

# What does the data look like? 
# Are any of the environmental factors higly correlated?
# Run a PCA?