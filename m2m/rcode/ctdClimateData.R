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

chokedDrops <- read.csv("m2m/input/chokedDrops.csv")
chokedData <- read.csv("m2m/input/chokedData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
choked <- merge(chokedDrops, chokedData, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))
choked$site <- "choked"

gooseSEDrops <- read.csv("m2m/input/gooseSEDrops.csv")
gooseSEData <- read.csv("m2m/input/gooseSEData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
gooseSE <- merge(gooseSEDrops, gooseSEData, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))
gooseSE$site <- "gooseSE"

gooseSWDrops <- read.csv("m2m/input/gooseSWDrops.csv")
gooseSWData <- read.csv("m2m/input/gooseSWData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
gooseSW <- merge(gooseSWDrops, gooseSWData, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))
gooseSW$site <- "gooseSW"

mcMullinNDrops <- read.csv("m2m/input/mcMullinNDrops.csv")
mcMullinNData <- read.csv("m2m/input/mcMullinNData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
mcMullinN <- merge(mcMullinNDrops, mcMullinNData, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))
mcMullinN$site <- "McMullinN"

pruthDrops <- read.csv("m2m/input/pruthDrops.csv")
pruthData <- read.csv("m2m/input/pruthData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
pruth <- merge(pruthDrops, pruthData, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))
pruth$site <- "pruth"

triquetDrops <- read.csv("m2m/input/triquetDrops.csv")
triquetData <- read.csv("m2m/input/triquetData.csv")

# shared columns: "Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude" 
triquet <- merge(triquetDrops, triquetData, by = c("Cast.PK", "Cruise", "Station", "Station.Longitude", "Station.Latitude"))
triquet$site <- "triquet"

ctd <- rbind(choked, gooseSE, gooseSW, mcMullinN, pruth, triquet)

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

# zoomed in maps:
#pdf("m2m/figures/siteMap2.pdf", height = 6, width = 6)
#choked
ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = ctdCoord, aes(x = Station.Longitude, y = Station.Latitude), size = 4, 
             shape = 23, fill = "purple4") +
   geom_label_repel(data = siteCoord,
                    aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  geom_label_repel(data = ctdCoord,
                    aes(x = Station.Longitude, y = Station.Latitude, label = Station),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-128.11, -128.13), ylim = c(51.67,51.69), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#Goose SE
ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = ctdCoord, aes(x = Station.Longitude, y = Station.Latitude), size = 4, 
             shape = 23, fill = "purple4") +
  # geom_label_repel(data = siteCoord,
  #                  aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  # geom_label_repel(data = ctdCoord,
  #                  aes(x = Station.Longitude, y = Station.Latitude, label = Station),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-128.44, -128.48), ylim = c(51.91,51.93), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#McMullins
ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = ctdCoord, aes(x = Station.Longitude, y = Station.Latitude), size = 4, 
             shape = 23, fill = "purple4") +
  geom_label_repel(data = siteCoord,
                   aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  geom_label_repel(data = ctdCoord,
                   aes(x = Station.Longitude, y = Station.Latitude, label = Station),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-128.11, -128.13), ylim = c(51.67,51.69), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#Triquet
ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = ctdCoord, aes(x = Station.Longitude, y = Station.Latitude), size = 4, 
             shape = 23, fill = "purple4") +
  geom_label_repel(data = siteCoord,
                   aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  geom_label_repel(data = ctdCoord,
                   aes(x = Station.Longitude, y = Station.Latitude, label = Station),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-128.11, -128.13), ylim = c(51.67,51.69), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#Pruth
ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = ctdCoord, aes(x = Station.Longitude, y = Station.Latitude), size = 4, 
             shape = 23, fill = "purple4") +
  geom_label_repel(data = siteCoord,
                   aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  geom_label_repel(data = ctdCoord,
                   aes(x = Station.Longitude, y = Station.Latitude, label = Station),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-128.11, -128.13), ylim = c(51.67,51.69), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = ctdCoord, aes(x = Station.Longitude, y = Station.Latitude), size = 4, 
             shape = 23, fill = "purple4") +
  geom_label_repel(data = siteCoord,
                   aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  geom_label_repel(data = ctdCoord,
                   aes(x = Station.Longitude, y = Station.Latitude, label = Station),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-128.11, -128.13), ylim = c(51.67,51.69), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
