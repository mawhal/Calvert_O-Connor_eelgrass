# started Nov 12, 2024 by DL
# Aim: to map the sites at which the microbial and macro community are observed as part of the M2M project

# sites based on summary from communityDataExplore.R
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

# but where can I find their lat longs?

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) 

#setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")
if(length(grep("deirdre", getwd())>0)) {setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {setwd("~/Documents/git/temp")
}

abiotic <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_abiotic_20200214.csv")
siteCoord <- unique(abiotic[, c("site","lat", "long")])
siteCoord <- subset(siteCoord, !is.na(lat))
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

(west <- data.frame(Longitude = c(-127.1686, -120.7816,-122.1417, - 119.8891), Latitude = c(54.7824, 49.0646,52.1417,50.8837)))
(east <- data.frame(Longitude = c(-72.1900,-74.0248,-71.09543597, -71.37388296), Latitude = c(42.5315, 45.9310,44.92466697,43.99837498)))

# (72.1900,74.0248) (42.5315, 45.9310)
pdf("figures/phenoTraitMap.pdf", height = 6, width = 6)
ggplot(data = world) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = west, aes(x = Longitude, y = Latitude), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = east, aes(x = Longitude, y = Latitude), size = 4, 
             shape = 23, fill = "chartreuse4") +
  coord_sf(xlim = c(-140.15, -54.12), ylim = c(38, 71), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()
