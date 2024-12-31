# started Nov 12, 2024 by DL
# Aim: to map the sites at which the microbial and macro community are observed as part of the M2M project

# sites based on summary from communityDataExplore.R
# but where can I find their lat longs?

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) 
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)

#setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")
if(length(grep("deirdre", getwd())>0)) {setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {setwd("~/Documents/git/temp")
}

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

#abiotic <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_abiotic_20200214.csv")
#abiotic <- read.csv("metadata/coordinates.csv")
abiotic <- read.csv("Data/R_Code_for_Data_Prep/master_data/(MASTER_merged_explanatory_20200)(214.csv")
siteCoord <- unique(abiotic[, c("site","lat", "long")])
siteCoord <- subset(siteCoord, !is.na(lat))
siteCoord <- siteCoord[siteCoord$site %in% datSites, ]

state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"))

pdf("m2m/figures/siteMap2.pdf", height = 6, width = 6)
ggplot(data = state_prov) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = siteCoord, aes(x = long, y = lat), size = 4, 
             shape = 23, fill = "chartreuse4") +
  # geom_label_repel(data = siteCoord,
  #                  aes(x = long, y = lat, label = site),hjust=1,vjust=0) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim = c(-129, -127), ylim = c(51.4,52.5), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()
