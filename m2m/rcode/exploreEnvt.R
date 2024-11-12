# Started Oct 9 2024 by D. Loughnan

# aim of this code is to explore the variation in the environmental variables across sites
# answer the question of whether there is enough to make it worth including in a model

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) 
library(cowplot)

#setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")
if(length(grep("deirdre", getwd())>0)) {  setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {   setwd("~/Documents/git/temp")
}

# climate data from hakai
c2015 <- read.csv("Data/abiotic/hakai_abiotic_2015.csv")
c2016 <- read.csv("Data/abiotic/hakai_abiotic_2016.csv")
c2017 <- read.csv("Data/abiotic/corrected_hakai_abiotic_2017.csv")
c2018 <- read.csv("Data/abiotic/corrected_hakai_abiotic_2018.csv")

clim2015 <- unique(c2015)
clim2016 <- unique(c2016)
clim2017 <- unique(c2017)
clim2018 <- unique(c2018)

clim2015$year <- "2015"
clim2016$year <- "2016"
clim2017$year <- "2017"
clim2018$year <- "2018"

colnames(clim2016)[colnames(clim2016) == "conductivity..us.cm."] <- "conductivity"
colnames(clim2016)[colnames(clim2016) == "depth..m."] <- "mean_depth"
colnames(clim2016)[colnames(clim2016) == "do..mg.l."] <- "dissolved_oxygen_concentration"
colnames(clim2016)[colnames(clim2016) == "do...."] <- "dissolved_oxygen"
colnames(clim2016)[colnames(clim2016) == "salinity..ppt."] <- "salinity"
colnames(clim2016)[colnames(clim2016) == "temperature..c."] <- "temperature"


colnames(clim2017)[colnames(clim2017) == "cond.uscm."] <- "conductivity"
colnames(clim2017)[colnames(clim2017) == "depth.m."] <- "mean_depth"
colnames(clim2017)[colnames(clim2017) == "do..mg.l."] <- "dissolved_oxygen_concentration"
colnames(clim2017)[colnames(clim2017) == "do...."] <- "dissolved_oxygen"
colnames(clim2017)[colnames(clim2017) == "salinity.ppt."] <- "salinity"
colnames(clim2017)[colnames(clim2017) == "temp.c."] <- "temperature"

colnames(clim2018)[colnames(clim2018) == "conductivity..spc."] <- "conductivity"
colnames(clim2018)[colnames(clim2018) == "sampling_depth..m."] <- "mean_depth"
colnames(clim2018)[colnames(clim2018) == "salinity..ppt."] <- "salinity"
colnames(clim2018)[colnames(clim2018) == "temp..c."] <- "temperature"

keep <- c("site", "date", "mean_depth", "depth.m", 
          "dissolved_oxygen_concentration", "dissolved_oxygen",  
          "salinity", "temperature", "year")
clim2015 <- clim2015[, c("site", "date", "mean_depth", 
                         #"dissolved_oxygen_concentration", "dissolved_oxygen",  
                         "salinity", "temperature", "year")]
clim2016 <- clim2016[, c("site", "date", "mean_depth", 
                         #"dissolved_oxygen_concentration", "dissolved_oxygen",  
                         "salinity", "temperature", "year")]
clim2017 <- clim2017[, c("site", "date", "mean_depth",  "salinity", "temperature", "year")]
clim2018 <- clim2018[, c("site", "date", "mean_depth",  "salinity", "temperature", "year")]

climVar <- rbind(clim2015, clim2016, clim2017, clim2018)

plot(climVar$mean_depth ~ as.factor((climVar$year)))
plot(climVar$salinity ~ as.factor((climVar$year)))
plot(climVar$temperature ~ as.factor((climVar$year)))

str(climVar)

annum <- sort(unique(climVar$year))
color <- c(rgb(72 / 255, 38 / 255, 119 / 255, alpha = 0.14), 
           rgb(149 / 255, 216 / 255, 64 / 255, alpha = 0.2),
           rgb(80 / 255, 38 / 255, 19 / 255, alpha = 0.14), 
           rgb(200 / 255, 16 / 255, 64 / 255, alpha = 0.2))

plot(1, type = 'n', xlim = c(5,20), ylim = c(0, 30))
for(i in 1:length(unique(climVar$year))){
  temp <- subset(climVar, year == annum[i] )
  hist(temp$temperature, add = T, col = color[i])
}

climTemp <- unique(climVar[, c(1,5,6)])

t2015 <- ggplot(subset(climVar, year == "2015"), aes(x = as.factor(site), y = temperature, col = "maroon")) +
  geom_point()  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
  #facet_wrap(vars(year))

t2016 <- ggplot(subset(climVar, year == "2016"), aes(x = as.factor(site), y = temperature)) +
  geom_point(aes(col = "cyan4"))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
#facet_wrap(vars(year))

t2017 <- ggplot(subset(climVar, year == "2017"), aes(x = as.factor(site), y = temperature)) +
  geom_point(aes(col = "goldenrod"))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
#facet_wrap(vars(year))

t2018 <- ggplot(subset(climVar, year == "2018"), aes(x = as.factor(site), y = temperature)) +
  geom_point(aes(col = "purple4"))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
#facet_wrap(vars(year))
pdf("m2m/figures/histTempYear2.pdf", width = 8, height = 4)
plot_grid(t2015, t2016, t2017, t2018, ncol = 4, nrow =1) 
dev.off()

ggplot(climVar, aes(x = temperature)) +
  geom_histogram(color="black", fill="white")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  facet_wrap(vars(year))

pdf("m2m/figures/histSalinityYear.pdf")
ggplot(climVar, aes(x = salinity)) +
  geom_histogram(color="black", fill="white")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  facet_wrap(vars(year))
dev.off()

ggplot(climVar, aes(x = as.factor(site), y = temperature)) +
  geom_point(aes(color= site))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  facet_wrap(vars(year))

c2015 <- unique(c2015)
c2015 <- c2015[order(c2015$latitude),]
siteOrder <- c2015$site

temp <- ggplot(c2015, aes(x = factor(site, level = siteOrder), y = temperature)) +
  geom_point(aes(color= site))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 

salinity <- ggplot(c2015, aes(factor(site, level = siteOrder), y = salinity)) +
  geom_point(aes(color= site))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 

dissO <- ggplot(c2015, aes(x = factor(site, level = siteOrder), y = dissolved_oxygen_concentration)) +
  geom_point(aes(color= site))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
range(c2015$dissolved_oxygen_concentration, na.rm = T)

bedArea <- ggplot(c2015, aes(x = factor(site, level = siteOrder), y = bed_area_m2)) +
  geom_point(aes(color= site))  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") 
range(c2015$bed_area_m2, na.rm = T)

pdf("m2m/figures/climVar2015.pdf", width = 12, height = 5)
plot_grid(temp, salinity, dissO, bedArea, ncol = 4, nrow =1) 
dev.off()
