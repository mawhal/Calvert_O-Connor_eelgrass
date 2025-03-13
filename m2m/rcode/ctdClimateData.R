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

# ggplot(ctdSub) +
#   geom_point(aes(x = year, y = Dissolved.O2..mL.L., col = Station)) +
#   facet_wrap(vars(site))

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

### Cleaning data to maximize the data available:
#1. subset depths to be ~ 2m; all stations have min Depth..m. less than but ~ 1m and greater than 2m
#aggregate(ctdSub["Depth..m."], ctdSub[c("Station","site")], FUN = min)

ctdSub$approxDepth <- round(ctdSub$Drop.depth..m., 0)
#ctd2 <- subset(ctdSub, Drop.depth..m.> 1.5 & Drop.depth..m.< 2.5 & Depth..m.> 1.9 & Depth..m.< 2.1)
ctd2 <- subset(ctdSub, Depth..m.> 1.9 & Depth..m.< 2.1 )
ctd2$dateType <- as.Date(ctd2$date)

# Pruth:
pruth <- subset(ctd2, site == "pruth")
sort(unique(pruth$year)) #"2015" "2016" "2017" "2018"
sort(unique(pruth$month))
# KC1 sampled all year for several years---great data

aggregate(pruth["date"], pruth[c("Station","year", "month")], FUN = length)
aggregate(pruth[c("date", "Temperature..deg.C.", "Dissolved.O2..mL.L.", "Salinity..PSU.")], 
          pruth[c("Station")], FUN = length)

#   Station  date      Temperature..deg.C. Dissolved.O2..mL.L. Salinity..PSU.
# 1     KC1  711                 711                 711            711
# 2    KC13   58                  58                  58             58
# 3     KC4   47                  47                  47             47

pdf("m2m/figures/pruthClimate/pruthClim.pdf", width = 12, height = 4)
ggplot(pruth) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(Station))

ggplot(pruth) +
  geom_point(aes(x = dateType, y = Dissolved.O2..mL.L., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(Station))

ggplot(pruth) +
  geom_point(aes(x = dateType, y = Salinity..PSU., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(Station))
dev.off()

# Choked:
choked <- subset(ctd2, site == "choked")

choked <- subset(ctd2, site == "choked")
sort(unique(choked$year)) #"2015" "2016" "2017" "2018"
sort(unique(choked$month)) # "03" "04" "05" "06" "07" "08" "09" "10"
# sampling: 1-4 times a month per year from March to October---seems pretty comparable in terms of sampling frequency
# KFPS04 has the most data, but by 7 points

# Station date Temperature..deg.C. Dissolved.O2..mL.L. Salinity..PSU.
# 1  KFPS04   47                  47                  47             47
# 2  KFPS08   40                  40                  40             40

aggregate(choked["date"], choked[c("Station","year", "month")], FUN = length)
aggregate(choked[c("date", "Temperature..deg.C.", "Dissolved.O2..mL.L.", "Salinity..PSU.")], 
          choked[c("Station")], FUN = length)

pdf("m2m/figures/chokedClimate/chokedClim.pdf", width = 12, height = 4)
par(mfrow = c(1,3))
ggplot(choked) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(choked) +
  geom_point(aes(x = dateType, y = Dissolved.O2..mL.L., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(choked) +
  geom_point(aes(x = dateType, y = Salinity..PSU., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Triquet:
triquet <- subset(ctd2, site == "triquet")
sort(unique(triquet$year)) #"2015" "2016" "2017" "2018"
sort(unique(triquet$month))
# "03" "04" "05" "06" "07" "08" "09"
# sampling done 1-4 times a month
# KFPC06 has the most data 

aggregate(triquet["date"], triquet[c("Station","year", "month")], FUN = length)
aggregate(triquet[c("date", "Temperature..deg.C.", "Dissolved.O2..mL.L.", "Salinity..PSU.")], 
          triquet[c("Station")], FUN = length)

# Station date Temperature..deg.C. Dissolved.O2..mL.L. Salinity..PSU.
# 1  KFPC06   33                  33                  33             33
# 2  MACRO1   30                  30                  30             30
# 3   SEA11    2                   2                   2              2

pdf("m2m/figures/triquetClimate/triquetClim.pdf", width = 12, height = 4)
ggplot(triquet) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(triquet) +
  geom_point(aes(x = dateType, y = Dissolved.O2..mL.L., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(triquet) +
  geom_point(aes(x = dateType, y = Salinity..PSU., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()


# Goose:
rmSite <- c("pruth", "choked", "triquet", "mcMullinN")
goose <- ctd2[!ctd2$site %in% rmSite,]

sort(unique(goose$year)) #"2015" "2016" "2017" "2018"
sort(unique(goose$month))
# "03" "04" "06" "07" "08" 
# sampling done once a month
# KFPC06 has the most data 

aggregate(goose["date"], goose[c("Station","year", "month")], FUN = length)
aggregate(goose[c("date", "Temperature..deg.C.", "Dissolved.O2..mL.L.", "Salinity..PSU.")], 
          goose[c("Station", "year")], FUN = length)

# Station date Temperature..deg.C. Dissolved.O2..mL.L. Salinity..PSU.
# 1  MACRO9    6                   6                   6              6
# 2   SEA06    3                   3                   3              3
# 3   SEA07    5                   5                   5              5
# I think we would want to combine data from all of them---Macro09 = early data

pdf("m2m/figures/gooseClimate/gooseClim.pdf", width = 12, height = 4)
ggplot(goose) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(goose) +
  geom_point(aes(x = dateType, y = Dissolved.O2..mL.L., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(goose) +
  geom_point(aes(x = dateType, y = Salinity..PSU., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

# McMullin:
mcmullin <- subset(ctd2, site == "mcMullinN")
sort(unique(mcmullin$year)) #"2016" "2017"
sort(unique(mcmullin$month))
#yikes, just two data points

aggregate(mcmullin["date"], mcmullin[c("Station","year", "month")], FUN = length)
aggregate(mcmullin[c("date", "Temperature..deg.C.", "Dissolved.O2..mL.L.", "Salinity..PSU.")], 
          mcmullin[c("Station")], FUN = length)


pdf("m2m/figures/mcMullinClimate/mcmullinNClim.pdf", width = 12, height = 4)
ggplot(mcmullin) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(mcmullin) +
  geom_point(aes(x = dateType, y = Dissolved.O2..mL.L., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(mcmullin) +
  geom_point(aes(x = dateType, y = Salinity..PSU., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

#############################################################################
# 2014 ERDDAP data:

erddapFiles <- list.files(path = "m2m/input/erddap", pattern =".csv" )
erddapFiles

erddapDat <- vector()

for(i in 1: length(erddapFiles)){
  temp <- read.csv(paste("m2m/input/erddap/", erddapFiles[i], sep = ""))
  temp$site <- erddapFiles[i]
  temp$site <- gsub(".csv","", temp$site)
  
  erddapDat <- rbind(erddapDat, temp)
}

head(erddapDat)
unique(erddapDat$site)

erddapDat$time <- gsub("T","-", erddapDat$time)
temp <- str_split_fixed(erddapDat$time, "-", 4)
erddapDat$year <- as.character(temp[,1])
erddapDat$month <- as.character(temp[,2])
erddapDat$day <- as.character(temp[,3])
erddapDat$dateType <- as.Date(paste(erddapDat$year, erddapDat$month, erddapDat$day, sep = "-"))

#1. choked
choked <- subset(erddapDat, site == "choked")
sort(unique(choked$year)) #"2014" "2015" "2016" "2017" "2018" "2019"
sort(unique(choked$month))
# "03" "04" "05" "06" "07" "08" "09"
# sampling done 1-4 times a month
# KFPC06 has the most data 

aggregate(choked["dateType"], choked[c("station","year", "month")], FUN = length)
aggregate(choked[c("dateType", "temperature", "dissolved_oxygen_ml_l", "salinity")], 
          choked[c("station")], FUN = length)

# station dateType temperature dissolved_oxygen_ml_l salinity
#1    KC12       91          91                    91       91

pdf("m2m/figures/chokedClimate/chokedClimERDDAP.pdf", width = 5, height = 4)
ggplot(choked, aes(x = dateType, y = temperature, col = station)) +
  geom_point() + 
  geom_smooth(method='lm') +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(choked, aes(x = dateType, y = dissolved_oxygen_ml_l, col = station)) +
  geom_point() +   
  geom_smooth(method='lm') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(choked, aes(x = dateType, y = salinity, col = station)) +
  geom_point() +   
  geom_smooth(method='lm') +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

# looking at ERDDAP data but using different depth thresholds:

range(choked$depth)

temp <- aggregate(choked[c("depth")], 
          choked[c("dateType","station")], FUN = min); range()

# depths interested in: 2, 3,4,5,6
minDepth <- c(1.8, 2.8, 3.8, 4.8, 5.8)
maxDepth <- c(2.2, 3.2, 4.2, 5.2, 6.2)


pdf("m2m/figures/stationSelection/chokedClimate/chokedClimERDDAPDepths.pdf", width = 12, height = 10)

for(i in 1:length(minDepth)){
  temp <- subset(choked, depth> minDepth[i] & depth < maxDepth[i] )
  
  t <- ggplot(temp, aes(x = dateType, y = temperature, col = station)) +
  geom_point() + 
  geom_smooth(method='lm') +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  annotate("text", x= temp$dateType[50], y=20, label= paste(minDepth[i], maxDepth[i], nrow(subset(temp, !is.na(temperature))), sep ="_")) + 
  facet_wrap(vars(station))
  plot(t)
}
dev.off()

# 2. Pruth
pruth <- subset(erddapDat, site == "pruth")
sort(unique(pruth$year)) #"2014" "2015" "2016" "2017" "2018" "2019"
sort(unique(pruth$month))
# "03" "04" "05" "06" "07" "08" "09"
# sampling done 1-4 times a month
# KFPC06 has the most data 

aggregate(pruth["dateType"], pruth[c("station","year", "month")], FUN = length)
aggregate(pruth[c("dateType", "temperature", "dissolved_oxygen_ml_l", "salinity")], 
          pruth[c("station")], FUN = length)

# station dateType temperature dissolved_oxygen_ml_l salinity
#1    KC12       91          91                    91       91

pdf("m2m/figures/pruthClimate/pruthClimERDDAP.pdf", width = 5, height = 4)
ggplot(pruth, aes(x = dateType, y = temperature, col = station)) +
  geom_point() + 
  geom_smooth(method='lm') +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(pruth, aes(x = dateType, y = dissolved_oxygen_ml_l, col = station)) +
  geom_point() +   
  geom_smooth(method='lm') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(pruth, aes(x = dateType, y = salinity, col = station)) +
  geom_point() +   
  geom_smooth(method='lm') +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

## look at depth 
range(pruth$depth)

pdf("m2m/figures/stationSelection/pruthClimate/pruthClimERDDAPDepths.pdf", width = 12, height = 10)

for(i in 1:length(minDepth)){
  temp <- subset(pruth, depth> minDepth[i] & depth < maxDepth[i] )
  
  t <- ggplot(temp, aes(x = dateType, y = temperature, col = station)) +
    geom_point() + 
    geom_smooth(method='lm') +   
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    annotate("text", x= temp$dateType[50], y=20, label= paste(minDepth[i], maxDepth[i], nrow(subset(temp, !is.na(temperature))), sep ="_")) + 
    facet_wrap(vars(station))
  plot(t)
}
dev.off()

triquet <- subset(erddapDat, site == "triquet")
sort(unique(triquet$year)) #"2014" "2015" "2016" "2017" "2018" 
sort(unique(triquet$month))
# "03" "04" "05" "06" "07" "08" "09"
# sampling done 1-4 times a month
# KFPC06 has the most data 

aggregate(triquet["dateType"], triquet[c("station","year", "month")], FUN = length)
aggregate(triquet[c("dateType", "temperature", "dissolved_oxygen_ml_l", "salinity")], 
          triquet[c("station")], FUN = length)

#   station dateType temperature dissolved_oxygen_ml_l salinity
# 1   FZH04       12          12                    12       12
# 2  KFPC04        7           7                     7        7
# 3  KFPC06       37          37                    37       37

pdf("m2m/figures/triquetClimate/triquetClimERDDAP.pdf", width = 5, height = 4)
ggplot(triquet, aes(x = dateType, y = temperature, col = station)) +
  geom_point() + 
  geom_smooth(method='lm') +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(triquet, aes(x = dateType, y = dissolved_oxygen_ml_l, col = station)) +
  geom_point() +   
  geom_smooth(method='lm') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(triquet, aes(x = dateType, y = salinity, col = station)) +
  geom_point() +   
  geom_smooth(method='lm') +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

## Look into different depths
range(triquet$depth)

pdf("m2m/figures/stationSelection/triquetClimate/triquetClimERDDAPDepths.pdf", width = 4, height = 5)

for(i in 1:length(minDepth)){
  temp <- subset(triquet, depth> minDepth[i] & depth < maxDepth[i] )
  
  t <- ggplot(temp, aes(x = dateType, y = temperature, col = station)) +
    geom_point() + 
    geom_smooth(method='lm') +   
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    annotate("text", x= temp$dateType[5], y=18, label= paste(minDepth[i], maxDepth[i], nrow(subset(temp, !is.na(temperature))), sep ="_"))
  plot(t)
}
dev.off()
#######################################################
# Goose & mcmullins (same data)
mcmullins_goose <- subset(erddapDat, site == "goose"); unique(mcmullins_goose$station)
sort(unique(mcmullins_goose$year)) #"2014" "2016" "2017" 
sort(unique(mcmullins_goose$month))
#  "05" "06" "07" "08" 
# sampling done 1-4 times a month
# KFPC06 has the most data 

aggregate(mcmullins_goose["dateType"], mcmullins_goose[c("station","year", "month")], FUN = length)
aggregate(mcmullins_goose[c("dateType", "temperature", "dissolved_oxygen_ml_l", "salinity")], 
          mcmullins_goose[c("station")], FUN = length)

#    station dateType temperature dissolved_oxygen_ml_l salinity
# 1   FZH10        2           2                     2        2
# 2  KFPC01        7           7                     7        7

pdf("m2m/figures/gooseClimate/mcmullins_gooseClimERDDAP.pdf", width = 5, height = 4)
ggplot(mcmullins_goose, aes(x = dateType, y = temperature, col = station)) +
  geom_point() + 
  geom_smooth(method='lm') +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(mcmullins_goose, aes(x = dateType, y = dissolved_oxygen_ml_l, col = station)) +
  geom_point() +   
  #  geom_smooth(method='lm') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(mcmullins_goose, aes(x = dateType, y = salinity, col = station)) +
  geom_point() +   
  #  geom_smooth(method='lm') +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

range(goose$depth)

pdf("m2m/figures/stationSelection/gooseClimate/gooseClimERDDAPDepths.pdf", width = 4, height = 5)

for(i in 1:length(minDepth)){
  temp <- subset(mcmullins_goose, depth> minDepth[i] & depth < maxDepth[i] )
  
  t <- ggplot(temp, aes(x = dateType, y = temperature, col = station)) +
    geom_point() + 
    geom_smooth(method='lm') +   
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    annotate("text", x= temp$dateType[5], y=18, label= paste(minDepth[i], maxDepth[i], nrow(subset(temp, !is.na(temperature))), sep ="_"))
  plot(t)
}
dev.off()

####################################################################################
# But how close are the stations that actually have data?
ctdLL <- ctd2[,c("site", "Station", "Station.Longitude","Station.Latitude")]
ctdLL$data <- "ctd"
names(ctdLL) <- c("site", "station", "longitude","latitude", "data")

erddapLL <- erddapDat[,c("site", "station", "longitude","latitude")]
erddapLL$data <- "erddap"
erddapLL$site[erddapLL$site == "mcmullins_goose"] <- "gooseSW"

coord <- read.csv("Data/abiotic/hakai_abiotic_2015.csv")
coord <- unique(coord[,c("region", "latitude", "longitude")])
names(coord) <- c("site",  "latitudeSite","longitudeSite")

coord$site[coord$site == "goose"] <- "gooseSW"
coord$site[coord$site == "mcmullins"] <- "mcMullinN"

compLL <- unique(rbind(ctdLL, erddapLL))
compLL <- unique(merge(compLL, coord, by = "site"))

compLL$distLat <- abs(compLL$latitude - compLL$latitudeSite) # 0.000220 0.026303
compLL$distLong <- abs(compLL$longitude - compLL$longitudeSite) # ranges from 0.00003 0.52953

####################################################################################

# Temperature vs sites
ggplot(ctd2) +
  geom_point(aes(x = dateType, y = Dissolved.O2..mL.L., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(site))

ggplot(ctd2) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(site))

ggplot(ctd2) +
  geom_point(aes(x = dateType, y = Salinity..PSU., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(site))
  
temp <- subset(ctd2, site == "mcMullinN")

ggplot(temp) +
  geom_point(aes(x = dateType, y = Temperature..deg.C., col = Station)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
subby <- ctd2[,c("Station", "site", "year")]
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

stations <- c("KFPS04","SEA06","SEA07","SEA08", "KC1", "KC4", "KFPC06","SEA11")

ctdStn <- ctdSub[ctdSub$Station %in% stations, ]

## What does the time series data look like per station?
ctdStn$dateType <- as.Date(ctdStn$date)

for(i in 1:length(sites)){
  temp <- subset(ctdStn, Station == stations[i])
  pdf(file = paste("m2m/figures/tempLogger/", sites[i], ".pdf", sep = ""), height = 3, width = 10)
  ggplot(temp) +
    geom_point(aes(x = dateType, y = WaterTemp)) + labs( x = "Time", y = "Temperature") +
    labs (title = sites[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  dev.off()
  
}

## Comparing relationships between parameters

pdf("m2m/figures/ctdClim_Latitue.pdf", width = 12, height = 4)
par(mfrow = c(1,3))
plot(ctdStn$Dissolved.O2..mL.L ~ ctdStn$Station.Latitude, xlab = "Latitude", ylab = "Dissolved O2")
plot(ctdStn$Salinity..PSU.~ ctdStn$Station.Latitude, xlab = "Latitude", ylab = "Salinity")
plot(ctdStn$Temperature..deg.C.~ ctdStn$Station.Latitude, xlab = "Latitude", ylab = "Temperature")
dev.off()

pdf("m2m/figures/ctdClim_Correlations.pdf", width = 12, height = 4)
par(mfrow = c(1,2))
plot(ctdStn$Dissolved.O2..mL.L~ ctdStn$Salinity..PSU., xlab = "Salinity", ylab = "Dissolved O2"); abline(a= 0, b =1)
plot(ctdStn$Salinity..PSU., ctdStn$Temperature..deg.C., xlab = "Salinity", ylab = "Temperature"); abline(a= 0, b =1)
dev.off()

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


## correlations between parameters:
# relationship between LAI and biomass
plot(seagrass$quadrat_biomass_g ~ seagrass$quadrat_lai, col = as.factor(seagrass$site), pch = 19); abline(0,1)
plot(seagrass$quadrat_biomass_g ~ seagrass$quadrat_shoot_density, col = as.factor(seagrass$site), pch = 19); abline(0,1)
plot(seagrass$quadrat_biomass_g ~ seagrass$quadrat_macroalgae_g, col = as.factor(seagrass$site), pch = 19); abline(0,1)
plot(seagrass$quadrat_biomass_g ~ seagrass$quadrat_microepiphyte_mg, col = as.factor(seagrass$site), pch = 19); abline(0,1)

plot(seagrass$quadrat_microepiphyte_mg ~ seagrass$quadrat_macroalgae_g, col = as.factor(seagrass$site), pch = 19); abline(0,1)
plot(seagrass$quadrat_microepiphyte_mg ~ seagrass$quadrat_lai, col = as.factor(seagrass$site), pch = 19); abline(0,1)
plot(seagrass$quadrat_macroalgae_g ~ seagrass$quadrat_lai, col = as.factor(seagrass$site), pch = 19); abline(0,1)


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