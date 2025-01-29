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
filesData <- list.files(path = "m2m/input/tempTimeSeries", pattern =".csv" )
filesData

tempData <- vector()

for(i in 1: length(filesData)){
temp <- read.csv(paste("m2m/input/tempTimeSeries/", filesData[i], sep = ""))
temp$site <- filesData[i]
temp$site <- gsub("Data.csv","", temp$site)

tempData <- rbind(tempData, temp)
}

tempData$site <- gsub("choked_inner","chokedinner", tempData$site)
tempData$site <- gsub("choked_sandspit","chokedsandspit", tempData$site)

splitName <- str_split_fixed(tempData$site, "_", 4)
names(splitName) <- c("project", "region", "subRegion", "aux")
tempData$site <- splitName[,2]

splitDate <- str_split_fixed(tempData$date, "-", 3)
names(splitDate) <- c("year", "month","day" )

seaData <- subset(tempData, year < 2019)

# Triquet Bay 2 has more data than 1; closest o Triquest S in the m2m data
# Triquet Seagrass = close to Triquest N

rmSite <- c("PruthBay01", "TriquetBay01", "TriquetKelp")
seaData <- seaData[!seaData$site %in% rmSite,]

#temp2 <- data.frame(t(sapply(ctd$date, function(x) substring(x, first=c(1,5,7), last=c(4,6,8))))); head(temp2)
# names(temp2) <- c("year", "month", "day")


# Calculate MAP and Season temps:
mat <- aggregate(seaData["WaterTemp"], seaData[c("site", "year")], FUN = mean, na.rm = T)
colnames(mat)[colnames(mat) == "WaterTemp"] <- "MAT"

# spring: March-May
spring <- c("Mar", "Apr", "May")
springT <- seaData[seaData$month %in% spring, ]
springTemp <- aggregate(springT["WaterTemp"], springT[c("site", "year")], FUN = mean, na.rm = T)
colnames(springTemp)[colnames(springTemp) == "WaterTemp"] <- "springTemp"

# June-Aug
summer <- c("Jun", "Jul", "Aug")
summerT <- seaData[seaData$month %in% summer, ]
summerTemp <- aggregate(summerT["WaterTemp"], summerT[c("site", "year")], FUN = mean, na.rm = T)
colnames(summerTemp)[colnames(summerTemp) == "WaterTemp"] <- "summerTemp"

# growing season: april-oct
growSea <- c( "Apr", "May", "Jun", "Jul", "Aug", "Sept")
aprilSept <- seaData[seaData$month %in% growSea, ]
growSeason <- aggregate(aprilSept["WaterTemp"], aprilSept[c("site", "year")], FUN = mean, na.rm = T)
colnames(growSeason)[colnames(growSeason) == "WaterTemp"] <- "growSeaTemp"

fin <- merge(mat, springTemp, all.x = T, by = c("site", "year"))
fin <- merge(fin, summerTemp, all.x = T, by = c("site", "year"))
fin <- merge(fin, growSeason, all.x = T, by = c("site", "year"))

write.csv(fin,"m2m/output/seasonalMatTemp.csv", row.names = F)
sites <- unique(seaData$site)

sumTemp <- unique(seaData[, c("year", "site")])
sumTemp$iter <- rep(1:nrow(subby))

for (s in 1:length(sites)){
  s <- 1
  temp <- subset(seaData, site == sites[s])
  years <- unique(temp$year)
  
  for(y in 1:length(years)){
    y <- 1
    tempY <- subset(temp, year == years[y])
    mat <- mean(tempY$WaterTemp)
  }
}

pdf("m2m/figures/tempLogger.pdf", width =15, height =10)
ggplot(seaData) +
  geom_point(aes(x = date, y = WaterTemp)) +   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_x_discrete(labels = everysecond(seaData$date))+
  facet_wrap(vars(site))
dev.off()

for(i in 1:length(sites)){
  temp <- subset(seaData, site == sites[i])
  climPlot <- ggplot(temp) +
    geom_point(aes(x = date, y = WaterTemp)) +
    xlab("Year") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  pdf(file = paste("m2m/figures/tempLogger/", sites[i], ".pdf", sep = ""), height = 3, width = 15)
  climPlot
  dev.off()
  
}


## if want to later---calculating julian date
# require(lubridate)
# x = as.Date('2010-06-10')
# yday(x)
