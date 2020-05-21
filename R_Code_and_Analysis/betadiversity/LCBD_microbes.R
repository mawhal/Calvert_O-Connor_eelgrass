####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will investigate beta diversity of grazer (and microbe)
###        community data
###
### started 3 March 2020 by Matt WHalen
### modified for microbes 4 March 2020 by Bia
###  
####################################################################


library(tidyverse)
library(adespatial)


## Goal calculate local contributions to beta diversity in different years (do we get the same answer?)

# steps - split data into years, make a list of community matrices, perform calculation

################################################
################# 16S dataset ##################
################################################

### read data
# read data
zostera_16S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", header=T)
 
zostera_16S <- zostera_16S[order(zostera_16S$year, zostera_16S$site),]
 names(zostera_16S)
 # split by year
comm.years <- split( zostera_16S, zostera_16S$year )
#lapply(comm.zero, dim)

# remove taxa that did not appear in a given year
# lapply( comm.years, function(z) = select this and run first
comm.zero <- lapply(comm.years, function(z) {
  tmp <- z[,-c(1:16)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:16)],df))
})

meta <- lapply(comm.zero, function(z) z[,7])
d <- lapply(comm.zero, function(z) z[,-c(1:16)])

# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBD <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBD <- data.frame( do.call( rbind, LCBD ) )
names(LCBD) <- c("LCBD","p.adj")
LCBD <- data.frame( zostera_16S[,1:13], LCBD )

#windows(3.5,4)
LCBD_16S <- ggplot( data=LCBD, aes( x=factor(year), y=factor(site), size=LCBD, fill=p.adj )) +
  geom_point(shape=21,alpha=0.5,position=position_jitter(width=0.1,height = 0.1)) +
  ylab( "Site" ) + xlab( "Year" ) + ggtitle("Prokaryotes") + theme_bw( )

ggsave("R_Code_and_Analysis/betadiversity/LCBD_16S.png", plot = LCBD_16S, width=250, height=200, units="mm",dpi=300)

################################################
################# 18S dataset ##################
################################################

### read data
# read data
zostera_18S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", header=T)

zostera_18S <- zostera_18S[order(zostera_18S$year, zostera_18S$site),]
names(zostera_18S)
# split by year
comm.years <- split( zostera_18S, zostera_18S$year )
#lapply(comm.zero, dim)

# remove taxa that did not appear in a given year
# lapply( comm.years, function(z) = select this and run first
comm.zero <- lapply(comm.years, function(z) {
  tmp <- z[,-c(1:13)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:13)],df))
})

meta <- lapply(comm.zero, function(z) z[,7])
d <- lapply(comm.zero, function(z) z[,-c(1:13)])

# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBD <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBD <- data.frame( do.call( rbind, LCBD ) )
names(LCBD) <- c("LCBD","p.adj")
LCBD <- data.frame( zostera_18S[,1:13], LCBD )

#windows(3.5,4)
LCBD_18S <- ggplot( data=LCBD, aes( x=factor(year), y=factor(site), size=LCBD, fill=p.adj )) +
  geom_point(shape=21,alpha=0.5,position=position_jitter(width=0.1,height = 0.1)) +
  ylab( "Site" ) + xlab( "Year" ) + ggtitle("Microeukaryotes") + theme_bw( )

ggsave("R_Code_and_Analysis/betadiversity/LCBD_18S.png", plot = LCBD_18S, width=250, height=200, units="mm",dpi=300)
