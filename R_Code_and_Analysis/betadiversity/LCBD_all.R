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
zostera_16S <- read_csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv")
zostera_16S[1:5,1:20]

# grab relevant columns, summarize means across samples within sites
# zostera_16S <- 
wide_16S <- zostera_16S %>% 
  select( -SampleID, -swab_id, -barcode_plate, -barcode_well, -region_year,
          -host_species, -host_type, -sample_type , -survey_type, -quadrat_id,
          -meso_shoot_id ) %>% 
  gather( "ASV", "read.count", -year, -region, -site_quadrat_id, -site ) %>% 
  group_by( year, region, site, ASV ) %>% 
  summarize( read.count=mean(read.count) ) %>% 
  spread( ASV, read.count )

# split by year
comm.years <- split( wide_16S, wide_16S$year )
#lapply(comm.years, dim)

# remove taxa that did not appear in a given year
# lapply( comm.years, function(z) = select this and run first
comm.zero <- lapply(comm.years, function(z) {
  tmp <- z[,-c(1:3)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:3)],df))
})
lapply(comm.zero, dim)

meta <- lapply(comm.zero, function(z) z[,1:3])
d <- lapply(comm.zero, function(z) z[,-c(1:3)])

# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBD16 <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBD16 <- data.frame( do.call( rbind, LCBD16 ), group="prokaryote" )
names(LCBD16) <- c("LCBD","p.adj","group")
LCBD16 <- data.frame( wide_16S[,1:3], LCBD16 )

#windows(3.5,4)
LCBD_16S <- ggplot( data=LCBD16, aes( x=factor(year), y=factor(site), size=LCBD, fill=p.adj )) +
  geom_point(shape=21,alpha=1 ) +
  ylab( "Site" ) + xlab( "Year" ) + ggtitle("Prokaryotes") + theme_bw( )
LCBD_16S
# ggsave("R_Code_and_Analysis/betadiversity/LCBD_16S.png", plot = LCBD_16S, width=250, height=200, units="mm",dpi=300)



################################################
################# 18S dataset ##################
################################################

### read data
# read data
zostera_18S <- read_csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv")
zostera_18S[1:5,1:20]

# grab relevant columns, summarize means across samples within sites
# zostera_18S <- 
wide_18S <- zostera_18S %>% 
  select( -SampleID, -host_type, -sample_type , -survey_type, 
          -meso_quadrat_id ) %>% 
  gather( "ASV", "read.count", -year, -region, -site_quadrat_id, -site ) %>% 
  group_by( year, region, site, ASV ) %>% 
  summarize( read.count=mean(read.count) ) %>% 
  spread( ASV, read.count )

# split by year
comm.years <- split( wide_18S, wide_18S$year )
#lapply(comm.years, dim)

# remove taxa that did not appear in a given year
# lapply( comm.years, function(z) = select this and run first
comm.zero <- lapply(comm.years, function(z) {
  tmp <- z[,-c(1:3)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:3)],df))
})
lapply(comm.zero, dim)

meta <- lapply(comm.zero, function(z) z[,1:3])
d <- lapply(comm.zero, function(z) z[,-c(1:3)])

# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBD18 <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBD18 <- data.frame( do.call( rbind, LCBD18 ), group="microeukaryote" )
names(LCBD18) <- c("LCBD","p.adj","group")
LCBD18 <- data.frame( wide_18S[,1:3], LCBD18 )

#windows(3.5,4)
LCBD_18S <- ggplot( data=LCBD18, aes( x=factor(year), y=factor(site), size=LCBD, fill=p.adj )) +
  geom_point(shape=21,alpha=1 ) +
  ylab( "Site" ) + xlab( "Year" ) + ggtitle("Prokaryotes") + theme_bw( )
LCBD_18S
# ggsave("R_Code_and_Analysis/betadiversity/LCBD_18S.png", plot = LCBD_18S, width=250, height=200, units="mm",dpi=300)






################################################
####### Invertebrates(macroeukaryotes) #########
################################################

# steps - split data into years, make a list of community matrices, perform calculation

### read data
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
  filter( is.na(remove), !is.na(taxon), site != "goose_north" )

# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean <- m.sum %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance))


# make a community dataset
m.meta <- m.mean %>% 
  spread( taxon, abundance, fill=0 ) 

# split by year
comm.years <- split( m.meta, m.meta$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:3)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:3)],df))
})
lapply(comm.zero, dim)

meta <- lapply(comm.zero, function(z) z[,1:3])
d <- lapply(comm.zero, function(z) z[,-c(1:3)])



# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBDin <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBDin <- data.frame( do.call( rbind, LCBDin ), group="macroeukaryote" )
names(LCBDin) <- c("LCBD","p.adj","group")
LCBDin <- data.frame( m.meta[,1:3], LCBDin )

# windows(3.5,4)
LCBD_inverts <- ggplot( data=LCBDin, aes( x=factor(year), y=factor(site), size=LCBD, fill=p.adj )) +
  geom_point(shape=21,alpha=1) + #,position=position_jitter(width=0.1,height = 0.1)) +
  ylab( "Site" ) + xlab( "Year" ) + ggtitle("Macroeukaryotes") + theme_bw( )
LCBD_inverts
# ggsave("R_Code_and_Analysis/betadiversity/LCBD_inverts.png", plot = LCBD_inverts, width=250, height=200, units="mm",dpi=300)



# combine into one figure
LCBDall <- bind_rows( LCBD16, LCBD18, LCBDin )
LCBDall$site[ LCBDall$site == "mcmullin_north" ] <- "mcmullins_north"
LCBDall$site[ LCBDall$site == "mcmullin_south" ] <- "mcmullins_south"
LCBDall$site <- factor( LCBDall$site )
LCBDall$group <- factor( LCBDall$group, levels = c("prokaryote","microeukaryote","macroeukaryote"))
LCBD_all <- ggplot( data=LCBDall, aes( x=factor(year), y=site, 
                                       size=LCBD, fill=p.adj )) +
  facet_wrap( ~group) +
  geom_point(shape=21,alpha=1) + #,position=position_jitter(width=0.1,height = 0.1)) +
  ylab( "Site" ) + xlab( "Year" ) + 
  theme_bw( ) + theme(legend.position="top")
LCBD_all
ggsave("R_Code_and_Analysis/betadiversity/LCBD_all.png", 
       plot = LCBD_all, width=8, height=3.5,dpi=300 )
