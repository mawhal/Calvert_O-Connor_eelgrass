####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will investigate beta diversity of grazer (and microbe)
###        community data
###
### started 3 March 2020 by Matt WHalen
###  
####################################################################


library(tidyverse)
library(adespatial)
library(vegan)


## Goal calculate local contributions to beta diversity in different years (do we get the same answer?)

# steps - split data into years, make a list of community matrices, perform calculation
#

### read data
# read data
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
  filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean <- m.sum %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance))


# make a community dataset
m.meta <- m.sum %>% 
  spread( taxon, abundance, fill=0 ) 

# split by year
comm.years <- split( m.meta, m.meta$group )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:6)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:6)],df))
})
lapply(comm.zero, dim)

meta <- lapply(comm.zero, function(z) z[,1:6])
d <- lapply(comm.zero, function(z) z[,-c(1:6)])



# Species and Local contributions to beta diversity
# transform community data
betas <- lapply( d, beta.div, "log.chord", nperm=999 )
signif = lapply( betas, function(z) which(z$p.LCBD <= 0.05) ) # Which are the significant LCBD indices?
nonsignif = lapply( betas, function(z) which(z$p.LCBD > 0.05) ) # Which are the non-significant LCBD indices?

# collate all LCBD info
LCBD <- lapply( betas, function(z) cbind(z$LCBD,z$p.adj) )
LCBD <- data.frame( do.call( rbind, LCBD ) )
names(LCBD) <- c("LCBD","p.adj")
LCBD <- data.frame( m.meta[,1:5], LCBD )

# windows(3.5,4)
LCBD_inverts <- ggplot( data=LCBD, aes( x=factor(year), y=factor(site), size=LCBD, fill=p.adj )) +
  geom_point(shape=21,alpha=0.5,position=position_jitter(width=0.1,height = 0.1)) +
  ylab( "Site" ) + xlab( "Year" ) + ggtitle("Macroeukaryotes") + theme_bw( )
ggsave("R_Code_and_Analysis/betadiversity/LCBD_inverts.png", plot = LCBD_inverts, width=250, height=200, units="mm",dpi=300)

## diversity partitioning
# community matrix, matrix with levels of sampling hierarchy


# start with 2015
d15 <- d[[2]]
m15 <- meta[[2]]
x <- lapply( meta, function(z) data.frame(z[,c(6,4,3)],whole=1) )
x15 <- x[[2]]
adipart( d15, x15, index="richness", weights="prop", relative=F, nsimul = 1000 )
multipart( d15, x15, index="renyi", weights="prop", relative=F, nsimul = 1000 )
adi1 <- adipart( d15, x15, index="shannon", weights="unif", relative=T, nsimul = 1000 )
multi1 <- multipart( d15, x15, index="renyi",  relative=T, nsimul = 1000 )
plot(x=rep(1,4),y=adi1$oecosimu$means[-c(2,3,4)], ylim=c(0,1))
add.plot <- data.frame( level=c("alpha","beta_sample","beta_site","beta_region"), add_rel=adi1$oecosimu$means[-c(2,3,4)] )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
mul.plot <- data.frame( level=c("alpha","beta_sample","beta_site","beta_region"), mul_rel=multi1$oecosimu$means[-c(2,3,4)] )
mul.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=1,fill=level,y=add_rel)) + geom_bar(stat='identity')


adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 1000 ), 
        d, x )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
                d, x )
multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2014:2017),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')
ggsave( "R_Code_and_Analysis/betadiversity/inverts_gamma_partition_additive_taxon3_shannon.png", width=4, height=3 )
add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write_csv(add.plot2,"R_Code_and_Analysis/output_data/adipart_invert_family.csv")
g1 <- ggplot(add.plot2 %>% filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2014:2017) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
