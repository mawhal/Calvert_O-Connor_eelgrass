####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will compare epifaunal communities in basic ways
###      by comparing summary patterns of richness, diversity, composition
###  code by Matt Whalen
###  started on   21 November 2019
###  
####################################################################



# load libraries
library( tidyverse )
library( vegan )
library( cowplot )

# read data
m <- read_csv( "output_data/O'Connor_hakai_seagrass_MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "output_data/O'Connor_hakai_seagrass_taxa_edit_20191114.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan

### select sites
# filter taxa and sites
mfilt <- m %>%
  filter( is.na(remove), !is.na(taxon4))
# # Four core sites
# mfilt <- m %>%
#   filter( is.na(remove), !is.na(taxon4),
#           site %in% c("inner choked","sandspit",
#                       "triquet north","triquet south") )


# summarize taxon counts per sample
m.sum <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon4 ) %>% 
  summarize( abundance=length(size) )



# make a community dataset
m.meta <- m.sum %>% 
  spread( taxon4, abundance, fill=0 )

meta <- data.frame(m.meta[,c(1:3)])
# meta$year <- factor( meta$year, ordered=T )
comm <- data.frame(m.meta[,-c(1:3)])
names(comm) <- make.cepnames( names(comm) )

# richness and diversity
meta$richness <- specnumber(comm)
meta$shannon <- diversity(comm)

# accumulation+collector curves
plot(specaccum( comm ))


# add geolocation
geo <- read_csv( "../metadata/coreen_hakai_sites_location_2015.csv" )
geo <- geo %>% 
  select( site, lat, long )

meta <- left_join( meta, geo )
# rename and reorder sites
meta$site <- gsub("inner choked","choked_inner",meta$site)
meta$site <- gsub("sandspit","choked_sandspit",meta$site)
meta$site <- gsub(" ","_",meta$site)
meta$site <- fct_reorder( meta$site, -meta$lat )


# comparable data across years
with(meta, table(year,site))
meta1415 <- meta %>% filter(year %in% 2014:2015)
comm1415 <- comm[meta$year %in% 2014:2015,]

# pick which community data to use
usecomm <- comm
usemeta <- meta

## figure size parameters
h <- 4
w <- 5
dpi <- 300

#### constrained ordination - rdRDA/capscale
rda1 <- capscale( usecomm ~ as.factor(year)+site, data=usemeta, distance="bray" )
plot(rda1)
summary(rda1)
centroids <- as.data.frame(summary(rda1)$centroids)
centroids$marker <- 2
centroids$region   <- substr(unlist(lapply( strsplit( rownames(centroids), split="_" ), function(z) z[1] )),5,20)
centroids$region[1:length(unique(usemeta$year))] <- "year"
centroids$marker[1:length(unique(usemeta$year))] <- 1
pointsrda <- as.data.frame(scores(rda1,scaling = 2)$sites)
plot(pointsrda, pch=19, cex=0.5)
  points( centroids, cex=centroids$marker, 
          col=as.numeric(as.factor(centroids$region)),
          pch=centroids$marker+1 )
pointsrda <- bind_cols(pointsrda,usemeta)
pointsrda$region <- unlist( lapply( strsplit(as.character(pointsrda$site),split="_"),function(z) z[1] ) )
pointsrda$region <- factor(pointsrda$region,
                           levels=c("mcmullins","goose","triquet","choked","pruth"),ordered=T)
sitetroids <- centroids[-c(1:length(unique(usemeta$year))),]
sitetroids$region <- factor(sitetroids$region,
                           levels=c("mcmullins","goose","triquet","choked","pruth"),ordered=T)
yeartroids <- centroids[1:length(unique(usemeta$year)),]
yeartroids$year   <- unique(usemeta$year)
R2 <- eigenvals(rda1)/sum(eigenvals(rda1))
ggplot( data=pointsrda, aes(x=CAP1,y=CAP2)) +
  stat_ellipse( aes(col=region), level=0.7 ) +
  geom_point(data=sitetroids, aes(fill=region), size=3,pch=21 )+
  geom_point( aes(col=region), alpha=0.5) +
  geom_text(data=yeartroids,aes(label=year),size=4,pch=21) +
  xlab(paste0('dbRDA1 (',round(R2[1],3)*100, '%)')) +
  ylab(paste0('dbRDA2 (',round(R2[2],3)*100, '%)')) +
  viridis::scale_color_viridis(discrete=T) +
  viridis::scale_fill_viridis(discrete=T) 
ggsave( "figs/epifauna_CAP_year+site.pdf",width = w, height = h, dpi = dpi )

#
#### unconstrained ordination - rdRDA/capscale
rda2 <- capscale( usecomm ~ 1, data=usemeta, distance="bray" )
plot(rda2)
summary(rda2)
pointsrda <- as.data.frame(scores(rda2,scaling = 2)$sites)
plot(pointsrda, pch=19, cex=0.5)

usemeta$region <- unlist( lapply( strsplit(as.character(usemeta$site),split="_"),function(z) z[1] ) )
usemeta$region <- factor( usemeta$region,
                           levels=c("mcmullins","goose","triquet","choked","pruth"),ordered=T)

pointsrda <- bind_cols(pointsrda,usemeta)
# get group centroids
sitetroids <- pointsrda %>% 
  group_by(site,region) %>% 
  summarize( MDS1=mean(MDS1),MDS2=mean(MDS2))
yeartroids <- pointsrda %>% 
  group_by(year) %>% 
  summarize( MDS1=mean(MDS1),MDS2=mean(MDS2))
R2 <- eigenvals(rda2)/sum(eigenvals(rda2))
ggplot( data=pointsrda, aes(x=MDS1,y=MDS2)) +
  stat_ellipse( aes(col=region), level=0.7 ) +
  geom_point(data=sitetroids, aes(fill=region), size=3,pch=21 )+
  geom_point( aes(col=region), alpha=0.5) +
  geom_text(data=yeartroids,aes(label=year),size=4,pch=21) +
  xlab(paste0('MDS1 (',round(R2[1],3)*100, '%)')) +
  ylab(paste0('MDS2 (',round(R2[2],3)*100, '%)')) +
  viridis::scale_color_viridis(discrete=T) +
  viridis::scale_fill_viridis(discrete=T) 
ggsave( "figs/epifauna_CAP_1.pdf",width = w, height = h, dpi = dpi )



# NMDS
mds <- metaMDS( usecomm, distance="bray", k=7 )
mds # very high stress with 2 axes 

# extract points for first two axes
usemeta <- cbind(usemeta, mds$points[,c(1,2)])

# match nmds to capscale plots
# get group centroids
sitetroids <- usemeta %>% 
  group_by(site,region) %>% 
  summarize( MDS1=mean(MDS1),MDS2=mean(MDS2))
yeartroids <- usemeta %>% 
  group_by(year) %>% 
  summarize( MDS1=mean(MDS1),MDS2=mean(MDS2))
ggplot( usemeta, aes(MDS1,MDS2)) + 
  stat_ellipse( aes(col=region), level=0.7 ) +
  geom_point(data=sitetroids, aes(fill=region), size=3,pch=21 )+
  geom_point( aes(col=region), alpha=0.5) + 
  geom_text(data=yeartroids,aes(label=year),size=4,pch=21) +
  xlab('MDS1') +
  ylab('MDS2') +
  viridis::scale_color_viridis(discrete=T) +
  viridis::scale_fill_viridis(discrete=T)

ggsave( "figs/epifauna_NMDS_year+site.pdf",width = w, height = h, dpi = dpi )
#