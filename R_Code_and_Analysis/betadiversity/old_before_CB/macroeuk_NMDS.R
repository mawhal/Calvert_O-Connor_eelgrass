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
m <- read_csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv" )
# replace spaces with periods for consistency and merging names
m$taxon <- gsub( " ", ".", m$taxon )
# bring in the updated data 
mtaxa_update <- read_csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )
# merge
m <- left_join( m,  mtaxa_update )
# replace periods with spaces for making simple names in vegan

### select sites
# filter taxa and sites
mfilt <- m %>%
  filter( !is.na(taxon2) )


# summarize taxon counts per sample
m.sum.fine <- mfilt %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, site, sample, taxon=taxon2 ) %>% 
  summarize( abundance=length(size) )
m.spread.fine <- m.sum.fine %>% 
  spread( taxon, abundance, fill=0 )
meta.fine <- data.frame(m.spread.fine[,c(1:3)])
comm.fine <- data.frame(m.spread.fine[,-c(1:3)])
names(comm.fine) <- make.cepnames( names(comm.fine) )
m.sum.family <- mfilt %>% 
  group_by( year, site, sample, taxon=taxon4 ) %>% 
  summarize( abundance=length(size) )
m.spread.family <- m.sum.family %>% 
  spread( taxon, abundance, fill=0 )
meta.family <- data.frame(m.spread.family[,c(1:3)])
comm.family <- data.frame(m.spread.family[,-c(1:3)])
comm.family <- comm.family[,names(comm.family)!="X.NA."]
names(comm.family) <- make.cepnames( names(comm.family) )
m.sum.coarse <- mfilt %>% 
  group_by( year, site, sample, taxon=taxon5 ) %>% 
  summarize( abundance=length(size) )
m.spread.coarse <- m.sum.coarse %>% 
  spread( taxon, abundance, fill=0 )
meta.coarse <- data.frame(m.spread.coarse[,c(1:3)])
comm.coarse <- data.frame(m.spread.coarse[,-c(1:3)])
comm.coarse <- comm.coarse[,names(comm.coarse)!="X.NA."]
names(comm.coarse) <- make.cepnames( names(comm.coarse) )

names(comm) <- make.cepnames( names(comm) )


# - first save the bray-curtis distances
meta1 <- meta %>%
  filter(year==2015) %>%
  unite(sample, site, sample, sep="_")
sample.names <- make.cepnames(meta1$sample)
write_csv( meta1, "output_data/2015_grazer_metadata.csv")
commdist <- vegdist( comm[meta$year==2015,], method = "bray" )
commdist <- as.matrix(commdist)
rownames(commdist) <- sample.names
colnames(commdist) <- sample.names

write_csv( data.frame(commdist), "output_data/2015_grazer_braycurtis.csv")




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


# nmds of comparable data


# NMDS
mds <- metaMDS( comm, distance="bray", k=7 )
mds # very high stress with 2 axes 

# extract points for first two axes
meta <- cbind(meta, mds$points[,c(1,2)])

# convex hull
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

meta$site2 <- factor( meta$site, ordered=T )
man.col <- viridis::viridis(5)
cols <- c( 'black','black',man.col[2],man.col[2],man.col[2],
           man.col[3],man.col[3],man.col[4],man.col[4],man.col[5] )
ggplot( meta, aes(MDS1,MDS2,col=site)) + 
  # stat_chull( fill=NA ) +
  stat_ellipse( aes(lty=site), lwd=1 ) +
  geom_point( size=3 ) +
  scale_linetype_manual( values=c(1,2,1,2,3,1,2,1,2,1)) +
  scale_color_manual( values=cols  )





## figure size parameters
h <- 4
w <- 5
dpi <- 300

#### constrained ordination - rdRDA/capscale
rda1 <- capscale( comm ~ as.factor(year)+site, data=meta, distance="bray" )
plot(rda1)
summary(rda1)
centroids <- as.data.frame(summary(rda1)$centroids)
centroids$marker <- 2
centroids$region   <- substr(unlist(lapply( strsplit( rownames(centroids), split="_" ), function(z) z[1] )),5,20)
centroids$region[1:4] <- "year"
centroids$marker[1:4] <- 1
pointsrda <- as.data.frame(scores(rda1,scaling = 2)$sites)
plot(pointsrda, pch=19, cex=0.5)
  points( centroids, cex=centroids$marker, 
          col=as.numeric(as.factor(centroids$region)),
          pch=centroids$marker+1 )
pointsrda <- bind_cols(pointsrda,meta)
pointsrda$region <- unlist( lapply( strsplit(as.character(pointsrda$site),split="_"),function(z) z[1] ) )
pointsrda$region <- factor(pointsrda$region,
                           levels=c("mcmullins","goose","triquet","choked","pruth"),ordered=T)
sitetroids <- centroids[-c(1:4),]
sitetroids$region <- factor(sitetroids$region,
                           levels=c("mcmullins","goose","triquet","choked","pruth"),ordered=T)
yeartroids <- centroids[1:4,]
yeartroids$year   <- 2014:2017
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
rda2 <- capscale( comm ~ 1, data=meta, distance="bray" )
plot(rda2)
summary(rda2)
pointsrda <- as.data.frame(scores(rda2,scaling = 2)$sites)
plot(pointsrda, pch=19, cex=0.5)

meta$region <- unlist( lapply( strsplit(as.character(meta$site),split="_"),function(z) z[1] ) )
meta$region <- factor( meta$region,
                           levels=c("mcmullins","goose","triquet","choked","pruth"),ordered=T)

pointsrda <- bind_cols(pointsrda,meta)
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

# match nmds to capscale plots
# get group centroids
sitetroids <- meta %>% 
  group_by(site,region) %>% 
  summarize( MDS1=mean(MDS1),MDS2=mean(MDS2))
yeartroids <- meta %>% 
  group_by(year) %>% 
  summarize( MDS1=mean(MDS1),MDS2=mean(MDS2))
ggplot( meta, aes(MDS1,MDS2)) + 
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