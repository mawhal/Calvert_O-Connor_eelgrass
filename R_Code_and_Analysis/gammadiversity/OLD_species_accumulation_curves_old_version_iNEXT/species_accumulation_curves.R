### Species accumulation curves ###
### Author: Bianca Trevizan Segovia ###
### Date created: May 20th, 2020 ###

library(vegan)
library(tidyverse)

#########################################
############ 16S prokaryotes ############
#########################################

### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", header=T)
names(microbes_16S_ASV)[1:16]

### Creating an object to store abundances only
abundances_16S <- microbes_16S_ASV %>% 
  dplyr::select(-(1:15))

# species accumulation curve for all years
curve_16S = specaccum(abundances_16S, method = "random", gamma = "chao1",
                      permutations = 100)
df.curve_16S <- data.frame(curve_16S$richness,curve_16S$sites,curve_16S$sd)
names(df.curve_16S) <- substring(names(df.curve_16S),11) # keep labels in columns from 16 character on
df.curve_16S$type <- "pooled"

#subset each year
abundances_16S_2015 <- microbes_16S_ASV %>% 
  dplyr::filter(year == "2015") %>% 
  dplyr::select(-(1:15))
abundances_16S_2016 <- microbes_16S_ASV %>% 
  dplyr::filter(year == "2016") %>% 
  dplyr::select(-(1:15))
abundances_16S_2017 <- microbes_16S_ASV %>% 
  dplyr::filter(year == "2017") %>% 
  dplyr::select(-(1:15))
abundances_16S_2018 <- microbes_16S_ASV %>% 
  dplyr::filter(year == "2018") %>% 
  dplyr::select(-(1:15))

# species accumulation curve for each year
curve_16S_2015 = specaccum(abundances_16S_2015, method = "random", gamma = "chao1")
curve_16S_2016 = specaccum(abundances_16S_2016, method = "random", gamma = "chao1")
curve_16S_2017 = specaccum(abundances_16S_2017, method = "random", gamma = "chao1")
curve_16S_2018 = specaccum(abundances_16S_2018, method = "random", gamma = "chao1")

#plot all and each year together
plot(curve_16S$sites, curve_16S$richness,
                       xlab="Number of samples",
                       ylab="Species Richness",
                       main="Prokaryotes")

plot(curve_16S_2015, add = TRUE, col = 2) 
plot(curve_16S_2016, add = TRUE, col = 3)
plot(curve_16S_2017, add = TRUE, col = 4)
plot(curve_16S_2018, add = TRUE, col = 5)


df.curve_16S_2015 <- data.frame(curve_16S_2015$richness,curve_16S_2015$sites,curve_16S_2015$sd)
names(df.curve_16S_2015) <- substring(names(df.curve_16S_2015),16) # keep labels in columns from 16 character on
df.curve_16S_2016 <- data.frame(curve_16S_2016$richness,curve_16S_2016$sites,curve_16S_2016$sd)
names(df.curve_16S_2016) <- substring(names(df.curve_16S_2016),16) 
df.curve_16S_2017 <- data.frame(curve_16S_2017$richness,curve_16S_2017$sites,curve_16S_2017$sd)
names(df.curve_16S_2017) <- substring(names(df.curve_16S_2017),16) 
df.curve_16S_2018 <- data.frame(curve_16S_2018$richness,curve_16S_2018$sites,curve_16S_2018$sd)
names(df.curve_16S_2018) <- substring(names(df.curve_16S_2018),16) 

df.curve_16S_2015$type <- "2015" 
df.curve_16S_2016$type <- "2016" 
df.curve_16S_2017$type <- "2017" 
df.curve_16S_2018$type <- "2018" 

# rbind everything together
df.curve_16S_all <- rbind(df.curve_16S_2015,df.curve_16S_2016,df.curve_16S_2017, df.curve_16S_2018)

curve_16S_all <- ggplot(df.curve_16S_all, aes(x=sites, y=richness))+
  facet_wrap(~type) +
  geom_point() +
  ggtitle("Prokaryotes") +
  labs(y = "Species Richness", x = "Number of samples") 
curve_16S_all

# ggsave("R_Code_and_Analysis/alphadiversity/specaccum_curve_prokaryotes.png", curve_16S_all, width=4, height=4.5,dpi=300)

# repeat with pooled samples
# rbind everything together
df.curve_16S_all <- rbind(df.curve_16S, df.curve_16S_2015,df.curve_16S_2016,df.curve_16S_2017, df.curve_16S_2018)

# smaller df for labelling years
labels_16S <- df.curve_16S_all %>% 
  group_by( type ) %>% 
  filter( sites==max(sites), type != "pooled" )

# library(rcompanion)
# groupwiseMean(richness ~ type, 
#               data   = df.curve_16S_all, 
#               conf   = 0.95, 
#               digits = 3)

df.curve_16S_all$ymin <- df.curve_16S_all$richness - 1.96*df.curve_16S_all$sd/sqrt(1)

df.curve_16S_all$ymax <- df.curve_16S_all$richness + 1.96*df.curve_16S_all$sd/sqrt(1)

curve_16S_all <- ggplot(df.curve_16S_all, aes(x=sites, y=richness, colour=type))+
  # facet_wrap(~type) +
  geom_line(aes(group=type), size=1.5) +
  geom_ribbon(aes( ymin=ymin, ymax=ymax, fill=type), alpha = 0.3)+
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")+
  geom_text( data=labels_16S,aes(x=sites,y=richness,label=type), adj=0, 
             nudge_y = c(0,0,40,0), nudge_x = c(2,2,1,2) ) +
  # ggtitle("Prokaryotes") +
  labs(y = "Species richness", x = "Number of samples", title = "Prokaryotes") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    panel.grid.major = element_blank(), #remove major grid
    panel.grid.minor = element_blank(), #remove minor grid
    axis.line = element_line(colour = "black"), #draw line in the axis)
    legend.position ="none",
    plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0)),
    panel.border = element_blank()) #remove lines outside the graph
curve_16S_all

# extrapolate
poolaccum_16S <- poolaccum(abundances_16S)
plot_poolaccum_16S <- plot(poolaccum_16S)
plot_poolaccum_16S
estimateR(abundances_16S)
 
# mean and range of chao estimator
meanchao <- apply( poolaccum_16S$chao, 1, mean)
rangechao <-  t(apply( poolaccum_16S$chao, 1, quantile, c(0.025,0.9755)))
chao <- as.data.frame(poolaccum_16S$chao)
chao$sample = 1:nrow(chao)
chao$mean <- meanchao
chao$min <- rangechao[,1]
chao$max <- rangechao[,2]

ggplot( chao, aes(y=mean, x=sample) ) + geom_line()+
  geom_line( aes(y=max) ) + 
  geom_line( aes(y=min) ) 

#############################################
############ 18S microeukaryotes ############
#############################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", header=T)
names(microbes_18S_ASV)[1:18]

### Creating an object to store abundances only
abundances_18S <- microbes_18S_ASV %>% 
  dplyr::select(-(1:15))

curve_18S = specaccum(abundances_18S, method = "random", gamma = "chao1",
                      permutations = 100)
df.curve_18S <- data.frame(curve_18S$richness,curve_18S$sites,curve_18S$sd)
names(df.curve_18S) <- substring(names(df.curve_18S),11) # keep labels in columns from 16 character on
df.curve_18S$type <- "pooled"

#subset each year
abundances_18S_2015 <- microbes_18S_ASV %>% 
  dplyr::filter(year == "2015") %>% 
  dplyr::select(-(1:15))
abundances_18S_2016 <- microbes_18S_ASV %>% 
  dplyr::filter(year == "2016") %>% 
  dplyr::select(-(1:15))
abundances_18S_2017 <- microbes_18S_ASV %>% 
  dplyr::filter(year == "2017") %>% 
  dplyr::select(-(1:15))
abundances_18S_2018 <- microbes_18S_ASV %>% 
  dplyr::filter(year == "2018") %>% 
  dplyr::select(-(1:15))

# species accumulation curve for each year
curve_18S_2015 = specaccum(abundances_18S_2015, method = "random", gamma = "chao1")
curve_18S_2016 = specaccum(abundances_18S_2016, method = "random", gamma = "chao1")
curve_18S_2017 = specaccum(abundances_18S_2017, method = "random", gamma = "chao1")
curve_18S_2018 = specaccum(abundances_18S_2018, method = "random", gamma = "chao1")

#plot all and each year together
plot(curve_18S$sites, curve_18S$richness,
     xlab="Number of samples",
     ylab="Species Richness",
     main="Prokaryotes")

plot(curve_18S_2015, add = TRUE, col = 2) 
plot(curve_18S_2016, add = TRUE, col = 3)
plot(curve_18S_2017, add = TRUE, col = 4)
plot(curve_18S_2018, add = TRUE, col = 5)


df.curve_18S_2015 <- data.frame(curve_18S_2015$richness,curve_18S_2015$sites,curve_18S_2015$sd)
names(df.curve_18S_2015) <- substring(names(df.curve_18S_2015),16) # keep labels in columns from 16 character on
df.curve_18S_2016 <- data.frame(curve_18S_2016$richness,curve_18S_2016$sites,curve_18S_2016$sd)
names(df.curve_18S_2016) <- substring(names(df.curve_18S_2016),16) 
df.curve_18S_2017 <- data.frame(curve_18S_2017$richness,curve_18S_2017$sites,curve_18S_2017$sd)
names(df.curve_18S_2017) <- substring(names(df.curve_18S_2017),16) 
df.curve_18S_2018 <- data.frame(curve_18S_2018$richness,curve_18S_2018$sites,curve_18S_2018$sd)
names(df.curve_18S_2018) <- substring(names(df.curve_18S_2018),16) 

df.curve_18S_2015$type <- "2015" 
df.curve_18S_2016$type <- "2016" 
df.curve_18S_2017$type <- "2017" 
df.curve_18S_2018$type <- "2018" 

# rbind everything together
df.curve_18S_all <- rbind(df.curve_18S_2015,df.curve_18S_2016,df.curve_18S_2017, df.curve_18S_2018)

curve_18S_all <- ggplot(df.curve_18S_all, aes(x=sites, y=richness, colour=type))+
  facet_wrap(~type) +
  geom_point() +
  ggtitle("Microeukaryotes") +
  labs(y = "Species Richness", x = "Number of samples") 
curve_18S_all

# ggsave("R_Code_and_Analysis/alphadiversity/specaccum_curve_microeukaryotes.png", curve_18S_all, width=250, height=200, units="mm",dpi=300)


## repeat with pooled samples
# rbind everything together
df.curve_18S_all <- rbind(df.curve_18S, df.curve_18S_2015,df.curve_18S_2016,df.curve_18S_2017, df.curve_18S_2018)

# smaller df for labelling years
labels_18S <- df.curve_18S_all %>% 
  group_by( type ) %>% 
  filter( sites==max(sites), type != "pooled" )

curve_18S_all <- ggplot(df.curve_18S_all, aes(x=sites, y=richness, colour=type))+
  # facet_wrap(~type) +
  geom_line(aes(group=type), size=1.5) +
  scale_colour_brewer(palette = "Dark2") +
  geom_text( data=labels_18S,aes(x=sites,y=richness,label=type), adj=0, 
             nudge_y = c(5,5,5,5), nudge_x = c(2,2,2,2) ) +
  labs(y = "", x = "Number of samples", title= "Microeukaryotes")  +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    panel.grid.major = element_blank(), #remove major grid
    panel.grid.minor = element_blank(), #remove minor grid
    axis.line = element_line(colour = "black"), #draw line in the axis)
    legend.position ="none",
    plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0)),
    panel.border = element_blank()) #remove lines outside the graph
curve_18S_all

# extrapolate
poolaccum_18S <- poolaccum(abundances_18S)
plot_poolaccum_18S <- plot(poolaccum_18S)
plot_poolaccum_18S



#################################################
############ Inverts macroeukaryotes ############
#################################################

#### Finest level ####

### Read table metadata and abundances
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

# filter taxa and sites # taxon 2 = finest
mfilt_finest <- m %>%
  select( year, group, region, site, sample, ID, taxon = taxon2, remove, size ) %>% 
  filter( is.na(remove), !is.na(taxon))

# summarize taxon counts per sample
m.sum_finest <- mfilt_finest %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, sample, ID, taxon ) %>% 
  summarize( abundance=length(size) )

# summarinze mean abundance per site
m.mean_finest <- m.sum_finest %>% 
  # unite( "ID", year,site,sample, remove=FALSE ) %>% 
  group_by( year, group, region, site, taxon ) %>% 
  summarize( abundance=mean(abundance)) 

# make a community dataset
m.meta_finest <- m.sum_finest %>% 
  pivot_wider( names_from = taxon, values_from = abundance, values_fill=0 ) # had to replace spread by pivot_wider (new tidyr) ### not working anymore

names(m.meta_finest)

m.meta_finest <- m.meta_finest %>% 
  ungroup( year, group, region, site )

abundances_inverts_finest <- m.meta_finest[,7:ncol(m.meta_finest)]

curve_inverts = specaccum(abundances_inverts_finest, method = "random", gamma = "chao1",
                      permutations = 100)
df.curve_inverts <- data.frame(curve_inverts$richness,curve_inverts$sites,curve_inverts$sd)
names(df.curve_inverts) <- substring(names(df.curve_inverts),15) # keep labels in columns from 16 character on
df.curve_inverts$type <- "pooled"

#subset each year
abundances_inverts_2014 <- m.meta_finest %>% 
  dplyr::filter(year == "2014")
abundances_inverts_2014 <- abundances_inverts_2014[,15:ncol(abundances_inverts_2014)]

abundances_inverts_2015 <- m.meta_finest %>% 
  dplyr::filter(year == "2015")
abundances_inverts_2015 <- abundances_inverts_2015[,15:ncol(abundances_inverts_2015)]

abundances_inverts_2016 <- m.meta_finest %>% 
  dplyr::filter(year == "2016")
abundances_inverts_2016 <- abundances_inverts_2016[,15:ncol(abundances_inverts_2016)]

abundances_inverts_2017 <- m.meta_finest %>% 
  dplyr::filter(year == "2017")
abundances_inverts_2017 <- abundances_inverts_2017[,15:ncol(abundances_inverts_2017)]

# species accumulation curve for each year
curve_inverts_2014 = specaccum(abundances_inverts_2014, method = "random", gamma = "chao1")
curve_inverts_2015 = specaccum(abundances_inverts_2015, method = "random", gamma = "chao1")
curve_inverts_2016 = specaccum(abundances_inverts_2016, method = "random", gamma = "chao1")
curve_inverts_2017 = specaccum(abundances_inverts_2017, method = "random", gamma = "chao1")

#plot all and each year together
plot(curve_inverts$sites, curve_inverts$richness,
     xlab="Number of samples",
     ylab="Species Richness",
     main="Prokaryotes")

plot(curve_inverts_2014, add = TRUE, col = 2) 
plot(curve_inverts_2015, add = TRUE, col = 3)
plot(curve_inverts_2016, add = TRUE, col = 4)
plot(curve_inverts_2017, add = TRUE, col = 5)

df.curve_inverts_2014 <- data.frame(curve_inverts_2014$richness,curve_inverts_2014$sites,curve_inverts_2014$sd)
names(df.curve_inverts_2014) <- substring(names(df.curve_inverts_2014),20) 
df.curve_inverts_2015 <- data.frame(curve_inverts_2015$richness,curve_inverts_2015$sites,curve_inverts_2015$sd)
names(df.curve_inverts_2015) <- substring(names(df.curve_inverts_2015),20) 
df.curve_inverts_2016 <- data.frame(curve_inverts_2016$richness,curve_inverts_2016$sites,curve_inverts_2016$sd)
names(df.curve_inverts_2016) <- substring(names(df.curve_inverts_2016),20) 
df.curve_inverts_2017 <- data.frame(curve_inverts_2017$richness,curve_inverts_2017$sites,curve_inverts_2017$sd)
names(df.curve_inverts_2017) <- substring(names(df.curve_inverts_2017),20) 

df.curve_inverts_2014$type <- "2014" 
df.curve_inverts_2015$type <- "2015" 
df.curve_inverts_2016$type <- "2016" 
df.curve_inverts_2017$type <- "2017" 

# rbind everything together
df.curve_inverts_all <- rbind(df.curve_inverts_2014, df.curve_inverts_2015,df.curve_inverts_2016,df.curve_inverts_2017)

curve_inverts_all <- ggplot(df.curve_inverts_all, aes(x=sites, y=richness))+
  facet_wrap(~type) +
  geom_point() +
  ggtitle("Macroeukaryotes") +
  labs(y = "Species Richness", x = "Number of samples") 
curve_inverts_all

# ggsave("R_Code_and_Analysis/alphadiversity/specaccum_curve_macroeukaryotes.png", curve_inverts_all, width=250, height=200, units="mm",dpi=300)

## repeat with pooled samples
# rbind everything together
df.curve_inverts_all <- rbind(df.curve_inverts, df.curve_inverts_2014, df.curve_inverts_2015,df.curve_inverts_2016,df.curve_inverts_2017)

# smaller df for labelling years
labels_inverts <- df.curve_inverts_all %>% 
  group_by( type ) %>% 
  filter( sites==max(sites), type != "pooled" )

curve_inverts_all <- ggplot(df.curve_inverts_all, aes(x=sites, y=richness, colour=type))+
  # facet_wrap(~type) +
  geom_line(aes(group=type),size=1.5) +
  scale_colour_manual(values=c("#E6AB02", "#1B9E77", "#D95F02", "#7570B3", "#66A61E")) +
  geom_text( data=labels_inverts,aes(x=sites,y=richness,label=type), adj=0, 
             nudge_y = c(0,0,1,2), nudge_x = c(2,2,2,-7) ) +
  labs(y = "", x = "Number of samples", title = "Macroeukaryotes") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    panel.grid.major = element_blank(), #remove major grid
    panel.grid.minor = element_blank(), #remove minor grid
    axis.line = element_line(colour = "black"), #draw line in the axis)
    legend.position ="none",
    plot.title = element_text(size=20, face = "bold", hjust = 0.1, margin = margin(t = 0, r =0, b = 20, l = 0)),
    panel.border = element_blank()) #remove lines outside the graph
curve_inverts_all

# extrapolate
poolaccum_inverts <- poolaccum(abundances_inverts_finest)
plot_poolaccum_inverts <- plot(poolaccum_inverts)
plot_poolaccum_inverts
library(iNEXT)
incidence_inverts_finest <- as.data.frame(abundances_inverts_finest)
incidence_inverts_finest <- ifelse( incidence_inverts_finest==0,0,1)
iNEXT(as.matrix(incidence_inverts_finest), datatype = "incidence_raw" )


# plot all curves together
cowplot::plot_grid(curve_16S_all,curve_18S_all,curve_inverts_all, labels = "AUTO", ncol = 3, label_x =.05, hjust = 1, label_size=20)
ggsave( "R_Code_and_Analysis/alphadiversity/accumulation_curves_all.png", 
        width=30, height=10 )


# all extrapolations
poolaccum_16S
poolaccum_18S
apply(poolaccum_inverts$chao, 1, mean)
