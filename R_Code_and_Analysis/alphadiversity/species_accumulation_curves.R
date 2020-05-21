### Species accumulation curves ###
### Author: Bianca Trevizan Segovia ###
### Date created: May 20th, 2020 ###

library(vegan)
library(dplyr)
library(ggplot2)

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
curve_16S = specaccum(abundances_16S, method = "random", 
                      permutations = 100)
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
curve_16S_2015 = specaccum(abundances_16S_2015, method = "random")
curve_16S_2016 = specaccum(abundances_16S_2016, method = "random")
curve_16S_2017 = specaccum(abundances_16S_2017, method = "random")
curve_16S_2018 = specaccum(abundances_16S_2018, method = "random")

#plot all and each year together
plot(curve_16S$sites, curve_16S$richness,
                       xlab="Number of Sites",
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
  labs(y = "Species Richness", x = "Number of Sites") 
curve_16S_all

ggsave("R_Code_and_Analysis/figs/specaccum_curve_prokaryotes.png", curve_16S_all, width=250, height=200, units="mm",dpi=300)

poolaccum_16S <- poolaccum(abundances_16S)
plot_poolaccum_16S <- plot(poolaccum_16S)
plot_poolaccum_16S

#############################################
############ 18S microeukaryotes ############
#############################################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", header=T)
names(microbes_18S_ASV)[1:18]

### Creating an object to store abundances only
abundances_18S <- microbes_18S_ASV %>% 
  dplyr::select(-(1:15))

curve_18S = specaccum(abundances_18S, method = "random", 
                      permutations = 100)
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
curve_18S_2015 = specaccum(abundances_18S_2015, method = "random")
curve_18S_2016 = specaccum(abundances_18S_2016, method = "random")
curve_18S_2017 = specaccum(abundances_18S_2017, method = "random")
curve_18S_2018 = specaccum(abundances_18S_2018, method = "random")

#plot all and each year together
plot(curve_18S$sites, curve_18S$richness,
     xlab="Number of Sites",
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

curve_18S_all <- ggplot(df.curve_18S_all, aes(x=sites, y=richness))+
  facet_wrap(~type) +
  geom_point() +
  ggtitle("Microeukaryotes") +
  labs(y = "Species Richness", x = "Number of Sites") 
curve_18S_all

ggsave("R_Code_and_Analysis/figs/specaccum_curve_microeukaryotes.png", curve_18S_all, width=250, height=200, units="mm",dpi=300)

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
  spread( taxon, abundance, fill=0 ) 

names(m.meta_finest)[1:16]

m.meta_finest <- m.meta_finest %>% 
  ungroup( year, group, region, site, taxon )

### Creating an object to store abundances only
abundances_inverts_finest <- m.meta_finest %>% 
  dplyr::select(-(1:6))

curve_inverts = specaccum(abundances_inverts_finest, method = "random", 
                      permutations = 100)
#subset each year
abundances_inverts_2014 <- m.meta_finest %>% 
  dplyr::filter(year == "2014") %>% 
  dplyr::select(-(1:15))
abundances_inverts_2015 <- m.meta_finest %>% 
  dplyr::filter(year == "2015") %>% 
  dplyr::select(-(1:15))
abundances_inverts_2016 <- m.meta_finest %>% 
  dplyr::filter(year == "2016") %>% 
  dplyr::select(-(1:15))
abundances_inverts_2017 <- m.meta_finest %>% 
  dplyr::filter(year == "2017") %>% 
  dplyr::select(-(1:15))

# species accumulation curve for each year
curve_inverts_2014 = specaccum(abundances_inverts_2014, method = "random")
curve_inverts_2015 = specaccum(abundances_inverts_2015, method = "random")
curve_inverts_2016 = specaccum(abundances_inverts_2016, method = "random")
curve_inverts_2017 = specaccum(abundances_inverts_2017, method = "random")

#plot all and each year together
plot(curve_inverts$sites, curve_inverts$richness,
     xlab="Number of Sites",
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
  labs(y = "Species Richness", x = "Number of Sites") 
curve_inverts_all

ggsave("R_Code_and_Analysis/figs/specaccum_curve_macroeukaryotes.png", curve_inverts_all, width=250, height=200, units="mm",dpi=300)

poolaccum_inverts <- poolaccum(abundances_inverts_finest)
plot_poolaccum_inverts <- plot(poolaccum_inverts)
plot_poolaccum_inverts
