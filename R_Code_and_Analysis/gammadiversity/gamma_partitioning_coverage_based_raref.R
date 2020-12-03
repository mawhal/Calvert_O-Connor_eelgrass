### ADIPART ###
### Author: Bianca Trevizan Segovia, modified using Matt's code ###
### Date created: December 3rd, 2020 ###

library(vegan)
library(tidyverse)
library(dplyr)

#############################
############ 16S ############
#############################

## ASV level ###

### Read table metadata and abundances
### Read table metadata and abundances
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_ASV)[1:16]

# split by year
comm.years <- split( microbes_16S_ASV, microbes_16S_ASV$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:15)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:15)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:15])
d <- lapply(comm.zero, function(z) z[,-c(1:15)])

x <- lapply( meta, function(z) data.frame(z[,c("site_quadrat_id", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

# multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
#                   d, x )
# multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
# multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
# lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2015:2018),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="R_Code_and_Analysis/output_data/adipart_16S_ASV_level.csv", row.names = F)
g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2015:2018) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/gammadiversity/microbes_16S_ASV_gamma_partition_additive_richness_prop.png", width=4, height=3 )


## FAMILY LEVEL adipart to compare ###
### Read table metadata and abundances
microbes_16S_family <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_family_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_family)[1:17]

# split by year
comm.years_family <- split( microbes_16S_family, microbes_16S_family$year )
lapply(comm.years_family, dim)

# remove taxa that did not appear in a given year
comm.zero_family <- lapply( comm.years_family, function(z) {
  tmp <- z[,-c(1:12)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:12)],df))
})

meta_family <- lapply(comm.zero_family, function(z) z[,1:12])
d_family <- lapply(comm.zero_family, function(z) z[,-c(1:12)])

x_family <- lapply( meta_family, function(z) data.frame(z[,c("site_quadrat_id", "site", "region")],whole=1) )

adds_family <- mapply( function(comm,x_family) adipart(comm,x_family, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d_family, x_family )
adds.mean_family <- lapply( adds_family[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs_family <- lapply(adds_family[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds_family[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

add.plot_family <- data.frame( year=rep(c(2015:2018),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs_family ), null=do.call( c, adds.mean_family ) )
add.plot_family$level <- factor( add.plot_family$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot_family, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2_family <- add.plot_family %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2_family, file="R_Code_and_Analysis/output_data/adipart_16S_family_level.csv", row.names = F)


# #############################
############ 18S ############
#############################

### Read table metadata and abundances
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_18S_ASV)[1:13]

# split by year
comm.years <- split( microbes_18S_ASV, microbes_18S_ASV$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:9)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:9)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:9])
d <- lapply(comm.zero, function(z) z[,-c(1:9)])

x <- lapply( meta, function(z) data.frame(z[,c("site_quadrat_id", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

# multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
#                   d, x )
# multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
# multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
# lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2015:2018),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="R_Code_and_Analysis/output_data/adipart_18S_ASV_level.csv", row.names=F)

g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2015:2018) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/gammadiversity/microbes_18S_ASV_gamma_partition_additive_richness_prop.png", width=4, height=3 )

# sampling design (alpha, beta1= within sites, b2=between sites, b3=between regions)

## FAMILY LEVEL adipart to compare ###
### Read table metadata and abundances
microbes_18S_family <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_family_level_1000_COVERAGE_RAREF.csv", header=T)
nrow(microbes_18S_family)
microbes_18S_family <- microbes_18S_family[rowSums(microbes_18S_family[-(1:9)]) !=0, ]
names(microbes_18S_family)[1:13]
# split by year
comm.years_family <- split( microbes_18S_family, microbes_18S_family$year )
lapply(comm.years_family, dim)

# remove taxa that did not appear in a given year
comm.zero_family <- lapply( comm.years_family, function(z) {
  tmp <- z[,-c(1:9)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:9)],df))
})

meta_family <- lapply(comm.zero_family, function(z) z[,1:9])
d_family <- lapply(comm.zero_family, function(z) z[,-c(1:9)])

x_family <- lapply( meta_family, function(z) data.frame(z[,c("site_quadrat_id", "site", "region")],whole=1) )

adds_family <- mapply( function(comm,x_family) adipart(comm,x_family, index="richness",weights="prop", 
                                                       relative=T, nsimul = 20 ), 
                       d_family, x_family )
adds.mean_family <- lapply( adds_family[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs_family <- lapply(adds_family[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds_family[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

add.plot_family <- data.frame( year=rep(c(2015:2018),each=4),
                               level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                               observed=do.call( c, adds.obs_family ), null=do.call( c, adds.mean_family ) )
add.plot_family$level <- factor( add.plot_family$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot_family, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2_family <- add.plot_family %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2_family, file="R_Code_and_Analysis/output_data/adipart_18S_family_level.csv", row.names = F)

#################################################
############ Inverts macroeukaryotes ############
#################################################

#### Finest level ####

### Read table metadata and abundances
inverts_finest <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")
names(inverts_finest)[1:12]

# split by year
comm.years <- split( inverts_finest, inverts_finest$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:7)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:7)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:7])
d <- lapply(comm.zero, function(z) z[,-c(1:7)])

x <- lapply( meta, function(z) data.frame(z[,c("ID", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

# multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
#                   d, x )
# multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
# multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
# lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2014:2017),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="R_Code_and_Analysis/output_data/adipart_inverts_finest_level.csv", row.names=F)

g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2014:2017) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/gammadiversity/inverts_finest_level_gamma_partition_additive_richness_prop.png", width=4, height=3 )

#### Family level ####

inverts_family <- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_family_1000_COVERAGE_RAREF.csv")
names(inverts_family)[1:12]
inverts_family <- inverts_family %>% 
  drop_na()  #remove samples with NAs (why did they drop out with CB rarefaction?)

# split by year
comm.years <- split( inverts_family, inverts_family$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:7)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:7)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:7])
d <- lapply(comm.zero, function(z) z[,-c(1:7)])

x <- lapply( meta, function(z) data.frame(z[,c("ID", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

# multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
#                   d, x )
# multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
# multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
# lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2014:2017),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="R_Code_and_Analysis/output_data/adipart_inverts_family_level.csv", row.names=F)

g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2014:2017) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/gammadiversity/inverts_family_level_gamma_partition_additive_richness_prop.png", width=4, height=3 )
