### ADIPART ###
### Author: Bianca Trevizan Segovia, modified using Matt's code ###
### Date created: April 7, 2020 ###
library(vegan)
library(tidyverse)
library(dplyr)

#############################
############ 16S ############
#############################

### Read table metadata and abundances
microbes_16S <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level.csv", header=T)
names(microbes_16S)[1:17]

# split by year
comm.years <- split( microbes_16S, microbes_16S$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:17)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:17)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:16])
d <- lapply(comm.zero, function(z) z[,-c(1:16)])

x <- lapply( meta, function(z) data.frame(z[,c("quadrat", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
                  d, x )
multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2015:2018),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="R_Code_and_Analysis/output_data/adipart_16S.csv", row.names = F)
g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2015:2018) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/figs/microbes_16S_gamma_partition_additive_richness_prop.png", width=4, height=3 )

#############################
############ 18S ############
#############################

### Read table metadata and abundances
microbes_18S <- read.csv("Data/data_parfrey/18S/18S_ASV_MASTER_Hakai_final.csv", header=T)

names(microbes_18S)[1:13]

# split by year
comm.years <- split( microbes_18S, microbes_18S$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:11)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:11)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:11])
d <- lapply(comm.zero, function(z) z[,-c(1:11)])

x <- lapply( meta, function(z) data.frame(z[,c("quadrat", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

multis <- mapply( function(comm,x) multipart(comm,x, index="tsallis",  relative=T, nsimul = 1000 ), 
                  d, x )
multis.mean <- lapply( multis[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
multis.obs <- lapply(multis[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( multis[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )


add.plot <- data.frame( year=rep(c(2015:2018),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/adipart_18S.csv", row.names=F)

g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2015:2018) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/figs/microbes_18S_gamma_partition_additive_richness_prop.png", width=4, height=3 )

# sampling design (alpha, beta1= within sites, b2=between sites, b3=between regions)
