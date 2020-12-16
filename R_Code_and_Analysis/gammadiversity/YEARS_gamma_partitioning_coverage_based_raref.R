### ADIPART ###
### Author: Bianca Trevizan Segovia, modified using Matt's code ###
### Date created: December 3rd, 2020 ###

library(vegan)
library(tidyverse)
library(dplyr)

#############################
############ 16S ############
#############################

### genus level ###

### Read table metadata and abundances
microbes_16S_genus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv", header=T)
names(microbes_16S_genus)[1:16]

# # split by year
# comm.years <- split( microbes_16S_genus, microbes_16S_genus$year )
# lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:12)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:12)],df))
})

meta <- lapply(comm.zero, function(z) z[,1:12])
d <- lapply(comm.zero, function(z) z[,-c(1:12)])

x <- lapply( meta, function(z) data.frame(z[,c("site_quadrat_id", "site", "region")],whole=1) )

adds <- mapply( function(comm,x) adipart(comm,x, index="richness",weights="prop", 
                                         relative=T, nsimul = 20 ), 
                d, x )
adds.mean <- lapply( adds[c(2,4,6,8)], function(z) z$means[-c(2,3,4)] )
adds.obs <- lapply(adds[c(1,3,5,7)], function(z) z[-c(2,3,4)] )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )

add.plot <- data.frame( year=rep(c(2015:2018),each=4),
                        level=rep(c("alpha","beta_sample","beta_site","beta_region"),4), 
                        observed=do.call( c, adds.obs ), null=do.call( c, adds.mean ) )
add.plot$level <- factor( add.plot$level, levels=c("alpha","beta_sample","beta_site","beta_region") )
ggplot( data=add.plot, aes(x=year,fill=level,y=observed)) + geom_bar(stat='identity')

add.plot2 <- add.plot %>% 
  gather( value="value", key = "key", - year, -level)
write.csv(add.plot2, file="R_Code_and_Analysis/output_data/adipart_16S_genus_level.csv", row.names = F)
g1 <- ggplot(add.plot2 %>% dplyr::filter(key == "observed"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4), labels=2015:2018) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=add.plot2 %>% dplyr::filter(key == "null"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.3, alpha=0.5)
ggsave( "R_Code_and_Analysis/gammadiversity/microbes_16S_genus_gamma_partition_additive_richness_prop.png", width=4, height=3 )
