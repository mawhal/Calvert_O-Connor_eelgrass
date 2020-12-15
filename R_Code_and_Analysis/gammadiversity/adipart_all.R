###
###

### comparing microbes at the genus level and inverts at the finest


### libraries
library(tidyverse)
library(ggplot2)

## read data
invert_finest <- read_csv("R_Code_and_Analysis/output_data/adipart_inverts_finest_level.csv")
invert_finest$type <- "macro"
m16_genus <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_genus_level.csv")
m16_genus$type <- "prok"
m18_genus <- read_csv("R_Code_and_Analysis/output_data/adipart_18S_genus_level.csv")
m18_genus$type <- "micro"
# macroeuk18_genus <- read_csv("R_Code_and_Analysis/output_data/adipart_macroeuk18S_ASV_level.csv")
# macroeuk18_genus$type <- "macroeuk18S"
  
## merge
d2 <- full_join(full_join(m18_genus,m16_genus),invert_finest)
# d2 <- full_join(d,macroeuk18_genus )
d2$level <- factor(d2$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d2 %>% filter(key == "observed" & type=="prok"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity", title = "Microbes genus Inverts finest taxonomic level") +
  scale_fill_viridis_d()
  # scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 <- g1 + geom_bar(data=d2 %>% filter(key == "observed" & type=="micro"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d2 %>% filter(key == "observed" & type=="micro"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d2 %>% filter(key == "observed" & type=="macro"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d2 %>% filter(key == "observed" & type=="macro"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 

# g2 <- g1 +geom_bar(data=d2 %>% filter(key == "observed" & type=="macroeuk18S"),
#                aes(x=year + 0.30, fill=level),
#                stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
#   geom_text( data=d2 %>% filter(key == "observed" & type=="macroeuk18S"), aes(x=year+0.30,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 )
# g2
  
ggsave( "R_Code_and_Analysis/gammadiversity/adipart_all_microbes_genus_inverts_finest.png", plot=g1,width=6, height=4)


#test with family level
m16_family <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_family_level.csv")
m16_family$type <- "prok"

m18_family <- read_csv("R_Code_and_Analysis/output_data/adipart_18S_family_level.csv")
m18_family$type <- "micro"

invert_family <- read_csv("R_Code_and_Analysis/output_data/adipart_inverts_family_level.csv")
invert_family$type <- "macro"

## merge
d <- full_join(full_join(m18_family,m16_family),invert_family)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="prok"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity", title = "Family level") +
  scale_fill_viridis_d()
# scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="micro"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="micro"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="macro"),
           aes(x=year + 0.15, fill=level),
           stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="macro"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "R_Code_and_Analysis/gammadiversity/adipart_all_family_level.png", width=6, height=4)
