###
###




### libraries
library(tidyverse)


## read data
invert <- read_csv("R_Code_and_Analysis/output_data/adipart_invert.csv")
invert$type <- "invert"
m16 <- read_csv("R_Code_and_Analysis/output_data/adipart_16S.csv")
m16$type <- "16S"
m18 <- read_csv("R_Code_and_Analysis/output_data/adipart_18S.csv")
m18$type <- "18S"
  
## merge
d <- full_join(full_join(m18,m16),invert)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="16S"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity") +
  scale_fill_viridis_d()
  # scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="18S"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="18S"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="invert"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="invert"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "figs/adipart_all.svg", width=6, height=4)
ggsave( "figs/adipart_all.png", width=6, height=4)


#test with family level 16S
m16_family <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_family_level.csv")
m16_family$type <- "16S_family"

## merge
d <- full_join(full_join(m18,m16_family),invert)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="16S_family"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity") +
  scale_fill_viridis_d()
# scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="18S"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="18S"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="invert"),
           aes(x=year + 0.15, fill=level),
           stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="invert"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "R_Code_and_Analysis/figs/adipart_all_16S_family_level.png", width=6, height=4)

#test with genus level 16S
m16_genus <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_genus_level.csv")
m16_genus$type <- "16S_genus"

## merge
d <- full_join(full_join(m18,m16_genus),invert)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="16S_genus"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity") +
  scale_fill_viridis_d()
# scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="18S"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="18S"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="invert"),
           aes(x=year + 0.15, fill=level),
           stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="invert"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "R_Code_and_Analysis/figs/adipart_all_16S_genus_level.png", width=6, height=4)

