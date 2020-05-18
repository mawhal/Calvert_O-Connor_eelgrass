###
###




### libraries
library(tidyverse)


## read data
invert_finest <- read_csv("R_Code_and_Analysis/output_data/adipart_invert_finest.csv")
invert_finest$type <- "macroeuk_finest"
m16_ASV <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_ASV_level.csv")
m16_ASV$type <- "prokaryotes_ASV"
m18_ASV <- read_csv("R_Code_and_Analysis/output_data/adipart_18S_ASV_level.csv")
m18_ASV$type <- "microeuk_ASV"
  
## merge
d <- full_join(full_join(m18_ASV,m16_ASV),invert_finest)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="prokaryotes_ASV"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity") +
  scale_fill_viridis_d()
  # scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="microeuk_ASV"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="microeuk_ASV"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="macroeuk_finest"),
              aes(x=year + 0.15, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="macroeuk_finest"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "R_Code_and_Analysis/figs/adipart_all_finest.svg", width=6, height=4)
ggsave( "R_Code_and_Analysis/figs/adipart_all_finest.png", width=6, height=4)


#test with family level
m16_family <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_family_level.csv")
m16_family$type <- "prokaryotes_family"

m18_family <- read_csv("R_Code_and_Analysis/output_data/adipart_18S_family_level.csv")
m18_family$type <- "microeuk_family"

invert_family <- read_csv("R_Code_and_Analysis/output_data/adipart_invert_family.csv")
invert_family$type <- "macroeuk_family"

## merge
d <- full_join(full_join(m18_family,m16_family),invert_family)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="prokaryotes_family"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity") +
  scale_fill_viridis_d()
# scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="microeuk_family"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="microeuk_family"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="macroeuk_family"),
           aes(x=year + 0.15, fill=level),
           stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="macroeuk_family"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "R_Code_and_Analysis/figs/adipart_all_family_level.png", width=6, height=4)


#test with genus level
m16_genus <- read_csv("R_Code_and_Analysis/output_data/adipart_16S_genus_level.csv")
m16_genus$type <- "prokaryotes_genus"

m18_genus <- read_csv("R_Code_and_Analysis/output_data/adipart_18S_genus_level.csv")
m18_genus$type <- "microeuk_genus"

invert_finest <- read_csv("R_Code_and_Analysis/output_data/adipart_invert_finest.csv")
invert_finest$type <- "macroeuk_finest"

## merge
d <- full_join(full_join(m18_genus,m16_genus),invert_finest)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="prokaryotes_genus"),
             aes(x=year - 0.15, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.15, col='black',alpha=0.75) + 
  geom_text( aes(y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  scale_x_continuous(limits=c(2013.9,2018.1), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity") +
  scale_fill_viridis_d()
# scale_fill_manual( values=rev(c("black","gray25","gray75","whitesmoke")) )
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="microeuk_genus"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.15, alpha=0.75, col="black") + 
  geom_text( data=d %>% filter(key == "observed" & type=="microeuk_genus"), aes(x=year+0,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) +
  geom_bar(data=d %>% filter(key == "observed" & type=="macroeuk_finest"),
           aes(x=year + 0.15, fill=level),
           stat="identity", position="stack", width=0.15, alpha=0.75, col='black') +
  geom_text( data=d %>% filter(key == "observed" & type=="macroeuk_finest"), aes(x=year+0.15,y=0.005,label=type), angle=90, hjust=0,size=3,fontface=3 ) 
ggsave( "R_Code_and_Analysis/figs/adipart_all_microbes_genus_macro_finest_level.png", width=6, height=4)

