###
###




### libraries
library(tidyverse)


## read data
invert <- read_csv("output_data/adipart_invert.csv")
invert$type <- "invert"
m16 <- read_csv("output_data/adipart_16S.csv")
m16$type <- "16S"
m18 <- read_csv("output_data/adipart_18S.csv")
m18$type <- "18S"


## merge
d <- full_join(full_join(m18,m16),invert)
d$level <- factor(d$level, levels=c("beta_region","beta_site","beta_sample","alpha"))

## plot
g1 <- ggplot(d %>% filter(key == "observed" & type=="16S"),
             aes(x=year - 0.2, y=value, fill=level)) +
  geom_bar(stat="identity", width=0.1) + 
  scale_x_continuous(breaks=c(1, 2, 3, 4, 5), labels=2014:2018) +
  labs(x="Year", y="Proportion of gamma diversity")
g1 + geom_bar(data=d %>% filter(key == "observed" & type=="18S"),
              aes(x=year + 0, fill=level),
              stat="identity", position="stack", width=0.1, alpha=1) + 
  geom_bar(data=d %>% filter(key == "observed" & type=="invert"),
              aes(x=year + 0.2, fill=level),
              stat="identity", position="stack", width=0.1, alpha=1)
