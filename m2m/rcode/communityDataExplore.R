# started oct 13, 2024 by D. Loughnan

# starting to explore the community data:
# 1. how many sites and what sites do we have community data from
# 2. Calculations of the alpha and beta and omega diversity 

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) 

#setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")
if(length(grep("deirdre", getwd())>0)) {  setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {   setwd("~/Documents/git/temp")
}