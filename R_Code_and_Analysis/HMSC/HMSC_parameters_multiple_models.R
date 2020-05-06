#############
### HMSC parameter estimates across models
### author: Matt Whalen
### started: 25 February 2020
###


### start with genus-level microbes only


# load libraries
library( Hmsc )
library( tidyverse )

###
# read models
# # 16S first
# load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2015_genus_chains_4_thin_100_samples_1000.Rdata")
# load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2016_genus_chains_4_thin_100_samples_1000.Rdata")
# load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2017_genus_chains_4_thin_100_samples_1000.Rdata")
# load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2018_genus_chains_4_thin_100_samples_1000.Rdata")
# 18S
load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2015_genus_chains_4_thin_100_samples_1000.Rdata")
load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2016_genus_chains_4_thin_100_samples_1000.Rdata")
load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2017_genus_chains_4_thin_100_samples_1000.Rdata")
load( "output_data/hmsc_models_microbes/16S_models/model_microbes_16S_2018_genus_chains_4_thin_100_samples_1000.Rdata")
 
# wrap these into a list
mlist <- list(mod_2015,mod_2016,mod_2017,mod_2018)


# get parameter estimates
postBeta = function(z){
  getPostEstimate(z, parName = "Beta")
}
betas <- lapply( mlist, postBeta )

# calculate mean and variance of parameter esimates
pbdf <- function(z){
  pbdf = data.frame( t(z$mean), taxon=colnames(z$mean) )
  names(pbdf) <- c("intercept","temp","sal","depth","area",
                 "biomass","LAI","microepi","macro", "taxon")
  return(pbdf)
}

pbdfs <- lapply( betas, pbdf )


coef_plot <- function(pbdf){
  coefs = pbdf %>%
  select( -taxon ) %>%
  gather( coefficient, value )
}

coefs <- lapply(pbdfs, coef_plot )

# combine
coeftab <- data.table::rbindlist(coefs, idcol = "index")

# add year
coeftab$year <- coeftab$index + 2014



### figures
windows(6,4)
ggplot( coeftab, aes(y=value, x=year)) +
  facet_wrap(~coefficient, scales="free_y") +
  geom_hline(yintercept=0) +
  stat_summary( fun.data=mean_cl_boot, geom='errorbar',
                size=1, width=0.1, position=position_dodge(width=0.9) )  +
  stat_summary( fun.y=mean, geom='point',
                size=3, position=position_dodge(width=0.9) )  +
  ylab('coefficient') + xlab('')


