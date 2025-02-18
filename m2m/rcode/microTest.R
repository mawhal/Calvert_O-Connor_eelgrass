rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

library(ggplot2) 
library(rstan)
library(shinystan)

#setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")
if(length(grep("deirdre", getwd())>0)) {  setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {   setwd("~/Documents/git/temp")
}

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

Nquad <- 20# quadrat
Nsite <- 25 # sites (choked, Triquet etc)
Nyear <- 25# number of year with microbial

# First making a data frame for the test micro data
Ntrt <- Nyear * Nsite * Nquad# total number of micro observations
Ntrt

#make a data frame 
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("iter"); trt.dat$iter <- 1:nrow(trt.dat)
trt.dat$rep <- c(1:Nquad)
trt.dat$year <- rep(1:Nyear, each = Nquad)
trt.dat$site <- rep(1:Nsite, each = Nyear*Nquad)

mu.grand <- 10 # the grand mean of the micro model
sigma.site <- 1 # we want to keep the variation across sites high

mu.site <- rnorm(Nsite, 0, sigma.site)
trt.dat$mu.site <- rep(mu.site + mu.grand, each = Nquad*Nyear) 

#now generating the effects of year
mu.year <- 4

# general variance
trt.var <- 1 
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

trt.dat$yMicroi <- trt.dat$mu.site +  trt.dat$trt.er + mu.year * trt.dat$year

micro_data <- list(yMicroi = trt.dat$yMicroi,
                   N = Ntrt,
                   n_site = Nsite,
                   micro_site = trt.dat$site,
                   yeari = trt.dat$year
                   )

mdl.micro <- stan('m2m/stan/microMdl.stan',
                  data = micro_data,
                  iter = 4000,
                  warmup = 3000,
                  chains = 4,
                  include = FALSE, pars = c("ypred_new","y_hat")
)

sumer <- summary(mdl.micro)$summary[c("mu_grand","b_year","sigma_site","sigma_Microy"),]; sumer


ssm <-  as.shinystan(mdl.micro)
launch_shinystan(ssm)

