rm(list=ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

library(ggplot2) 
library(rstan)
library(shinystan)

if(length(grep("deirdre", getwd())>0)) {  setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else if
(length(grep("temp", getwd())>0)) {   setwd("~/Documents/git/temp")
}

Nquad <- 15# quadrat
Nsite <- 25 # sites (choked, Triquet etc)
Nyear <- 15# number of year with microbial

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

### Adding climate variables #####

temperature <- rnorm(Nsite*Nyear, 5, 2)
trt.dat$temp <- rep(temperature, each = Nquad) 

mu.force = -10 
sigma.force = 1
alpha.force.sp <- rnorm(Nsite*Nyear, mu.force, sigma.force)
trt.dat$alphaForceSp <- rep(alpha.force.sp, each = Nquad)

salinity <- rnorm(Nsite*Nyear, 10, 2)
trt.dat$salin <- rep(salinity, each = Nquad) 

mu.salin = -3 
sigma.salin = 5
alpha.salin.sp <- rnorm(Nsite*Nyear, mu.salin, sigma.salin)
trt.dat$alphasalinSp <- rep(alpha.salin.sp, each = Nquad)

# muTemp.site <- 12
# sigmaTemp.site <- 3
# alphaTemp.site  <- rnorm(Nsite, muTemp.site, sigmaTemp.site)
# trt.dat$alphaTemp.site <- rep(alphaTemp.site, each = Nquad) 
# 
# betaTraitxforce <- 0 #interaction between trait and phenology
# 
# # beta.temp <- alphaTemp.site + alpha.site * betaTraitxforce
# # beta.temp.site <- rep(beta.temp,)
# # pheno.dat$beta.force.sp <- rep(beta.force.sp, each = nphen)
# 
# #Generate the cue values (ie the F that gets multipled with betaForcing{sp})
# mu.Temp <- 1 # This is a big forcing effect.  turned down to 5
# sigma.Temp <- 1
# temp.i <- rnorm(Nsite, mu.Temp, sigma.Temp)  # predictor frocing, forcei in stan
# trt.dat$temp.i <- rep(temp.i, each = Nquad*Nyear)
# 


# general variance
trt.var <- 1 
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

#trt.dat$yMicroi <- trt.dat$mu.site +  trt.dat$trt.er + mu.year * trt.dat$year + trt.dat$alphaForceSp * trt.dat$temp + trt.dat$alphasalinSp * trt.dat$salin
trt.dat$yMicroi <- trt.dat$mu.site +  trt.dat$trt.er + mu.year * trt.dat$year + mu.force * trt.dat$temp + mu.salin * trt.dat$salin

micro_data <- list(yMicroi = trt.dat$yMicroi,
                   N = Ntrt,
                   n_site = Nsite,
                   micro_site = trt.dat$site,
                   yeari = trt.dat$year,
                   tempi = trt.dat$temp,
                   salini = trt.dat$salin)


mdl.micro <- stan('m2m/stan/macroMdl2.stan',
                  data = micro_data,
                  iter = 4000,
                  warmup = 3000,
                  chains = 4,
                  include = FALSE, pars = c("ypred_new","y_hat")
)

sumer <- summary(mdl.micro)$summary[c("mu_grand","b_year","mu_force","mu_salin","sigma_site","sigma_Microy"),]; sumer

# Nrep <- 25# rep per trait
# Npop <- 8
# Nspp <- 25# number of species with traits (making this 20 just for speed for now)
# 
# # First making a data frame for the test trait data
# Nmacro <- Nspp * Npop * Nrep# total number of traits observations
# Nmacro
# 
# #make a dataframe for height
# macro.data <- data.frame(matrix(NA, Nmacro, 1))
# names(macro.data) <- c("rep")
# macro.data$rep <- c(1:Nrep)
# macro.data$species <- rep(1:Nspp, each = Nrep)
# macro.data$pop <- rep(1:Npop, each = Nspp*Nrep)
# 
# 
# lati <- seq(1:8)
# macro.data$lat <- rep(lati, each = Nrep*Nspp)
# macro.data$lat <- as.numeric(macro.data$lat)
# 
# temperature <- rnorm(Nspp, 5, 2)
# macro.data$temp <- rep(temperature, each = Nrep*Npop)
# 
# muTemp <- 12
# sigmaTemp <- 3
# alphaTemp <- rnorm(Nspp, muTemp, sigmaTemp)
# macro.data$alphaTemp <- rep(alphaTemp, each = Nrep)
# 
# mu.year <- 4
# 
# mu.site <- 10
# sigma.species <- 0.5 # we want to keep the variaiton across spp. high
# mu.site <- rnorm(Nspp, mu.site, sigma.species)
# macro.data$mu.site <- rep(mu.site, each = Nrep) #adding ht data for ea. sp
# 
# # general variance
# macro.var <- 1 #sigma_traity in the stan code
# macro.data$macro.er <- rnorm(Nmacro, 0, macro.var)
# 
# # generate yhat - heights -  for this first macro model
# for (i in 1:Nmacro){
#   macro.data$yTraiti[i] <- 
#     macro.data$mu.site[i] + macro.data$macro.er[i] + mu.year * macro.data$lat[i] + muTemp * macro.data$temp[i]
# }
# 
# all.data <- list(yTraiti = macro.data$yTraiti,
#                  N = Nmacro,
#                  n_spec = Nspp,
#                  trait_species = as.numeric(as.factor(macro.data$species)),
#                  lati = macro.data$lat,
#                  tempi = macro.data$temp)
# 
# all.data <- list(yMicroi = macro.data$yTraiti,
#                    N = Nmacro,
#                    n_site = Nspp,
#                    micro_site = as.numeric(as.factor(macro.data$species)),
#                    yeari = macro.data$lat,
#                  tempi = macro.data$temp)
# 
# mdl <- stan("m2m/stan/macroMdl2.stan",
#                         data = all.data,
#                     iter = 4000, warmup = 3000, chains=4,
#                      include = FALSE, pars = c("y_hat"))
# 
# sumer <- summary(mdl)$summary[c("mu_grand","b_year", "mu_force","sigma_site","sigma_Microy", "sigmaTempS"),]; sumer
# 
# 
# sumer <- summary(mdl)$summary[c("muSp","b_year", "mu_force","sigma_sp", "sigma_traity", "sigmaTempS"),]; sumer
############ macro part #######################################