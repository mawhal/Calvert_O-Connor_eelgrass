# code started Feb 27, 2025 by D. Loughnan
# aim of this code is to develop test data to jointly model changes in microbial diversity and how it impacts changes in macro diversity.
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

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

Nquad <- 20# quadrat
Nsite <- 30 # sites (choked, Triquet etc)
Nyear <- 30# number of year with microbial

# First making a data frame for the test micro data
Nmicro <- Nyear * Nsite * Nquad# total number of micro observations
Nmicro

#make a data frame 
micro.dat <- data.frame(matrix(NA, Nmicro, 1))
names(micro.dat) <- c("iter"); micro.dat$iter <- 1:nrow(micro.dat)
micro.dat$rep <- c(1:Nquad)
micro.dat$site <- rep(1:Nsite, each = Nquad)
micro.dat$year <- rep(1:Nyear, each = Nquad*Nsite)

mu.grand <- 10 # the grand mean of the micro model
sigma.site <- 1 # we want to keep the variation across sites high

mu.site <- rnorm(Nsite, 0, sigma.site)
micro.dat$mu.site <- rep(mu.site + mu.grand, each = Nquad) 

#now generating the effects of year
yeari <- rnorm(Nsite*Nyear, 10, 1)
micro.dat$yeari <- (rep(yeari, each = Nquad))

mu.year <- 4
sigma.year = 1
alpha.year.site <- rnorm(Nsite, mu.year, sigma.year)
micro.dat$alphayearSite <- rep(alpha.year.site, each = Nquad)

# general variance
micro.var <- 1 
micro.dat$micro.er <- rnorm(Nmicro, 0, micro.var)

#micro.dat$yMicroi <- micro.dat$mu.site +  micro.dat$micro.er + micro.dat$alphayearSite * micro.dat$year

for (i in 1:Nmicro){
  micro.dat$yMicroi2[i] <-  micro.dat$mu.site[i] +  micro.dat$micro.er[i] + micro.dat$alphayearSite[i]  * micro.dat$yeari[i] 
    }

micro_data <- list(yMicroi = micro.dat$yMicroi,
                   N = Nmicro,
                   n_sites = Nsite,
                   micro_sites = micro.dat$site,
                   yeari = micro.dat$yeari
                   )

mdl.micro <- stan('m2m/stan/microMdl.stan',
                  data = micro_data,
                  iter = 4000,
                  warmup = 3000,
                  chains = 4,
                  include = FALSE, pars = c("ypred_new","y_hat")
)

sumer <- summary(mdl.micro)$summary[c("mu_grand","muYear", "sigmaYear","sigma_site","sigma_Microy"),]; sumer

#                  mean      se_mean          sd      2.5%       25%        50%        75%      97.5%     n_eff      Rhat
# mu_grand  10   10.0862238 0.0065285949 0.206007884 9.6774653 9.9533429 10.0839305 10.2190603 10.4924307  995.6994 1.0024579
# muYear    4    4.0026940 0.0020930133 0.161715515 3.6746745 3.8971607  4.0040823  4.1104909  4.3154459 5969.7948 0.9992787
# sigmaYear  1   0.9078530 0.0015954714 0.126805947 0.7021279 0.8166935  0.8944491  0.9822854  1.1887270 6316.8591 0.9997024
# sigma_site  1  0.9999804 0.0031526018 0.166878513 0.7147537 0.8855282  0.9828262  1.1027991  1.3586601 2801.9643 1.0004308
# sigma_Microy 1 0.9891385 0.0000541771 0.005236842 0.9790957 0.9855876  0.9890563  0.9926463  0.9997088 9343.4533 0.9993346

ssm <-  as.shinystan(mdl.micro)
launch_shinystan(ssm)

