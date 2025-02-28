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
Nsite <- 25 # sites (choked, Triquet etc)
Nyear <- 25# number of year with microbial

# First making a data frame for the test micro data
Nmicro <- Nyear * Nsite * Nquad# total number of micro observations
Nmicro

#make a data frame 
m2m.dat <- data.frame(matrix(NA, Nmicro, 1))
names(m2m.dat) <- c("iter"); m2m.dat$iter <- 1:nrow(m2m.dat)
m2m.dat$rep <- c(1:Nquad)
m2m.dat$site <- rep(1:Nsite, each = Nquad)
m2m.dat$year <- rep(1:Nyear, each = Nquad*Nsite)

mu.grand <- 10 # the grand mean of the micro model
sigma.site <- 1 # we want to keep the variation across sites high

mu.site <- rnorm(Nsite, 0, sigma.site)
m2m.dat$mu.site <- rep(mu.site + mu.grand, each = Nquad) 

#now generating the effects of year
yeari <- rnorm(Nsite*Nyear, 10, 1)
m2m.dat$yeari <- (rep(yeari, each = Nquad))

mu.year <- 4
sigma.year = 1
alpha.year.site <- rnorm(Nsite, mu.year, sigma.year)
m2m.dat$alphayearSite <- rep(alpha.year.site, each = Nquad)

# general variance
micro.var <- 1 
m2m.dat$micro.er <- rnorm(Nmicro, 0, micro.var)

#m2m.dat$yMicroi <- m2m.dat$mu.site +  m2m.dat$micro.er + m2m.dat$alphayearSite * m2m.dat$year

for (i in 1:Nmicro){
  m2m.dat$yMicroi[i] <-  m2m.dat$mu.site[i] +  m2m.dat$micro.er[i] + m2m.dat$alphayearSite[i]  * m2m.dat$yeari[i] 
}

#######################################################
Nmacro <- Nquad*Nyear*Nsite 
# 
# m2m.dat <- data.frame(matrix(NA, Nmacro, 2))
# names(m2m.dat) <- c("rep","sites")
# m2m.dat$quad <- c(1:Nquad)
# m2m.dat$sites <- rep(c(1:Nsite), each = Nquad)
# m2m.dat$year <- rep(1:Nyear, each = Nsite*Nquad)

tempi <- rnorm(Nsite*Nyear, 5, 2)
m2m.dat$tempi <- rep(tempi, each = Nquad)
#m2m.dat$tempi <- rnorm(Nmacro, 5,2)

salini <- rnorm(Nsite*Nyear, 8, 3)
m2m.dat$salini <- rep(salini, each = Nquad)
#m2m.dat$salini <- rnorm(Nmacro, 8,2)

biomassi <- rnorm(Nsite*Nyear, 6, 4)
m2m.dat$biomassi <- rep(rep(biomassi, each = Nquad))
#m2m.dat$biomassi <- rnorm(Nmacro, 6,4)

macroAi <- rnorm(Nsite*Nyear, 15, 5)
m2m.dat$macroAi <- rep(rep(macroAi, each = Nquad))
#m2m.dat$macroAi <- rnorm(Nmacro, 15,5)

bedAreai <- rnorm(Nsite*Nyear, 10, 1)
m2m.dat$bedAreai <- rep(rep(bedAreai, each = Nquad))
#m2m.dat$bedAreai <- rnorm(Nmacro, 10,1)

depthi <- rnorm(Nsite*Nyear, 10, 1)
m2m.dat$depthi <- rep(rep(depthi, each = Nquad))
#m2m.dat$depthi <- rnorm(Nmacro, 5,2)

distEi <- rnorm(Nsite*Nyear, 10, 1)
m2m.dat$distEi <- rep(rep(distEi, each = Nquad))
#m2m.dat$distEi <- rnorm(Nmacro, 3,4)

mu.year = -4
sigma.year = 1
alpha.year.site <- rnorm(Nsite, mu.year, sigma.year)
m2m.dat$alphayearSite <- rep(alpha.year.site, each = Nquad)

mu.temp = -1
sigma.temp = 1
alpha.temp.site <- rnorm(Nsite, mu.temp, sigma.temp)
m2m.dat$alphatempSite <- rep(alpha.temp.site, each = Nquad)

mu.salin = 2
sigma.salin = 2
alpha.salin.site <- rnorm(Nsite, mu.salin, sigma.salin)
m2m.dat$alphasalinSite <- rep(alpha.salin.site, each = Nquad)

mu.biomass = 4
sigma.biomass = 3
alpha.biomass.site <- rnorm(Nsite, mu.biomass, sigma.biomass)
m2m.dat$alphabiomassSite <- rep(alpha.biomass.site, each = Nquad)

mu.macroA = 6
sigma.macroA = 4
alpha.macroA.site <- rnorm(Nsite, mu.macroA, sigma.macroA)
m2m.dat$alphamacroASite <- rep(alpha.macroA.site, each = Nquad)

mu.bedArea = 8
sigma.bedArea = 5
alpha.bedArea.site <- rnorm(Nsite, mu.bedArea, sigma.bedArea)
m2m.dat$alphabedAreaSite <- rep(alpha.bedArea.site, each = Nquad)

mu.depth = 10
sigma.depth = 6
alpha.depth.site <- rnorm(Nsite, mu.depth, sigma.depth)
m2m.dat$alphadepthSite <- rep(alpha.depth.site, each = Nquad)

mu.distE = 12
sigma.distE = 7
alpha.distE.site <- rnorm(Nsite, mu.distE, sigma.distE)
m2m.dat$alphadistESite <- rep(alpha.distE.site, each = Nquad)

########
mu.macro.site = 80
sigma.macro.site = 10
alphamacroSite <- rnorm(Nsite, mu.macro.site, sigma.macro.site)
m2m.dat$alphamacroSite <- rep(alphamacroSite, each = Nquad)

sigma.macro.y = 3
m2m.dat$eMacro <- rnorm(Nmacro, 0, sigma.macro.y)

##########
betaMicroxV1 <- 3
betaMicroxV2 <- 2
betaMicroxV3 <- 1

for (i in 1:Nmacro){
  m2m.dat$betaV1[i] <-  m2m.dat$alphayearSite[i] + (betaMicroxV1 *  m2m.dat$mu.site[i])
  
  m2m.dat$betaV2[i]<- m2m.dat$alphatempSite[i] + (betaMicroxV2*  m2m.dat$mu.site[i])
  
  m2m.dat$betaV3[i] <-m2m.dat$alphasalinSite[i] + (betaMicroxV3* m2m.dat$mu.site[i])
}


for (i in 1:Nmacro){
  m2m.dat$yMacroi[i] <-  m2m.dat$alphamacroSite[i] + m2m.dat$eMacro[i] + 
    #m2m.dat$betaV1[i] * m2m.dat$yeari[i] + m2m.dat$betaV2[i] * m2m.dat$tempi[i]  +  m2m.dat$betaV3[i] * m2m.dat$salini[i]   
    + m2m.dat$alphabiomassSite[i] * m2m.dat$biomassi[i] + m2m.dat$alphamacroASite[i] * m2m.dat$macroAi[i] +   m2m.dat$alphabedAreaSite[i] * m2m.dat$bedAreai[i] + m2m.dat$alphadepthSite[i] * m2m.dat$depthi[i] +  m2m.dat$alphadistESite[i] * m2m.dat$distEi[i] 
}

all.data <- list(
  yMicroi = m2m.dat$yMicroi,
  N = nrow(m2m.dat),
  n_sites = Nsite,
  micro_sites = m2m.dat$site,
  yeari = m2m.dat$yeari,
  ########################
  Nmacro = nrow(m2m.dat),
  macro_sites = as.numeric(as.factor(m2m.dat$site)),
  sites = as.numeric(as.factor(m2m.dat$site)),
  yMacroi = m2m.dat$yMacroi,
  Var1i = m2m.dat$yeari,
  Var2i = m2m.dat$tempi,
  Var3i = m2m.dat$salini,
  Var4i = m2m.dat$biomassi,
  Var5i = m2m.dat$macroAi,
  Var6i = m2m.dat$bedAreai,
  Var7i = m2m.dat$depthi,
  Var8i = m2m.dat$distEi)

mdl <- stan('m2m/stan/sandBoxStan.stan',
                  data =all.data,
                  iter = 4000,
                  warmup = 3000,
                  chains = 4,
                  include = FALSE, pars = c("ypred_new","y_hat")
)

sumer <- summary(mdl)$summary[c("mu_grand","muYear", "sigmaYear","sigma_site","sigma_Microy",
                                # "muVar1",
                                # "muVar2",
                                # "muVar3",
                                "muVar4",
                                "muVar5",
                                "muVar6",
                                "muVar7",
                                "muVar8",
                                "muMacroSites",
                                # "sigmaVar1", 
                                # "sigmaVar2", 
                                # "sigmaVar3",  
                                "sigmaVar4", 
                                "sigmaVar5",
                                "sigmaVar6",
                                "sigmaVar7",
                                "sigmaVar8",
                                "sigmaMacroSites",
                                "sigmaMacro_y"),]; sumer

