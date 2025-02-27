# Started December 15, 2024: DL
#the aim of this test data is to double check that we are confident in the SSD values

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
require(rstan)
#require(truncnorm)

Nquad <- 5# N quad
Nyear <- 10# N year
Nsite <- 50# N sites

Nmacro <- Nquad*Nyear*Nsite 

macro.dat <- data.frame(matrix(NA, Nmacro, 2))
names(macro.dat) <- c("rep","sites")
macro.dat$quad <- c(1:Nquad)
macro.dat$sites <- rep(c(1:Nsite), each = Nquad)
macro.dat$year <- rep(1:Nyear, each = Nsite*Nquad)

yeari <- rnorm(Nsite*Nyear, 10, 1)
macro.dat$yeari <- (rep(yeari, each = Nquad))
#macro.dat$yeari <- rnorm(Nmacro, 10, 1)

tempi <- rnorm(Nsite*Nyear, 5, 2)
macro.dat$tempi <- rep(tempi, each = Nquad)
#macro.dat$tempi <- rnorm(Nmacro, 5,2)

salini <- rnorm(Nsite*Nyear, 8, 3)
macro.dat$salini <- rep(salini, each = Nquad)
#macro.dat$salini <- rnorm(Nmacro, 8,2)

biomassi <- rnorm(Nsite*Nyear, 6, 4)
macro.dat$biomassi <- rep(rep(biomassi, each = Nquad))
#macro.dat$biomassi <- rnorm(Nmacro, 6,4)

macroAi <- rnorm(Nsite*Nyear, 15, 5)
macro.dat$macroAi <- rep(rep(macroAi, each = Nquad))
#macro.dat$macroAi <- rnorm(Nmacro, 15,5)

bedAreai <- rnorm(Nsite*Nyear, 10, 1)
macro.dat$bedAreai <- rep(rep(bedAreai, each = Nquad))
#macro.dat$bedAreai <- rnorm(Nmacro, 10,1)

depthi <- rnorm(Nsite*Nyear, 10, 1)
macro.dat$depthi <- rep(rep(depthi, each = Nquad))
#macro.dat$depthi <- rnorm(Nmacro, 5,2)

distEi <- rnorm(Nsite*Nyear, 10, 1)
macro.dat$distEi <- rep(rep(distEi, each = Nquad))
#macro.dat$distEi <- rnorm(Nmacro, 3,4)

mu.year = -4
sigma.year = 1
alpha.year.site <- rnorm(Nsite, mu.year, sigma.year)
macro.dat$alphayearSite <- rep(alpha.year.site, each = Nquad)

mu.temp = -1
sigma.temp = 1
alpha.temp.site <- rnorm(Nsite, mu.temp, sigma.temp)
macro.dat$alphatempSite <- rep(alpha.temp.site, each = Nquad)

mu.salin = 2
sigma.salin = 2
alpha.salin.site <- rnorm(Nsite, mu.salin, sigma.salin)
macro.dat$alphasalinSite <- rep(alpha.salin.site, each = Nquad)

mu.biomass = 4
sigma.biomass = 3
alpha.biomass.site <- rnorm(Nsite, mu.biomass, sigma.biomass)
macro.dat$alphabiomassSite <- rep(alpha.biomass.site, each = Nquad)

mu.macroA = 6
sigma.macroA = 4
alpha.macroA.site <- rnorm(Nsite, mu.macroA, sigma.macroA)
macro.dat$alphamacroASite <- rep(alpha.macroA.site, each = Nquad)

mu.bedArea = 8
sigma.bedArea = 5
alpha.bedArea.site <- rnorm(Nsite, mu.bedArea, sigma.bedArea)
macro.dat$alphabedAreaSite <- rep(alpha.bedArea.site, each = Nquad)

mu.depth = 10
sigma.depth = 6
alpha.depth.site <- rnorm(Nsite, mu.depth, sigma.depth)
macro.dat$alphadepthSite <- rep(alpha.depth.site, each = Nquad)

mu.distE = 12
sigma.distE = 7
alpha.distE.site <- rnorm(Nsite, mu.distE, sigma.distE)
macro.dat$alphadistESite <- rep(alpha.distE.site, each = Nquad)

########
mu.macro.site = 80
sigma.macro.site = 10
alphamacroSite <- rnorm(Nsite, mu.macro.site, sigma.macro.site)
macro.dat$alphamacroSite <- rep(alphamacroSite, each = Nquad)

sigma.macro.y = 3
macro.dat$eMacro <- rnorm(Nmacro, 0, sigma.macro.y)

for (i in 1:Nmacro){
  macro.dat$yMu[i] <-  macro.dat$alphamacroSite[i] + macro.dat$alphatempSite[i] * macro.dat$tempi[i]  +   macro.dat$alphayearSite[i] * macro.dat$yeari[i]  +  macro.dat$alphasalinSite[i] * macro.dat$salini[i]   + macro.dat$alphabiomassSite[i] * macro.dat$biomassi[i] + macro.dat$alphamacroASite[i] * macro.dat$macroAi[i] +   macro.dat$alphabedAreaSite[i] * macro.dat$bedAreai[i] + macro.dat$alphadepthSite[i] * macro.dat$depthi[i] +  macro.dat$alphadistESite[i] * macro.dat$distEi[i]
}

macro.dat$ymacroi <- macro.dat$yMu + macro.dat$eMacro

all.data <- list(
  n_sites = Nsite,
  Nmacro = nrow(macro.dat),
  macro_sites = as.numeric(as.factor(macro.dat$sites)),
  sites = as.numeric(as.factor(macro.dat$sites)),
  yMacroi = macro.dat$ymacroi,
  Var1i = macro.dat$yeari,
  Var2i = macro.dat$tempi,
  Var3i = macro.dat$salini,
  Var4i = macro.dat$biomassi,
  Var5i = macro.dat$macroAi,
  Var6i = macro.dat$bedAreai,
  Var7i = macro.dat$depthi,
  Var8i = macro.dat$distEi)

mdl <- stan("m2m/stan/macroMdl.stan",
            data = all.data,
            iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)


sumer <- data.frame(summary(mdl)$summary[c("muVar1Sp","muVar2Sp",
                                           "muVar3Sp","muVar4Sp",
                                           "muVar5Sp",
                                           "muVar6Sp",
                                           "muVar7Sp", "muVar8Sp",
                                           "muMacroSites",
                                           "sigmaVar1Sp", "sigmaVar2Sp", 
                                           "sigmaVar3Sp",  
                                           "sigmaVar4Sp", 
                                           "sigmaVar5Sp",
                                           "sigmaVar6Sp",
                                           "sigmaVar7Sp",
                                           "sigmaVar8Sp",
                                           "sigmaMacroSites",
                                           "sigmaMacro_y"),])

mdlparam <- data.frame(var = c("muVar1Sp",
                               "muVar2Sp",
                               "muVar3Sp", 
                               "muVar4Sp",
                               "muVar5Sp",
                               "muVar6Sp",
                               "muVar7Sp",
                               "muVar8Sp",
                               "muMacroSites",
                               "sigmaVar1Sp",
                               "sigmaVar2Sp", 
                               "sigmaVar3Sp",  
                               "sigmaVar4Sp", 
                               "sigmaVar5Sp",
                               "sigmaVar6Sp",
                               "sigmaVar7Sp",
                               "sigmaVar8Sp",
                               "sigmaMacroSites",
                               "sigmaMacro_y"), 
                       param = c(-4, -1, 2, 4,
                                 6, 
                                 8, 10, 12, 
                                 80, 1,1,
                                 2, 3, 4,
                                 5, 6, 7,
                                 10,3))
output <- cbind(mdlparam, sumer); output

# increasing the number of sites makes the esti better
# var param       mean      se_mean         sd      X2.5.       X25.       X50.       X75.     X97.5.    n_eff      Rhat
# muVar1Sp               muVar1Sp    -4 -4.0023767 0.0026930889 0.18352711 -4.3525220 -4.1249375 -4.0028314 -3.8834188 -3.6377232 4644.073 1.0000532
# muVar2Sp               muVar2Sp    -1 -0.9746716 0.0015744385 0.11790327 -1.2086029 -1.0519996 -0.9743727 -0.8989400 -0.7393154 5607.900 0.9997644
# muVar3Sp               muVar3Sp     2  1.9786049 0.0035468970 0.28718328  1.4171466  1.7847336  1.9806937  2.1665156  2.5436466 6555.731 0.9995232
# muVar4Sp               muVar4Sp     4  3.5684797 0.0053478998 0.41340302  2.7341751  3.2931297  3.5704062  3.8467887  4.3746121 5975.590 0.9994095
# muVar5Sp               muVar5Sp     6  5.9891053 0.0063330463 0.50646236  5.0226872  5.6457619  5.9795336  6.3338018  7.0002092 6395.419 0.9995152
# muVar6Sp               muVar6Sp     8  8.4527672 0.0093955559 0.73008927  6.9994691  7.9684444  8.4599630  8.9377680  9.8643685 6038.193 0.9995373
# muVar7Sp               muVar7Sp    10 10.0948914 0.0126654286 1.00445504  8.1113880  9.4376902 10.1054960 10.7468964 12.0879546 6289.574 0.9996285
# muVar8Sp               muVar8Sp    12 12.1628155 0.0140729451 1.11253018 10.0061238 11.4310504 12.1661783 12.9004764 14.3576516 6249.620 0.9996243
# muMacroSites       muMacroSites    80 79.3694290 0.0300546466 1.73862010 75.8155804 78.2513147 79.3520895 80.5249583 82.7871960 3346.464 1.0004537
# sigmaVar1Sp         sigmaVar1Sp     1  1.2165160 0.0023767128 0.14567101  0.9654889  1.1179757  1.2042189  1.3064342  1.5337094 3756.582 1.0006642
# sigmaVar2Sp         sigmaVar2Sp     1  0.7831381 0.0013099170 0.09048269  0.6269771  0.7200242  0.7763018  0.8394386  0.9845740 4771.375 0.9996443
# sigmaVar3Sp         sigmaVar3Sp     2  2.0330000 0.0027808546 0.21596387  1.6531419  1.8782239  2.0180772  2.1699707  2.5003745 6031.226 0.9997777
# sigmaVar4Sp         sigmaVar4Sp     3  2.9307296 0.0043135089 0.31411872  2.4009310  2.7092071  2.8961929  3.1317208  3.6077536 5303.057 0.9993331
# sigmaVar5Sp         sigmaVar5Sp     4  3.5751086 0.0051938612 0.36114772  2.9667593  3.3199337  3.5431074  3.8033633  4.3706331 4834.918 0.9995633
# sigmaVar6Sp         sigmaVar6Sp     5  5.2418839 0.0074112311 0.54352239  4.3221037  4.8522943  5.1968504  5.5848382  6.4380993 5378.413 0.9993760
# sigmaVar7Sp         sigmaVar7Sp     6  7.0911610 0.0101562992 0.74431202  5.7517493  6.5749577  7.0356736  7.5331018  8.6916908 5370.801 0.9992330
# sigmaVar8Sp         sigmaVar8Sp     7  7.8151316 0.0105000783 0.82492498  6.3891709  7.2380248  7.7449791  8.3342669  9.6423316 6172.255 0.9995617
# sigmaMacroSites sigmaMacroSites    10 10.2687278 0.0391630075 1.58804943  7.4189378  9.1506947 10.1817408 11.3074740 13.5644587 1644.281 1.0007174
# sigmaMacro_y       sigmaMacro_y     3  3.0482497 0.0007204463 0.04797135  2.9545088  3.0150757  3.0473229  3.0807222  3.1431339 4433.642 1.0002165# 



