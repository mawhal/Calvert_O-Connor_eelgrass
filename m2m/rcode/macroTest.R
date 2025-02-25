# Started December 15, 2024: DL
#the aim of this test data is to double check that we are confident in the SSD values

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
require(rstan)
#require(truncnorm)

Nrep <- 10# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 80# number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Npop * Nrep# total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp, each = Nrep)
trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep)
trt.dat$tran <- rep(1:Ntran, each = 4*Nrep*Nspp)

lati <- rnorm(8, 50, 5)
trt.dat$lat <- rep(lati, each = Nrep*Nspp)
trt.dat$lat <- as.numeric(trt.dat$lat)
trt.dat$latZ <- (trt.dat$lat-mean(trt.dat$lat,na.rm=TRUE))/(sd(trt.dat$lat,na.rm=TRUE))

mu.tranE <- 5
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

mu.tranlat  <- 2
trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$latZ)

sigma.species <- 10 # we want to keep the variation across spp. high

mu.trtsp <- rnorm(Nspp, 50, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp
#hist(mu.trtsp)
# general variance
trt.var <- 5 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <- 
    trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i] + mu.tranlat * (trt.dat$dumE[i]*trt.dat$latZ[i])  
}

all.data <- list(yTraiti = trt.dat$yTraiti,
                 N = Ntrt,
                 n_spec = Nspp,
                 trait_species = as.numeric(as.factor(trt.dat$species)),
                 lati = trt.dat$latZ,
                 tranE = as.numeric(trt.dat$dumE)
)


# mdl <- stan("stan/modelDevelopment/justDummyIntTrait.stan",
#   data = all.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )

# sumer <- summary(mdl)$summary
# 
# bTranE <- sumer[grep("b_tranE", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# bTranLat <- sumer[grep("b_tranlat", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# 
# sigma_sp <- sumer[grep("sigma_sp", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# sigma_traity <- sumer[grep("sigma_traity", rownames(sumer)), c("mean","2.5%","25%","50%", "75%","97.5%")]
# 
# 
# mdl.out <- data.frame( "Parameter" = c("bTranE","bTranLat","sigma_sp","sigma_traity"),
#   "Test.data.values" = c( mu.tranE, mu.tranlat, sigma.species, trt.var) ,
#   "Estiamte"= c(bTranE[1], bTranLat[1],  sigma_sp[1], sigma_traity[1]),
#   "2.5"= c(bTranE[2], bTranLat[2],  sigma_sp[2], sigma_traity[2]),
#   "25"= c( bTranE[3], bTranLat[3],  sigma_sp[3], sigma_traity[3]),
#   "50"= c( bTranE[4], bTranLat[4],  sigma_sp[4], sigma_traity[4]),
#   "75"= c( bTranE[5], bTranLat[5],  sigma_sp[5], sigma_traity[5]),
#   "97.5"= c(bTranE[6], bTranLat[6],  sigma_sp[6], sigma_traity[6]))
# 
# mdl.out

# 20 reps, 40 spp
# Parameter Test.data.values   Estiamte       X2.5        X25        X50        X75     X97.5
# 1       bTranE              5.0 4.99576221 4.97106368 4.98668986 4.99569333 5.00464410 5.0206824
# 2     bTranLat              2.0 1.99877138 1.98263847 1.99282010 1.99871669 2.00444280 2.0157245
# 3     sigma_sp              0.1 0.09051132 0.06862914 0.08174978 0.08946797 0.09844331 0.1175569
# 4 sigma_traity              0.5 0.50371058 0.49502363 0.50045952 0.50361544 0.50684106 0.5126634

##### Phenology test data ###########################
Nrep <- 5# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 30# number of species with traits (making this 20 just for speed for now)

Nchill <- 8 # sm high low, mp high low
Nphoto <- 2 # high low
Nforce <- 2 # high amd low

# Nrep <- 20
# Nspp <- 15

Nph <- Nrep*Npop*Nspp*Nchill #bc N force and N chill are completely confounded

pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nrep)
pheno.dat$species <- rep(c(1:Nspp), each = Nrep)
pheno.dat$pop <- rep(1:Npop, each = Nspp*Nrep)

# chilli <- rnorm(Nspp*Npop, 10, 1)
# pheno.dat$chilli <- rep(rep(chilli, each = Nrep*Npop))
pheno.dat$chilli <- rnorm(Nph, 10, 1)

# forcei <- rnorm(Nchill, 5, 2)
# pheno.dat$forcei <- rep(rep(forcei, each = Nspp*Nrep))
pheno.dat$forcei <- rnorm(Nph, 5,2)

# salini <- rnorm(Nchill, 8, 3)
# pheno.dat$salini <- rep(rep(salini, each = Nspp*Nrep))
pheno.dat$salini <- rnorm(Nph, 8,2)

# biomassi <- rnorm(Nchill, 6, 4)
# pheno.dat$biomassi <- rep(rep(biomassi, each = Nspp*Nrep))
pheno.dat$biomassi <- rnorm(Nph, 6,4)

# macroAi <- rnorm(Nchill, 15, 5)
# pheno.dat$macroAi <- rep(rep(macroAi, each = Nspp*Nrep))
pheno.dat$macroAi <- rnorm(Nph, 15,5)

# bedAreai <- rnorm(Nchill, 10, 1)
# pheno.dat$bedAreai <- rep(rep(bedAreai, each = Nspp*Nrep))
pheno.dat$bedAreai <- rnorm(Nph, 10,1)

# depthi <- rnorm(Nchill, 10, 1)
# pheno.dat$depthi <- rep(rep(depthi, each = Nspp*Nrep))
pheno.dat$depthi <- rnorm(Nph, 5,2)

# distEi <- rnorm(Nchill, 10, 1)
# pheno.dat$distEi <- rep(rep(distEi, each = Nspp*Nrep))
pheno.dat$distEi <- rnorm(Nph, 3,4)

mu.chill = -4
sigma.chill = 1
alpha.chill.sp <- rnorm(Nspp, mu.chill, sigma.chill)
pheno.dat$alphaChillSp <- rep(alpha.chill.sp, each = Nrep)

mu.force = -1
sigma.force = 1
alpha.force.sp <- rnorm(Nspp, mu.force, sigma.force)
pheno.dat$alphaForceSp <- rep(alpha.force.sp, each = Nrep)

########
mu.salin = 2
sigma.salin = 2
alpha.salin.sp <- rnorm(Nspp, mu.salin, sigma.salin)
pheno.dat$alphasalinSp <- rep(alpha.salin.sp, each = Nrep)

mu.biomass = 4
sigma.biomass = 3
alpha.biomass.sp <- rnorm(Nspp, mu.biomass, sigma.biomass)
pheno.dat$alphabiomassSp <- rep(alpha.biomass.sp, each = Nrep)

# 5
mu.macroA = 6
sigma.macroA = 4
alpha.macroA.sp <- rnorm(Nspp, mu.macroA, sigma.macroA)
pheno.dat$alphamacroASp <- rep(alpha.macroA.sp, each = Nrep)

mu.bedArea = 8
sigma.bedArea = 5
alpha.bedArea.sp <- rnorm(Nspp, mu.bedArea, sigma.bedArea)
pheno.dat$alphabedAreaSp <- rep(alpha.bedArea.sp, each = Nrep)

mu.depth = 10
sigma.depth = 6
alpha.depth.sp <- rnorm(Nspp, mu.depth, sigma.depth)
pheno.dat$alphadepthSp <- rep(alpha.depth.sp, each = Nrep)

mu.distE = 12
sigma.distE = 7
alpha.distE.sp <- rnorm(Nspp, mu.distE, sigma.distE)
pheno.dat$alphadistESp <- rep(alpha.distE.sp, each = Nrep)

########
mu.pheno.sp = 80
sigma.pheno.sp = 10
alphaPhenoSp <- rnorm(Nspp, mu.pheno.sp, sigma.pheno.sp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = Nrep)

sigma.pheno.y = 3
pheno.dat$ePhen <- rnorm(Nph, 0, sigma.pheno.y)


#pheno.datTrait <- merge(pheno.dat, unique(trt.dat[,c("species","mu.trtsp")]), by = "species")
#head(pheno.datTrait,50)

for (i in 1:Nph){
  pheno.dat$yMu[i] <-  pheno.dat$alphaPhenoSp[i] + pheno.dat$alphaForceSp[i] * pheno.dat$forcei[i]  +   pheno.dat$alphaChillSp[i] * pheno.dat$chilli[i]  +  pheno.dat$alphasalinSp[i] * pheno.dat$salini[i]   + pheno.dat$alphabiomassSp[i] * pheno.dat$biomassi[i] + pheno.dat$alphamacroASp[i] * pheno.dat$macroAi[i] +   pheno.dat$alphabedAreaSp[i] * pheno.dat$bedAreai[i] + pheno.dat$alphadepthSp[i] * pheno.dat$depthi[i] +  pheno.dat$alphadistESp[i] * pheno.dat$distEi[i]
}

pheno.dat$yPhenoi <- pheno.dat$yMu + pheno.dat$ePhen

all.data <- list(
  n_spec = Nspp,
  Nph = nrow(pheno.dat),
  phenology_species = as.numeric(as.factor(pheno.dat$species)),
  species2 = as.numeric(as.factor(pheno.dat$species)),
  yPhenoi = pheno.dat$yPhenoi,
  Var1i = pheno.dat$chilli,
  Var2i = pheno.dat$forcei,
  Var3i = pheno.dat$salini,
  Var4i = pheno.dat$biomassi,
  Var5i = pheno.dat$macroAi,
  Var6i = pheno.dat$bedAreai,
  Var7i = pheno.dat$depthi,
  Var8i = pheno.dat$distEi)
# started running at 1:58 50 s transition time
mdl <- stan("m2m/stan/justDummyInt_2.stan",
            data = all.data,
            iter = 4000, warmup = 3000, chains=4,
            include = FALSE, pars = c("y_hat")
)


sumer <- data.frame(summary(mdl)$summary[c("muVar1Sp","muVar2Sp",
                                           "muVar3Sp","muVar4Sp",
                                           "muVar5Sp",
                                           "muVar6Sp",
                                           "muVar7Sp", "muVar8Sp",
                                           "muPhenoSp",
                                           "sigmaVar1Sp", "sigmaVar2Sp", 
                                           "sigmaVar3Sp",  
                                           "sigmaVar4Sp", 
                                           "sigmaVar5Sp",
                                           "sigmaVar6Sp",
                                           "sigmaVar7Sp",
                                           "sigmaVar8Sp",
                                           "sigmaPhenoSp",
                                           "sigmapheno_y"),])

mdlparam <- data.frame(var = c("muVar1Sp","muVar2Sp",
                               "muVar3Sp", "muVar4Sp",
                               "muVar5Sp",
                               "muVar6Sp",
                               "muVar7Sp",
                               "muVar8Sp",
                               "muPhenoSp",
                               "sigmaVar1Sp",  "sigmaVar2Sp", 
                               "sigmaVar3Sp",  "sigmaVar4Sp", 
                               "sigmaVar5Sp", 
                               "sigmaVar6Sp", 
                               "sigmaVar7Sp", 
                               "sigmaVar8Sp", 
                               "sigmaPhenoSp",
                               "sigmapheno_y"), 
                       param = c(-4, -1, 2, 4,
                                 6, 
                                 8, 10, 12, 
                                 80, 1,1,
                                 2, 3, 4,
                                 5, 6, 7,
                                 10,3))
output <- cbind(mdlparam, sumer); output

# Parameter Test.data.values    Estiamte        X2.5         X55         X50         X75      X97.5
# 1         mutran                5   4.7180702   4.1419345   4.5146800   4.7135841   4.9221097   5.296150
# 2      mutranLat                2   2.2290689   1.5452195   1.9988924   2.2325799   2.4692511   2.892048
# 3     mu_forcesp              -10 -10.0168913 -12.1073613 -10.7438022 -10.0506820  -9.2684616  -7.945515
# 4     mu_chillsp              -14 -14.5300854 -16.9042080 -15.3608747 -14.5675223 -13.6647878 -12.109585
# 5     mu_photosp               -5  -5.4247022  -6.8634576  -5.9516105  -5.4564071  -4.8874566  -3.895768
# 6     mu_phenosp               80  80.3587164  77.9590554  79.5359570  80.3723926  81.1714601  82.732046
# 7   sigma_traity                5   4.9130902   4.8270385   4.8829401   4.9132615   4.9420710   5.002455
# 8       sigma_sp               10  10.7699891   9.2335413  10.1295446  10.7084796  11.3305787  12.674339
# 9  sigma_forcesp                1   1.0349833   0.7823935   0.9384598   1.0308520   1.1249996   1.309711
# 10 sigma_chillsp                1   0.8102007   0.3266061   0.6802449   0.8198023   0.9574059   1.204739
# 11 sigma_photosp                1   0.8482172   0.6871503   0.7842336   0.8428351   0.9076939   1.038829
# 12  sigma_phenoy                3   3.0212577   2.9838130   3.0084186   3.0208992   3.0343502   3.058598
# 13         betaF                3   3.0042213   2.9635365   2.9907854   3.0046815   3.0178851   3.044582
# 14         betaC               -4  -3.9976865  -4.0426522  -4.0135536  -3.9975950  -3.9818976  -3.950773
# 15         betaP               -2  -1.9979220  -2.0258830  -2.0081868  -1.9974009  -1.9876068  -1.971152

# postLMA<- data.frame(rstan::extract(mdl))
# 
# pdf("betaTraitChillPostPrior.pdf")
# hist(postLMA$betaTraitxForce, main ="betaTraitxChill", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
# hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
# abline(v =0, col="red", lwd=3, lty=2)
# dev.off()
# 
# muTraitSp <- sumer[grep("muSp", rownames(sumer))]
# pdf("latitudeMdlEstivsSimNoGrand.pdf")
# plot(muTraitSp ~ mu.trtsp, xlab = "simulated muTraitSp", ylab = "mdl estimated muTraitSp")
# abline(0,1)
# dev.off()
# 