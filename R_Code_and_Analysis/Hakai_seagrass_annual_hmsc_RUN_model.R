####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code will prepare data for running HMSC models 
###      on grazer community data
###  code by Matt Whalen
###  started on   18 February 2020
###  
####################################################################



# load libraries
library(Hmsc)


# read models
models <- list.files( path = getwd(), pattern = "*.Rdata")
# pick a model
mselect <- models[1]
# load it - will bring "x" into the environment, which has everything we need to run a model
load( mselect ) 

### running the model
thin = 100
samples = 1000
transient = 100
nChains = 4 
set.seed(1)

# Run MCMC chains
mod <- sampleMcmc( x, samples = samples , transient = transient, 
                  thin = thin, verbose = 10, nChains = nChains) # nParallel = 3 when running more than 1 Chain for the real model = running more than one chain allows us to check if all arer doing the same thing


# saving
model <- unlist(lapply( strsplit( mselect, "_" ), function(z) paste(z[1],z[2],sep="_") ))
ModelDir <- getwd()
filename = file.path(ModelDir, paste(as.character(model),
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples),".Rdata",sep = ""))

save(mod,file=filename)
