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
models <- list.files( path = paste0(getwd(),"/output_data"), pattern = "*specification.Rdata")
# pick a model
mselect <- models[1]
# load it - will bring "x" into the environment, which has everything we need to run a model
load( paste0("output_data/", mselect) ) 

### running the model
thin = 10
samples = 100
transient = 10
nChains = 1 
nParallel = nChains
set.seed(1)

# Run MCMC chains
mod <- sampleMcmc( x, samples = samples , transient = transient, 
                  thin = thin, verbose = 10, nChains = nChains,
                  nParallel = nParallel ) 


# saving
model <- unlist(lapply( strsplit( mselect, "_" ), function(z) paste(z[1],z[2],sep="_") ))
ModelDir <- paste0(getwd(),"/output_data")
filename = file.path(ModelDir, paste(as.character(model),
                                     "_chains_",as.character(nChains),
                                     "_thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples),".Rdata",sep = ""))

save(mod,file=filename)
