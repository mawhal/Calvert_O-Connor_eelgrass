####################################################################
###  Hakai + O'Connor Seagrass Mesograzer Data From Calvert Island 
###  Dataset begins in 2014 
###  Samples from several sites (not always the same)
###
###  This code calculates two dimensional locations from oceanographic
###     distance
###
###  code by Matt Whalen
###  started on   04 February 2020
###  
####################################################################

# load libraries
library( tidyverse )


## Note, this script does not calculate oceangraphic distance
## rather, it takes these distances and tries to coerce them 
## into a form that can be used by HMSC (x and y components)

## start with lat/long data - 
latlong <- read_csv( "../metadata/coreen_hakai_sites_location_2015.csv" ) %>% select( lat,long )

# calculate a distance matrix (Euclidean)
dists <- dist(latlong)

# use distance matrix in eigenvalue decomposition
decomp <- eigen(dists)

plot(latlong)
plot(decomp$vectors[,1:2])


gram <- function(X){
  dims = dim(X)
  M = matrix(nrow=dims[1],ncol=dims[2])
  
}