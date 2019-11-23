## Test HMSC run at Seagrass Retreat
# started 23 November 2019

# load libraries
library(Hmsc)
library(tidyverse)

### read data
# response data
Ygrazer
Ymicrobe
Yfish

# explanatory data
X
# filter rows for each response

# random effects
# quad, site, space? Do we need quad and site if we have space?



### building the model

## formula for fixed effects
XFormula = ~()

mgrazer <- Hmsc( Y = Ygrazer, 
           XData = Xgrazer, XFormula = XFormula,
           distr = "poisson",
           studyDesign = studyDesign, 
           ranLevels = list( space=rLspace ) )


### running the model
thin = 1
samples = 100
nChains = 2
set.seed(1)
