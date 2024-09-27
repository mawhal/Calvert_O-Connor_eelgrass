# started Sept 26, 2024

# Aim of this code is to just play around with Richness

require(Richness)
require(vegan)

data("BCI")

estimateRichness(BCI)

estimateRichness(BCI, meanStates = TRUE, Apx_detectP_terms = TRUE)

estimateRichness(BCI, boot = TRUE, numBoot = 10)


data("pyrifos")
week <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
communitylist <- split(pyrifos, f = week)
estimateRichness(communitylist)
# bootstrapping on a list returns a tidy data frame
require(tidyverse)
temp <- estimateRichness(communitylist, boot = TRUE, numBoot = 10)
