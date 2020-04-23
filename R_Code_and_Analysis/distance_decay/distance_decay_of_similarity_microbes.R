### DISTANCE DECAY OF SIMILARITY ###
### Author: Bianca Trevizan Segovia ###
### Date created: April 7, 2020 ###

library(vegan)
library(stats)
library(ggplot2)
library(readr)
library(reshape2)
library(dplyr)

#############################
############ 16S ############
#############################

### Read table metadata and abundances
microbes_16S <- read.csv("Data/data_parfrey/16S/16S_ASV_MASTER_Hakai_final.csv", header=T)
names(microbes_16S)[1:17]

# split by year
comm.years <- split( microbes_16S, microbes_16S$year )
lapply(comm.years, dim)

# remove taxa that did not appear in a given year
comm.zero <- lapply( comm.years, function(z) {
  tmp <- z[,-c(1:16)]
  w <- colSums(tmp)>0
  df <- tmp[, w ]
  return( data.frame(z[,c(1:16)],df))
})

lat_long <- lapply(comm.zero, function(z) z[,9:10])
d <- lapply(comm.zero, function(z) z[,-c(1:16)])


spe_hell_16S <- mapply( function(d) decostand(d, "hellinger"), 
                d )

### CREATE DISTANCE MATRIX FOR SPECIES
spe_dist_16S <- lapply( spe_hell_16S[c(2,4,6,8)], function(d) vegdist(spe_hell_16S,method="bray"), 
                        d )
lapply( adds[c(2,4,6,8)], function(z) data.frame(z[c("z","pval")]) )



### Spatial data - coordinates lat long
spa_coord_16S <- microbes_16S %>% 
  dplyr::select(c("LAT" , "LONG"))

### CREATE DISTANCE MATRIX FOR LATLONG
spa_dist_16S <- vegdist(spa_coord_16S, method="eu")
spa_dist_16S

## MANTEL
mantel_16S <-mantel(spe_dist_16S, spa_dist_16S, permutations=9999, method = "pearson") 
mantel_16S
out_mantel_16S <- capture.output(mantel_16S)
write.csv(as.data.frame(out_mantel_16S), file = "R_Code_and_Analysis/output_data/mantel_distance_decay_output/mantel_distance_decay_16S.csv", row.names=F)

### convert into data frame for ggplot
df_spe_dist_16S.vector <- melt(as.matrix(spe_dist_16S), varnames = c("row", "col"))
names(df_spe_dist_16S.vector)[names(df_spe_dist_16S.vector) == 'value'] <- 'dissimilarity_16S'

df_spa_dist_16S.vector <- melt(as.matrix(spa_dist_16S), varnames = c("row", "col"))
names(df_spa_dist_16S.vector)[names(df_spa_dist_16S.vector) == 'value'] <- 'spatial_distance_16S'
as.data.frame(df_spe_dist_16S.vector)
as.data.frame(df_spa_dist_16S.vector)

### make a single data frame for dissimilarity, similarity and spatial_distance
all_data_16S <- bind_cols(df_spe_dist_16S.vector, df_spa_dist_16S.vector) # bind_cols combine two data frames

### Transforming into similarities = 1 - Bray
all_data_16S$similarity_16S <- 1 - all_data_16S$dissimilarity_16S
head(all_data_16S)

# Fit regression line
require(stats)
reg<-lm(similarity_16S ~ spatial_distance_16S, data = all_data_16S)
reg #

coeff=coefficients(reg) 
# Equation of the line : ###If you want to add equation to the graph
eq = paste0("y = ", round(coeff[2],2), "*x + ", round(coeff[1],2)) ### print eq to grab equation and put it on my_text.2
# Plot -  use coefficient values from result above
### use this below **** geom_abline(intercept =  0.2947272 , slope =  -0.2071026   , color="red")

# Define and add annotation -------------------------------------
library(grid)
my_text <- bquote(paste("r = 0.324, P< 0.001")) # mantel_16S 
my_grob = grid.text(my_text, x=0.8,  y=0.05, gp=gpar(col="black", fontsize=11, fontface="bold"))

my_text.2 <- bquote(paste("y = 0.29 - 0.21*x")) ###change order to be y = a + bx
my_grob.2 = grid.text(my_text.2, x=0.8,  y=0.1, gp=gpar(col="black", fontsize=11, fontface="bold"))

###Scatterplot with ggplot
p <- ggplot(all_data_16S, aes(x=spatial_distance_16S, y=dissimilarity_16S)) + geom_point(shape=1, size=2.5, color= "palegreen3")+
  #annotate("text", x = 50, y = 1.00, label = "B", size = 6, colour = "black") +
  theme_bw() + 
  theme (axis.title.x = element_text(size=18), #font size of x title
         axis.title.y = element_text(size=18), #font size of y title
         axis.text = element_text(size = 14), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank()) #remove lines outside the graph

p2 <- p + ggtitle("") +
  xlab("spatial distance") + ylab("16S species composition dissimilarity") +    
  geom_abline(intercept =  0.6856377, slope = 0.2252979 , color="red" , linetype = "dashed", size=1.5 ) +
  annotation_custom(my_grob) + annotation_custom(my_grob.2) 
p2

pdf("/Users/bia/PostDoc/projects/Hakai_Quadra_data_retreat/DISTANCE_DECAY/dist_decay_16S.pdf",width=8)
print(p2)
dev.off()

library(betapart) #try DDS with baselga's model exponential 
DDS.betapart.exp <- decay.model(spe_dist_16S, spa_dist_16S, model.type="exponential", y.type="similarities", perm=100)
plot.decay(DDS.betapart.exp, col=rgb(0,0,0,0.5))
plot.decay(DDS.betapart.exp, col="red", remove.dots=TRUE, add=TRUE)



############### ************* ###############
###############     18S       ###############
############### ************* ###############

### Read table metadata and abundances
# *** this one contains mcmullin 2016 data!!!
microbes_18S <- read.csv("~/PostDoc/projects/Hakai_Quadra_data_retreat/18S_retreat/MASTER_18S_DISTANCE_DECAY.csv", header=T)
names(microbes_18S)
nrow(microbes_18S)
#remove rows with NAs
microbes_18S_no_NAs <- microbes_18S %>% drop_na()
nrow(microbes_18S_no_NAs)

# what samples did we loose? 
loose <- anti_join(microbes_18S, microbes_18S_no_NAs)

### get only abundance
abund_18S <- microbes_18S_no_NAs %>% dplyr::select(-(1:10))
names(abund_18S)
#View(abund_18S)
#Hellinger pre-transformation of the species data
spe_hell_18S <- decostand(abund_18S, "hellinger")
spe_hell_18S

### CREATE DISTANCE MATRIX FOR SPECIES
spe_dist_18S <-vegdist(spe_hell_18S,method="bray")
spe_dist_18S
#View(df_spe_dist_18S.vector)

# df_spe_dist_18S<- data.frame(dissimilarity=spe_dist_18S[upper.tri(spe_dist_18S, diag = FALSE)])
# df_spe_dist_18S #get row number that contain values to remove NAs
# df_spe_dist_18S.vector <- df_spe_dist_18S[1:2016,]

### Spatial data - coordinates lat long
spa_coord_18S <- microbes_18S_no_NAs %>% 
  dplyr::select(c("LAT" , "LONG"))

### CREATE DISTANCE MATRIX FOR LATLONG
spa_dist_18S<-vegdist(spa_coord_18S, method="eu")
spa_dist_18S

## MANTEL
mantel_18S <-mantel(spe_dist_18S, spa_dist_18S,permutations=9999, method = "pearson") 
mantel_18S
out_mantel_18S <- capture.output(mantel_18S)
write.csv(as.data.frame(out_mantel_18S), file = "/Users/bia/PostDoc/projects/Hakai_Quadra_data_retreat/DISTANCE_DECAY/out_mantel_distance_decay_18S.csv", row.names=F)

###ggplot only workd with data frames, so we need to convert this data into data frame
df_spe_dist_18S.vector <- melt(as.matrix(spe_dist_18S), varnames = c("row", "col"))
names(df_spe_dist_18S.vector)[names(df_spe_dist_18S.vector) == 'value'] <- 'dissimilarity_18S'
df_spe_dist_18S.vector

df_spa_dist_18S.vector <- melt(as.matrix(spa_dist_18S), varnames = c("row", "col"))
names(df_spa_dist_18S.vector)[names(df_spa_dist_18S.vector) == 'value'] <- 'spatial_distance_18S'
df_spa_dist_18S.vector

as.data.frame(df_spe_dist_18S.vector)
as.data.frame(df_spa_dist_18S.vector)

###make a single data frame for dissimilarity, similarity and spatial_distance
all_data_18S <- bind_cols(df_spe_dist_18S.vector, df_spa_dist_18S.vector) # bind_cols combine two data frames = CAUTION = use only if they are in the same order
attach(all_data_18S)

### Transforming into similarities = 1 - Bray
all_data_18S$similarity_18S <- 1 - all_data_18S$dissimilarity_18S
head(all_data_18S)

# Fit regression line
require(stats)
reg<-lm(dissimilarity_18S ~ spatial_distance_18S, data = all_data_18S)
reg #

coeff=coefficients(reg) 
# Equation of the line : ###If you want to add equation to the graph
eq = paste0("y = ", round(coeff[2],2), "*x + ", round(coeff[1],2)) ### print eq to grab equation and put it on my_text.2
# Plot -  use coefficient values from result above
### use this below **** geom_abline(intercept = 0.7129, slope = 0.3333 , color="red")*** #linetype="dotted" if you want dotted

# Define and add annotation -------------------------------------
library(grid)
my_text <- bquote(paste("r = 0.42, P< 0.001"))
my_grob = grid.text(my_text, x=0.8,  y=0.05, gp=gpar(col="black", fontsize=11, fontface="bold"))
### use this below ***annotation_custom(my_grob)***

my_text.2 <- bquote(paste("y = 0.71 + 0.33*x")) ###change order to be y = a + bx
my_grob.2 = grid.text(my_text.2, x=0.8,  y=0.1, gp=gpar(col="black", fontsize=11, fontface="bold"))
### use this below ***annotation_custom(my_grob)***

###Scatterplot with ggplot
p <- ggplot(all_data_18S, aes(x=spatial_distance_18S, y=dissimilarity_18S)) + geom_point(shape=1, size=2.5, color= "palegreen4")+
  #annotate("text", x = 50, y = 1.00, label = "B", size = 6, colour = "black") +
  theme_bw() + 
  theme (axis.title.x = element_text(size=18), #font size of x title
         axis.title.y = element_text(size=18), #font size of y title
         axis.text = element_text(size = 14), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank()) #remove lines outside the graph

p2 <- p + ggtitle("") +
  xlab("spatial distance") + ylab("18S species composition dissimilarity") +    
  geom_abline(intercept =  0.7129, slope = 0.3333 , color="red" , linetype = "dashed", size=1.5 ) +
  annotation_custom(my_grob) + annotation_custom(my_grob.2) 
p2

pdf("/Users/bia/PostDoc/projects/Hakai_Quadra_data_retreat/DISTANCE_DECAY/dist_decay_18S.pdf",width=8)
print(p2)
dev.off()

library(betapart) #try DDS with baselga's model exponential 
DDS.betapart.exp <- decay.model(spe_dist_18S, spa_dist_18S, model.type="exponential", y.type="similarities", perm=100)
plot.decay(DDS.betapart.exp, col=rgb(0,0,0,0.5))
plot.decay(DDS.betapart.exp, col="red", remove.dots=TRUE, add=TRUE)


###### DIFFERENCES IN SLOPE ######
library(simba)
# zostera vs seawater # results returned positive = first relationship exhibits the steeper slope (zostera steeper)
diffslope_16S_18S <-diffslope(spe_dist_16S, spa_dist_16S, spe_dist_18S, spa_dist_18S,permutations = 1000, ic = FALSE,
                              resc.x = FALSE, resc.y = TRUE, trace=FALSE)

out_diffslope_16S_18S <- capture.output(diffslope_16S_18S)
write.csv(as.data.frame(out_diffslope_16S_18S), file = "/Users/bia/PostDoc/projects/Hakai_Quadra_data_retreat/DISTANCE_DECAY/out_diffslope_16S_18S.csv",row.names=F)


###diffslope Calculate the difference in slope or intercept of two regression lines
#Description
#The function can be used to calculate the difference in slope between two datasets containing each
#two vectors. Follows an idea of Nekola & White (1999) for calculating the statistical inference of the
#difference in slope between two regression lines. The plot method allows easy plotting of the actual difference in slope against the distribution of permuted values.

###diffslope: As the function was initially build to easily calculate the difference in slope between
#the regression lines of distance decay plots, the independent vectors are meant to contain distance
#values whereas the dependent vectors should represent similarity values.For each permutation run the rows are interchanged
#randomly between the two data.frames and the difference in slope calculated thereafter is calculated and collected into a vector. The p-value is
#then computed as the ratio between the number of cases where the differences in slope exceed the
#difference in slope of the inital configuration and the number of permutations.

#If the difference in slope returns negative, the slope (distance decay) of the second relationship is
#less pronounced, if it returns positive, the second relationship exhibits a stronger distance decay
#(slope) than the first. This holds for distance decay relationships. If y increases with x, it is vice versa.
## email answer by simba creator = You are absolutely right, when you have a positive DDR the result returns positive if the first relationship exhibits the steeper slope.

##################################
