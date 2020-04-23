### Plots for HMSC model visualization microbes 16S ###
### Author: Keila Stark ###
### Date created: February, 2020 ###
### Modified by: Bia Segovia ###
### Date modified: February 24, 2020 ###

# load libraries
library(Hmsc)
library(tidyverse)
library(viridis)
library(corrplot)
library(vegan) 

####################################
######### 18S GENUS 2015 ###########
####################################

# load model and check the name in the environment, replace that in the m(models)
load("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_models/model_microbes_18S_2015_genus_chains_4_thin_100_samples_1000.Rdata")


# MCMC trace plot inspection / model fit----------------------------------------------

mpost = convertToCodaObject(mod_2015)
plot( mpost$Beta )

#I realize this is ugly or overwhelming but I really can't think of any other way to plot this


## Assess model fit
preds = computePredictedValues(mod_2015)
MF <- evaluateModelFit(hM = mod_2015, predY = preds)
# The R2 and RMSE 


# Visualizing significance and directionality of beta parameters (I think Whalen wrote this code) ----------

## parameter estimates
postBeta = getPostEstimate(mod_2015, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth_h_h", "area", "biomass",
                              "lai", "microepi", "macro" ),
                            
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth_h", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
                          levels = colnames(mod_2015$Y)[order(colSums(mod_2015$Y),decreasing = TRUE)],
                          ordered = TRUE)

par_est <- ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")


pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/par_est_genus_18S_2015.pdf")
par_est
dev.off()

# Variance partitioning ---------------------------------------------------

VP = computeVariancePartitioning(mod_2015) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP) 

# rowMeans(VP$vals)
# ^ this line gives mean variance explained across species, you could manually write these percentages in to the labels in the new dataframe below so it's easier to see the importance of the variables. The package version automatically prints these, but I don't know how to replicate this in tidyverse

VP.df <- as.data.frame(VP$vals) %>%
  mutate(effect = factor(c("temperature","salinity","depth_h","area",
                           "biomass","LAI","microepiphytes","macroalgae",
                           "quadrat","region","site"),
                         levels = rev(c("temperature","salinity","depth_h","area",
                                        "biomass","LAI","microepiphytes","macroalgae",
                                        "quadrat","region","site")),
                         ordered = TRUE)) %>%
  gather(key = species, value = variance, -effect) %>%
  group_by(species) %>%
  mutate(tempR2 = variance[effect == "temperature"])


hold <- VP.df %>% filter(effect == "temperature") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(mod_2015$Y)[order(colSums(mod_2015$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(mod_2015$Y))


varpart_18S_genus_2015 <- ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "Variance explained")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "Taxa")

ggsave("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/varpart_18S_genus_2015.pdf", plot = varpart_18S_genus_2015 , width=250, height=200, units="mm",dpi=300)

# Example code with different colours (generated from viridis) and filler text if you want to display the values for mean % variance explained:

#par(mar=c(8,8,2.5,13))
#ggplot(VP.df1, aes(y = Variance, x = taxon, fill = factor(effect, levels = c("temperature (x%)","salinity (x%)","depth_h (x%)","area (x%)","biomass (x%)","LAI (x%)","microepiphytes (x%)","macroalgae (x%)","random: quadrat (x%)","random: region (x%)","random: site (x%)"))))+geom_bar(stat = "identity", color = 1)+ theme_classic()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))+ scale_fill_manual(values = c("#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF", "#FCA636FF","#FCA636FF","#FCA636FF","#FCA636FF", "#E16462FF" ,"#B12A90FF" ,"#6A00A8FF" , "#0D0887FF"), name = "Mean variance explained across all species")+ geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+ geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+ scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+  xlab(label = "Species")

# I would suggest picking colours that reflect the different classes of variables (abiotic, biotic, random)


# Correlation plots -------------------------------------------------------

OmegaCor = computeAssociations(mod_2015) # calculate covariance matrix from omega parameter
supportLevel = 0.95 #set 95% credible interval for values so those that don't have high enough statistical confidence print as blank/white

#create covariance matrices from omega parameter estimates

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot2 = ((OmegaCor[[2]]$support>supportLevel)
           + (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean

toPlot3 = ((OmegaCor[[3]]$support>supportLevel)
           + (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean


## These 3 create the corrplots withoout re-ordering species names:

# corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[1]), mar=c(0,0,0,0)) 

# corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[2]), mar=c(0,0,1,0))

# corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[3]), mar=c(0,0,1,0)) 

## These 3 create the corrplots re-ordering according to Ward D2 cluster analysis so co-occurrences groupings come up:

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2015_quadrat.pdf")
corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2015$rLNames[1]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2015_region.pdf")
corrplot_18S_genus_2015_region <-corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2015$rLNames[2]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2015_site.pdf")
corrplot_18S_genus_2015_site <- corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method = "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2015$rLNames[3]), mar=c(0,0,1,0))
dev.off()


####################################
######### 18S FAMILY 2015 ###########
####################################

# load model and check the name in the environment, replace that in the m(models)
load("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_models/model_microbes_18S_2015_family_chains_4_thin_100_samples_1000.Rdata")


# MCMC trace plot inspection / model fit----------------------------------------------

mpost = convertToCodaObject(mod_2015)
plot( mpost$Beta )

#I realize this is ugly or overwhelming but I really can't think of any other way to plot this


## Assess model fit
preds = computePredictedValues(mod_2015)
MF <- evaluateModelFit(hM = mod_2015, predY = preds)
# The R2 and RMSE 


# Visualizing significance and directionality of beta parameters (I think Whalen wrote this code) ----------

## parameter estimates
postBeta = getPostEstimate(mod_2015, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth_h_h", "area", "biomass",
                              "lai", "microepi", "macro" ),
                            
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth_h", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
                          levels = colnames(mod_2015$Y)[order(colSums(mod_2015$Y),decreasing = TRUE)],
                          ordered = TRUE)

par_est <- ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")


pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/par_est_family_18S_2015.pdf")
par_est
dev.off()

# Variance partitioning ---------------------------------------------------

VP = computeVariancePartitioning(mod_2015) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP) 

# rowMeans(VP$vals)
# ^ this line gives mean variance explained across species, you could manually write these percentages in to the labels in the new dataframe below so it's easier to see the importance of the variables. The package version automatically prints these, but I don't know how to replicate this in tidyverse

VP.df <- as.data.frame(VP$vals) %>%
  mutate(effect = factor(c("temperature","salinity","depth_h","area",
                           "biomass","LAI","microepiphytes","macroalgae",
                           "quadrat","region","site"),
                         levels = rev(c("temperature","salinity","depth_h","area",
                                        "biomass","LAI","microepiphytes","macroalgae",
                                        "quadrat","region","site")),
                         ordered = TRUE)) %>%
  gather(key = species, value = variance, -effect) %>%
  group_by(species) %>%
  mutate(tempR2 = variance[effect == "temperature"])


hold <- VP.df %>% filter(effect == "temperature") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(mod_2015$Y)[order(colSums(mod_2015$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(mod_2015$Y))


varpart_18S_family_2015 <- ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "Variance explained")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "Taxa")

ggsave("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/varpart_18S_family_2015.pdf", plot = varpart_18S_family_2015 , width=250, height=200, units="mm",dpi=300)

# Example code with different colours (generated from viridis) and filler text if you want to display the values for mean % variance explained:

#par(mar=c(8,8,2.5,13))
#ggplot(VP.df1, aes(y = Variance, x = taxon, fill = factor(effect, levels = c("temperature (x%)","salinity (x%)","depth_h (x%)","area (x%)","biomass (x%)","LAI (x%)","microepiphytes (x%)","macroalgae (x%)","random: quadrat (x%)","random: region (x%)","random: site (x%)"))))+geom_bar(stat = "identity", color = 1)+ theme_classic()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))+ scale_fill_manual(values = c("#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF", "#FCA636FF","#FCA636FF","#FCA636FF","#FCA636FF", "#E16462FF" ,"#B12A90FF" ,"#6A00A8FF" , "#0D0887FF"), name = "Mean variance explained across all species")+ geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+ geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+ scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+  xlab(label = "Species")

# I would suggest picking colours that reflect the different classes of variables (abiotic, biotic, random)


# Correlation plots -------------------------------------------------------

OmegaCor = computeAssociations(mod_2015) # calculate covariance matrix from omega parameter
supportLevel = 0.95 #set 95% credible interval for values so those that don't have high enough statistical confidence print as blank/white

#create covariance matrices from omega parameter estimates

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot2 = ((OmegaCor[[2]]$support>supportLevel)
           + (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean

toPlot3 = ((OmegaCor[[3]]$support>supportLevel)
           + (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean


## These 3 create the corrplots withoout re-ordering species names:

# corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[1]), mar=c(0,0,0,0)) 

# corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[2]), mar=c(0,0,1,0))

# corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[3]), mar=c(0,0,1,0)) 

## These 3 create the corrplots re-ordering according to Ward D2 cluster analysis so co-occurrences groupings come up:

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2015_quadrat.pdf")
corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2015$rLNames[1]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2015_region.pdf")
corrplot_18S_family_2015_region <-corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2015$rLNames[2]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2015_site.pdf")
corrplot_18S_family_2015_site <- corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method = "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2015$rLNames[3]), mar=c(0,0,1,0))
dev.off()



####################################
######### 18S GENUS 2016 ###########
####################################

# load model and check the name in the environment, replace that in the m(models)
load("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_models/model_microbes_18S_2016_genus_chains_4_thin_100_samples_1000.Rdata")


# MCMC trace plot inspection / model fit----------------------------------------------

mpost = convertToCodaObject(mod_2016)
plot( mpost$Beta )

#I realize this is ugly or overwhelming but I really can't think of any other way to plot this


## Assess model fit
preds = computePredictedValues(mod_2016)
MF <- evaluateModelFit(hM = mod_2016, predY = preds)
# The R2 and RMSE 


# Visualizing significance and directionality of beta parameters (I think Whalen wrote this code) ----------

## parameter estimates
postBeta = getPostEstimate(mod_2016, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth_h_h", "area", "biomass",
                              "lai", "microepi", "macro" ),
                            
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth_h", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
                          levels = colnames(mod_2016$Y)[order(colSums(mod_2016$Y),decreasing = TRUE)],
                          ordered = TRUE)

par_est <- ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")


pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/par_est_genus_18S_2016.pdf")
par_est
dev.off()

# Variance partitioning ---------------------------------------------------

VP = computeVariancePartitioning(mod_2016) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP) 

# rowMeans(VP$vals)
# ^ this line gives mean variance explained across species, you could manually write these percentages in to the labels in the new dataframe below so it's easier to see the importance of the variables. The package version automatically prints these, but I don't know how to replicate this in tidyverse

VP.df <- as.data.frame(VP$vals) %>%
  mutate(effect = factor(c("temperature","salinity","depth_h","area",
                           "biomass","LAI","microepiphytes","macroalgae",
                           "quadrat","region","site"),
                         levels = rev(c("temperature","salinity","depth_h","area",
                                        "biomass","LAI","microepiphytes","macroalgae",
                                        "quadrat","region","site")),
                         ordered = TRUE)) %>%
  gather(key = species, value = variance, -effect) %>%
  group_by(species) %>%
  mutate(tempR2 = variance[effect == "temperature"])


hold <- VP.df %>% filter(effect == "temperature") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(mod_2016$Y)[order(colSums(mod_2016$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(mod_2016$Y))


varpart_18S_genus_2016 <- ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "Variance explained")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "Taxa")

ggsave("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/varpart_18S_genus_2016.pdf", plot = varpart_18S_genus_2016 , width=250, height=200, units="mm",dpi=300)

# Example code with different colours (generated from viridis) and filler text if you want to display the values for mean % variance explained:

#par(mar=c(8,8,2.5,13))
#ggplot(VP.df1, aes(y = Variance, x = taxon, fill = factor(effect, levels = c("temperature (x%)","salinity (x%)","depth_h (x%)","area (x%)","biomass (x%)","LAI (x%)","microepiphytes (x%)","macroalgae (x%)","random: quadrat (x%)","random: region (x%)","random: site (x%)"))))+geom_bar(stat = "identity", color = 1)+ theme_classic()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))+ scale_fill_manual(values = c("#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF", "#FCA636FF","#FCA636FF","#FCA636FF","#FCA636FF", "#E16462FF" ,"#B12A90FF" ,"#6A00A8FF" , "#0D0887FF"), name = "Mean variance explained across all species")+ geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+ geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+ scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+  xlab(label = "Species")

# I would suggest picking colours that reflect the different classes of variables (abiotic, biotic, random)


# Correlation plots -------------------------------------------------------

OmegaCor = computeAssociations(mod_2016) # calculate covariance matrix from omega parameter
supportLevel = 0.95 #set 95% credible interval for values so those that don't have high enough statistical confidence print as blank/white

#create covariance matrices from omega parameter estimates

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot2 = ((OmegaCor[[2]]$support>supportLevel)
           + (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean

toPlot3 = ((OmegaCor[[3]]$support>supportLevel)
           + (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean


## These 3 create the corrplots withoout re-ordering species names:

# corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[1]), mar=c(0,0,0,0)) 

# corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[2]), mar=c(0,0,1,0))

# corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[3]), mar=c(0,0,1,0)) 

## These 3 create the corrplots re-ordering according to Ward D2 cluster analysis so co-occurrences groupings come up:

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2016_quadrat.pdf")
corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2016$rLNames[1]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2016_region.pdf")
corrplot_18S_genus_2016_region <-corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2016$rLNames[2]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2016_site.pdf")
corrplot_18S_genus_2016_site <- corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method = "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2016$rLNames[3]), mar=c(0,0,1,0))
dev.off()



####################################
######### 18S FAMILY 2016 ###########
####################################

# load model and check the name in the environment, replace that in the m(models)
load("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_models/model_microbes_18S_2016_family_chains_4_thin_100_samples_1000.Rdata")


# MCMC trace plot inspection / model fit----------------------------------------------

mpost = convertToCodaObject(mod_2016)
plot( mpost$Beta )

#I realize this is ugly or overwhelming but I really can't think of any other way to plot this


## Assess model fit
preds = computePredictedValues(mod_2016)
MF <- evaluateModelFit(hM = mod_2016, predY = preds)
# The R2 and RMSE 


# Visualizing significance and directionality of beta parameters (I think Whalen wrote this code) ----------

## parameter estimates
postBeta = getPostEstimate(mod_2016, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth_h_h", "area", "biomass",
                              "lai", "microepi", "macro" ),
                            
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth_h", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
                          levels = colnames(mod_2016$Y)[order(colSums(mod_2016$Y),decreasing = TRUE)],
                          ordered = TRUE)

par_est <- ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")


pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/par_est_family_18S_2016.pdf")
par_est
dev.off()

# Variance partitioning ---------------------------------------------------

VP = computeVariancePartitioning(mod_2016) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP) 

# rowMeans(VP$vals)
# ^ this line gives mean variance explained across species, you could manually write these percentages in to the labels in the new dataframe below so it's easier to see the importance of the variables. The package version automatically prints these, but I don't know how to replicate this in tidyverse

VP.df <- as.data.frame(VP$vals) %>%
  mutate(effect = factor(c("temperature","salinity","depth_h","area",
                           "biomass","LAI","microepiphytes","macroalgae",
                           "quadrat","region","site"),
                         levels = rev(c("temperature","salinity","depth_h","area",
                                        "biomass","LAI","microepiphytes","macroalgae",
                                        "quadrat","region","site")),
                         ordered = TRUE)) %>%
  gather(key = species, value = variance, -effect) %>%
  group_by(species) %>%
  mutate(tempR2 = variance[effect == "temperature"])


hold <- VP.df %>% filter(effect == "temperature") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(mod_2016$Y)[order(colSums(mod_2016$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(mod_2016$Y))


varpart_18S_family_2016 <- ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "Variance explained")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "Taxa")

ggsave("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/varpart_18S_family_2016.pdf", plot = varpart_18S_family_2016 , width=250, height=200, units="mm",dpi=300)

# Example code with different colours (generated from viridis) and filler text if you want to display the values for mean % variance explained:

#par(mar=c(8,8,2.5,13))
#ggplot(VP.df1, aes(y = Variance, x = taxon, fill = factor(effect, levels = c("temperature (x%)","salinity (x%)","depth_h (x%)","area (x%)","biomass (x%)","LAI (x%)","microepiphytes (x%)","macroalgae (x%)","random: quadrat (x%)","random: region (x%)","random: site (x%)"))))+geom_bar(stat = "identity", color = 1)+ theme_classic()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))+ scale_fill_manual(values = c("#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF", "#FCA636FF","#FCA636FF","#FCA636FF","#FCA636FF", "#E16462FF" ,"#B12A90FF" ,"#6A00A8FF" , "#0D0887FF"), name = "Mean variance explained across all species")+ geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+ geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+ scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+  xlab(label = "Species")

# I would suggest picking colours that reflect the different classes of variables (abiotic, biotic, random)


# Correlation plots -------------------------------------------------------

OmegaCor = computeAssociations(mod_2016) # calculate covariance matrix from omega parameter
supportLevel = 0.95 #set 95% credible interval for values so those that don't have high enough statistical confidence print as blank/white

#create covariance matrices from omega parameter estimates

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot2 = ((OmegaCor[[2]]$support>supportLevel)
           + (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean

toPlot3 = ((OmegaCor[[3]]$support>supportLevel)
           + (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean


## These 3 create the corrplots withoout re-ordering species names:

# corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[1]), mar=c(0,0,0,0)) 

# corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[2]), mar=c(0,0,1,0))

# corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[3]), mar=c(0,0,1,0)) 

## These 3 create the corrplots re-ordering according to Ward D2 cluster analysis so co-occurrences groupings come up:

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2016_quadrat.pdf")
corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2016$rLNames[1]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2016_region.pdf")
corrplot_18S_family_2016_region <-corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2016$rLNames[2]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2016_site.pdf")
corrplot_18S_family_2016_site <- corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method = "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2016$rLNames[3]), mar=c(0,0,1,0))
dev.off()




####################################
######### 18S GENUS 2018 ###########
####################################

# load model and check the name in the environment, replace that in the m(models)
load("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_models/model_microbes_18S_2018_genus_chains_4_thin_100_samples_1000.Rdata")


# MCMC trace plot inspection / model fit----------------------------------------------

mpost = convertToCodaObject(mod_2018)
plot( mpost$Beta )

#I realize this is ugly or overwhelming but I really can't think of any other way to plot this


## Assess model fit
preds = computePredictedValues(mod_2018)
MF <- evaluateModelFit(hM = mod_2018, predY = preds)
# The R2 and RMSE 


# Visualizing significance and directionality of beta parameters (I think Whalen wrote this code) ----------

## parameter estimates
postBeta = getPostEstimate(mod_2018, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth_h_h", "area", "biomass",
                              "lai", "microepi", "macro" ),
                            
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth_h", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
                          levels = colnames(mod_2018$Y)[order(colSums(mod_2018$Y),decreasing = TRUE)],
                          ordered = TRUE)

par_est <- ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")


pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/par_est_genus_18S_2018.pdf")
par_est
dev.off()

# Variance partitioning ---------------------------------------------------

VP = computeVariancePartitioning(mod_2018) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP) 

# rowMeans(VP$vals)
# ^ this line gives mean variance explained across species, you could manually write these percentages in to the labels in the new dataframe below so it's easier to see the importance of the variables. The package version automatically prints these, but I don't know how to replicate this in tidyverse

VP.df <- as.data.frame(VP$vals) %>%
  mutate(effect = factor(c("temperature","salinity","depth_h","area",
                           "biomass","LAI","microepiphytes","macroalgae",
                           "quadrat","region","site"),
                         levels = rev(c("temperature","salinity","depth_h","area",
                                        "biomass","LAI","microepiphytes","macroalgae",
                                        "quadrat","region","site")),
                         ordered = TRUE)) %>%
  gather(key = species, value = variance, -effect) %>%
  group_by(species) %>%
  mutate(tempR2 = variance[effect == "temperature"])


hold <- VP.df %>% filter(effect == "temperature") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(mod_2018$Y)[order(colSums(mod_2018$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(mod_2018$Y))


varpart_18S_genus_2018 <- ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "Variance explained")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "Taxa")

ggsave("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/varpart_18S_genus_2018.pdf", plot = varpart_18S_genus_2018 , width=250, height=200, units="mm",dpi=300)

# Example code with different colours (generated from viridis) and filler text if you want to display the values for mean % variance explained:

#par(mar=c(8,8,2.5,13))
#ggplot(VP.df1, aes(y = Variance, x = taxon, fill = factor(effect, levels = c("temperature (x%)","salinity (x%)","depth_h (x%)","area (x%)","biomass (x%)","LAI (x%)","microepiphytes (x%)","macroalgae (x%)","random: quadrat (x%)","random: region (x%)","random: site (x%)"))))+geom_bar(stat = "identity", color = 1)+ theme_classic()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))+ scale_fill_manual(values = c("#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF", "#FCA636FF","#FCA636FF","#FCA636FF","#FCA636FF", "#E16462FF" ,"#B12A90FF" ,"#6A00A8FF" , "#0D0887FF"), name = "Mean variance explained across all species")+ geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+ geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+ scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+  xlab(label = "Species")

# I would suggest picking colours that reflect the different classes of variables (abiotic, biotic, random)


# Correlation plots -------------------------------------------------------

OmegaCor = computeAssociations(mod_2018) # calculate covariance matrix from omega parameter
supportLevel = 0.95 #set 95% credible interval for values so those that don't have high enough statistical confidence print as blank/white

#create covariance matrices from omega parameter estimates

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot2 = ((OmegaCor[[2]]$support>supportLevel)
           + (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean

toPlot3 = ((OmegaCor[[3]]$support>supportLevel)
           + (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean


## These 3 create the corrplots withoout re-ordering species names:

# corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[1]), mar=c(0,0,0,0)) 

# corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[2]), mar=c(0,0,1,0))

# corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[3]), mar=c(0,0,1,0)) 

## These 3 create the corrplots re-ordering according to Ward D2 cluster analysis so co-occurrences groupings come up:

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2018_quadrat.pdf")
corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2018$rLNames[1]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2018_region.pdf")
corrplot_18S_genus_2018_region <-corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2018$rLNames[2]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_genus_2018_site.pdf")
corrplot_18S_genus_2018_site <- corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method = "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2018$rLNames[3]), mar=c(0,0,1,0))
dev.off()




####################################
######### 18S FAMILY 2018 ###########
####################################

# load model and check the name in the environment, replace that in the m(models)
load("~/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_models/model_microbes_18S_2018_family_chains_4_thin_100_samples_1000.Rdata")


# MCMC trace plot inspection / model fit----------------------------------------------

mpost = convertToCodaObject(mod_2018)
plot( mpost$Beta )

#I realize this is ugly or overwhelming but I really can't think of any other way to plot this


## Assess model fit
preds = computePredictedValues(mod_2018)
MF <- evaluateModelFit(hM = mod_2018, predY = preds)
# The R2 and RMSE 


# Visualizing significance and directionality of beta parameters (I think Whalen wrote this code) ----------

## parameter estimates
postBeta = getPostEstimate(mod_2018, parName = "Beta")

pos.neg <- data.frame(pos = c(postBeta$support), neg = c(postBeta$supportNeg))
pos.neg[pos.neg< 0.95] <- 0
pos.neg$neg <- -pos.neg$neg
pos.neg$value <- pos.neg$pos + pos.neg$neg
pos.neg$parameter <- factor(c("intercept", "temperature",  "salinity",
                              "depth_h_h", "area", "biomass",
                              "lai", "microepi", "macro" ),
                            
                            levels = c("intercept", "temperature",  "salinity",
                                       "depth_h", "area", "biomass",
                                       "lai", "microepi", "macro" ),
                            ordered = TRUE)
pos.neg$species <- factor(rep(colnames(postBeta$mean), each = 9),
                          levels = colnames(mod_2018$Y)[order(colSums(mod_2018$Y),decreasing = TRUE)],
                          ordered = TRUE)

par_est <- ggplot(pos.neg, aes(y = parameter, x = species, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", "high" = "red", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 90,size=7.5))+
  xlab(label = "")+
  ylab(label = "")


pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/par_est_family_18S_2018.pdf")
par_est
dev.off()

# Variance partitioning ---------------------------------------------------

VP = computeVariancePartitioning(mod_2018) #, group = c(1,1,1,2,2,3,4,4),groupnames=c("temperature","dispersal","week", "dispersal * week"))
# plotVariancePartitioning(m, VP = VP) 

# rowMeans(VP$vals)
# ^ this line gives mean variance explained across species, you could manually write these percentages in to the labels in the new dataframe below so it's easier to see the importance of the variables. The package version automatically prints these, but I don't know how to replicate this in tidyverse

VP.df <- as.data.frame(VP$vals) %>%
  mutate(effect = factor(c("temperature","salinity","depth_h","area",
                           "biomass","LAI","microepiphytes","macroalgae",
                           "quadrat","region","site"),
                         levels = rev(c("temperature","salinity","depth_h","area",
                                        "biomass","LAI","microepiphytes","macroalgae",
                                        "quadrat","region","site")),
                         ordered = TRUE)) %>%
  gather(key = species, value = variance, -effect) %>%
  group_by(species) %>%
  mutate(tempR2 = variance[effect == "temperature"])


hold <- VP.df %>% filter(effect == "temperature") %>% arrange(desc(tempR2))

VP.df$species <- factor(VP.df$species, 
                        levels = colnames(mod_2018$Y)[order(colSums(mod_2018$Y),decreasing = TRUE)], 
                        ordered = TRUE)

R2.df <- data.frame(R2 = round(MF$SR2,1), species = colnames(mod_2018$Y))


varpart_18S_family_2018 <- ggplot(VP.df,aes(y = variance, x = species, fill = effect))+
  geom_bar(stat = "identity", color = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("darkred", "maroon","pink", "darkorange2","orange","darkgreen","forestgreen","lightgreen", "chartreuse", "gray25","gray75","whitesmoke"), name = "Variance explained")+
  geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+
  geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+
  scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+
  xlab(label = "Taxa")

ggsave("/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/varpart_18S_family_2018.pdf", plot = varpart_18S_family_2018 , width=250, height=200, units="mm",dpi=300)

# Example code with different colours (generated from viridis) and filler text if you want to display the values for mean % variance explained:

#par(mar=c(8,8,2.5,13))
#ggplot(VP.df1, aes(y = Variance, x = taxon, fill = factor(effect, levels = c("temperature (x%)","salinity (x%)","depth_h (x%)","area (x%)","biomass (x%)","LAI (x%)","microepiphytes (x%)","macroalgae (x%)","random: quadrat (x%)","random: region (x%)","random: site (x%)"))))+geom_bar(stat = "identity", color = 1)+ theme_classic()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))+ scale_fill_manual(values = c("#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF", "#FCA636FF","#FCA636FF","#FCA636FF","#FCA636FF", "#E16462FF" ,"#B12A90FF" ,"#6A00A8FF" , "#0D0887FF"), name = "Mean variance explained across all species")+ geom_text(data = R2.df, aes(y = -0.02, fill = NULL, label = R2), size = 2)+ geom_point(data = R2.df, aes(y = -0.06, fill = NULL, size = R2))+ scale_size_continuous(breaks = seq(0.15,0.60,by = 0.15))+  xlab(label = "Species")

# I would suggest picking colours that reflect the different classes of variables (abiotic, biotic, random)


# Correlation plots -------------------------------------------------------

OmegaCor = computeAssociations(mod_2018) # calculate covariance matrix from omega parameter
supportLevel = 0.95 #set 95% credible interval for values so those that don't have high enough statistical confidence print as blank/white

#create covariance matrices from omega parameter estimates

toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot2 = ((OmegaCor[[2]]$support>supportLevel)
           + (OmegaCor[[2]]$support<(1-supportLevel))>0)*OmegaCor[[2]]$mean

toPlot3 = ((OmegaCor[[3]]$support>supportLevel)
           + (OmegaCor[[3]]$support<(1-supportLevel))>0)*OmegaCor[[3]]$mean


## These 3 create the corrplots withoout re-ordering species names:

# corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[1]), mar=c(0,0,0,0)) 

# corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200),  tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[2]), mar=c(0,0,1,0))

# corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), tl.cex=.75, tl.col="black", title=paste("random effect level:", mod$rLNames[3]), mar=c(0,0,1,0)) 

## These 3 create the corrplots re-ordering according to Ward D2 cluster analysis so co-occurrences groupings come up:

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2018_quadrat.pdf")
corrplot(toPlot, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2018$rLNames[1]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2018_region.pdf")
corrplot_18S_family_2018_region <-corrplot(toPlot2, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method= "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2018$rLNames[2]), mar=c(0,0,1,0))
dev.off()

pdf(file="/Users/bia/PostDoc/projects/Calvert_O-Connor_eelgrass/R_Code_and_Analysis/output_data/hmsc_models_microbes/18S_figures/corrplot_18S_family_2018_site.pdf")
corrplot_18S_family_2018_site <- corrplot(toPlot3, method = "color",col=colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust" , hclust.method = "ward.D2", tl.cex=.75, tl.col="black", title=paste("random effect level:", mod_2018$rLNames[3]), mar=c(0,0,1,0))
dev.off()