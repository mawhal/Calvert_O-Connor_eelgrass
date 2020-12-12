### Taxa shared across regions ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 03, 2020 ###

#### load packages ####
library(dplyr)
library(phyloseq)
library(tidyverse)
library(ggvenn)

########################################
############ ACROSS REGIONS ############
########################################

###############################################
############ 16S prokaryotes GENUS ############
###############################################

###load phyloseq object 
prokary_data_notrarefied <- readRDS("Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")

prokary_data_notrarefied <- prokary_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
prokary.shared <- prokary_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(prokary.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(prokary.shared)[otu >= 1] <- 1 

### Select seawater data ###
prokary.shared_16_17 <- prokary.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
prokary_groupC = subset_samples(prokary.shared_16_17, region == "choked")
prokary_groupG = subset_samples(prokary.shared_16_17, region == "goose")
prokary_groupP = subset_samples(prokary.shared_16_17, region == "pruth")
prokary_groupT = subset_samples(prokary.shared_16_17, region == "triquet")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
prokary_taxa_sums_grC <- as.data.frame(filter_taxa(prokary_groupC, 
                                                 function(x) sum(x))) 
# do sums from presence absence OTU table
prokary_taxa_sums_grG <- as.data.frame(filter_taxa(prokary_groupG, 
                                                 function(x) sum(x))) 
# do sums from presence absence OTU table
prokary_taxa_sums_grP <- as.data.frame(filter_taxa(prokary_groupP, 
                                                 function(x) sum(x))) 
# do sums from presence absence OTU table
prokary_taxa_sums_grT <- as.data.frame(filter_taxa(prokary_groupT, 
                                                 function(x) sum(x))) 

### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_prokary_C <- row.names(
  prokary_taxa_sums_grC)[which(prokary_taxa_sums_grC[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_prokary_G <- row.names(
  prokary_taxa_sums_grG)[which(prokary_taxa_sums_grG[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_prokary_P <- row.names(
  prokary_taxa_sums_grP)[which(prokary_taxa_sums_grP[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_prokary_T <- row.names(
  prokary_taxa_sums_grT)[which(prokary_taxa_sums_grT[,1] >= 2)] 

### Venn Diagram ###
prokary_regions <- list("Choked"=keep_prokary_C,
                    "Goose"=keep_prokary_G,
                    "Pruth"=keep_prokary_P,
                    "Triquet"=keep_prokary_T)
venn_prokary <- ggvenn(
  prokary_regions,
  fill_color = c("#e41a1c","#377eb8","#4daf4a", "#984ea3"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_prokary <- venn_prokary +
  annotate("text", x = 0, y = 1.6, label = "Prokaryotes", size = 8)

venn_prokary

ggsave("R_Code_and_Analysis/varpart/Prokaryotes_taxa_plot.tiff", plot = venn_prokary, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")


##################################################
############ 18S microeukaryotes GENUS ############
##################################################

###load phyloseq object 
microeuk_data_notrarefied <- readRDS("Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds")
#View(as.data.frame(tax_table(microeuk_data_notrarefied)))

microeuk_data_notrarefied <- microeuk_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
microeuk.shared <- microeuk_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(microeuk.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(microeuk.shared)[otu >= 1] <- 1 

### Select seawater data ###
microeuk.shared_16_17 <- microeuk.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
microeuk_groupC = subset_samples(microeuk.shared_16_17, region == "choked")
microeuk_groupG = subset_samples(microeuk.shared_16_17, region == "goose")
microeuk_groupP = subset_samples(microeuk.shared_16_17, region == "pruth")
microeuk_groupT = subset_samples(microeuk.shared_16_17, region == "triquet")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
microeuk_taxa_sums_grC <- as.data.frame(filter_taxa(microeuk_groupC, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
microeuk_taxa_sums_grG <- as.data.frame(filter_taxa(microeuk_groupG, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
microeuk_taxa_sums_grP <- as.data.frame(filter_taxa(microeuk_groupP, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
microeuk_taxa_sums_grT <- as.data.frame(filter_taxa(microeuk_groupT, 
                                                   function(x) sum(x))) 

### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_microeuk_C <- row.names(
  microeuk_taxa_sums_grC)[which(microeuk_taxa_sums_grC[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_microeuk_G <- row.names(
  microeuk_taxa_sums_grG)[which(microeuk_taxa_sums_grG[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_microeuk_P <- row.names(
  microeuk_taxa_sums_grP)[which(microeuk_taxa_sums_grP[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_microeuk_T <- row.names(
  microeuk_taxa_sums_grT)[which(microeuk_taxa_sums_grT[,1] >= 2)] 

### Venn Diagram ###
microeuk_regions <- list("Choked"=keep_microeuk_C,
                        "Goose"=keep_microeuk_G,
                        "Pruth"=keep_microeuk_P,
                        "Triquet"=keep_microeuk_T)
venn_microeuk <- ggvenn(
  microeuk_regions,
  fill_color = c("#e41a1c","#377eb8","#4daf4a", "#984ea3"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_microeuk <- venn_microeuk +
  annotate("text", x = 0, y = 1.6, label = "Microeukaryotes", size = 8)

venn_microeuk

ggsave("R_Code_and_Analysis/varpart/microeukotes_taxa_plot.tiff", plot = venn_microeuk, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")

#############################
############ Inverts ########
#############################

###load phyloseq object 
inverts_data_notrarefied <- readRDS("Data/macro_eukaryotes/invert_phyloseq.rds")

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
inverts.shared <- inverts_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(inverts.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(inverts.shared)[otu >= 1] <- 1 

### Select seawater data ###
inverts.shared_16_17 <- inverts.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
inverts_groupC = subset_samples(inverts.shared_16_17, region == "choked")
inverts_groupG = subset_samples(inverts.shared_16_17, region == "goose")
inverts_groupP = subset_samples(inverts.shared_16_17, region == "pruth")
inverts_groupT = subset_samples(inverts.shared_16_17, region == "triquet")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
inverts_taxa_sums_grC <- as.data.frame(filter_taxa(inverts_groupC, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
inverts_taxa_sums_grG <- as.data.frame(filter_taxa(inverts_groupG, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
inverts_taxa_sums_grP <- as.data.frame(filter_taxa(inverts_groupP, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
inverts_taxa_sums_grT <- as.data.frame(filter_taxa(inverts_groupT, 
                                                   function(x) sum(x))) 

### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_inverts_C <- row.names(
  inverts_taxa_sums_grC)[which(inverts_taxa_sums_grC[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_inverts_G <- row.names(
  inverts_taxa_sums_grG)[which(inverts_taxa_sums_grG[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_inverts_P <- row.names(
  inverts_taxa_sums_grP)[which(inverts_taxa_sums_grP[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_inverts_T <- row.names(
  inverts_taxa_sums_grT)[which(inverts_taxa_sums_grT[,1] >= 2)] 

### Venn Diagram ###
inverts_regions <- list("Choked"=keep_inverts_C,
                        "Goose"=keep_inverts_G,
                        "Pruth"=keep_inverts_P,
                        "Triquet"=keep_inverts_T)
venn_inverts <- ggvenn(
  inverts_regions,
  fill_color = c("#e41a1c","#377eb8","#4daf4a", "#984ea3"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_inverts <- venn_inverts +
  annotate("text", x = 0, y = 1.6, label = "Inverts", size = 8)

venn_inverts

ggsave("R_Code_and_Analysis/varpart/inverts_taxa_plot.tiff", plot = venn_inverts, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")



#################################
############ Macroeuk18S ########
#################################

###load phyloseq object 
macro18S_data_notrarefied <- readRDS("Data/micro_eukaryotes/macro18S_filtered.rds")

macro18S_data_notrarefied <- macro18S_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
macro18S.shared <- macro18S_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(macro18S.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(macro18S.shared)[otu >= 1] <- 1 

### Select seawater data ###
macro18S.shared_16_17 <- macro18S.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
macro18S_groupC = subset_samples(macro18S.shared_16_17, region == "choked")
macro18S_groupG = subset_samples(macro18S.shared_16_17, region == "goose")
macro18S_groupP = subset_samples(macro18S.shared_16_17, region == "pruth")
macro18S_groupT = subset_samples(macro18S.shared_16_17, region == "triquet")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
macro18S_taxa_sums_grC <- as.data.frame(filter_taxa(macro18S_groupC, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
macro18S_taxa_sums_grG <- as.data.frame(filter_taxa(macro18S_groupG, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
macro18S_taxa_sums_grP <- as.data.frame(filter_taxa(macro18S_groupP, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
macro18S_taxa_sums_grT <- as.data.frame(filter_taxa(macro18S_groupT, 
                                                   function(x) sum(x))) 

### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_macro18S_C <- row.names(
  macro18S_taxa_sums_grC)[which(macro18S_taxa_sums_grC[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_macro18S_G <- row.names(
  macro18S_taxa_sums_grG)[which(macro18S_taxa_sums_grG[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_macro18S_P <- row.names(
  macro18S_taxa_sums_grP)[which(macro18S_taxa_sums_grP[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_macro18S_T <- row.names(
  macro18S_taxa_sums_grT)[which(macro18S_taxa_sums_grT[,1] >= 2)] 

### Venn Diagram ###
macro18S_regions <- list("Choked"=keep_macro18S_C,
                        "Goose"=keep_macro18S_G,
                        "Pruth"=keep_macro18S_P,
                        "Triquet"=keep_macro18S_T)
venn_macro18S <- ggvenn(
  macro18S_regions,
  fill_color = c("#e41a1c","#377eb8","#4daf4a", "#984ea3"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_macro18S <- venn_macro18S +
  annotate("text", x = 0, y = 1.6, label = "Macro18S", size = 8)

venn_macro18S

ggsave("R_Code_and_Analysis/varpart/macro18S_taxa_plot.tiff", plot = venn_macro18S, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")

######################################
############ ACROSS YEARS ############
######################################


#########################################
############ 16S prokaryotes GENUS ############
#########################################

###load phyloseq object 
prokary_data_notrarefied <- readRDS("Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")

prokary_data_notrarefied <- prokary_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
prokary.shared <- prokary_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(prokary.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(prokary.shared)[otu >= 1] <- 1 

### Select seawater data ###
prokary.shared_16_17 <- prokary.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
prokary_group16 = subset_samples(prokary.shared_16_17, year == "2016")
prokary_group17 = subset_samples(prokary.shared_16_17, year == "2017")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
prokary_taxa_sums_16 <- as.data.frame(filter_taxa(prokary_group16, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
prokary_taxa_sums_17 <- as.data.frame(filter_taxa(prokary_group17, 
                                                   function(x) sum(x))) 
### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_prokary_16 <- row.names(
  prokary_taxa_sums_16)[which(prokary_taxa_sums_16[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_prokary_17 <- row.names(
  prokary_taxa_sums_17)[which(prokary_taxa_sums_17[,1] >= 2)] 

### Venn Diagram ###
prokary_years <- list("2016"=keep_prokary_16,
                      "2017"=keep_prokary_17)
venn_prokary_years <- ggvenn(
  prokary_years,
  fill_color = c("#e41a1c","#377eb8"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_prokary_years <- venn_prokary_years +
  annotate("text", x = 0, y = 1.6, label = "Prokaryotes_16_17", size = 8)

venn_prokary_years

ggsave("R_Code_and_Analysis/varpart/Prokaryotes_16_17_taxa_plot.tiff", plot = venn_prokary_years, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")


#########################################
############ 18S microeukaryotes ########
#########################################

###load phyloseq object 
microeuk_data_notrarefied <- readRDS("Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds")

microeuk_data_notrarefied <- microeuk_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
microeuk.shared <- microeuk_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(microeuk.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(microeuk.shared)[otu >= 1] <- 1 

### Select seawater data ###
microeuk.shared_16_17 <- microeuk.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
microeuk_group16 = subset_samples(microeuk.shared_16_17, year == "2016")
microeuk_group17 = subset_samples(microeuk.shared_16_17, year == "2017")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
microeuk_taxa_sums_16 <- as.data.frame(filter_taxa(microeuk_group16, 
                                                  function(x) sum(x))) 
# do sums from presence absence OTU table
microeuk_taxa_sums_17 <- as.data.frame(filter_taxa(microeuk_group17, 
                                                  function(x) sum(x))) 
### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_microeuk_16 <- row.names(
  microeuk_taxa_sums_16)[which(microeuk_taxa_sums_16[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_microeuk_17 <- row.names(
  microeuk_taxa_sums_17)[which(microeuk_taxa_sums_17[,1] >= 2)] 

### Venn Diagram ###
microeuk_years <- list("2016"=keep_microeuk_16,
                      "2017"=keep_microeuk_17)
venn_microeuk_years <- ggvenn(
  microeuk_years,
  fill_color = c("#e41a1c","#377eb8"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_microeuk_years <- venn_microeuk_years +
  annotate("text", x = 0, y = 1.6, label = "Microeukaryotes_16_17", size = 8)

venn_microeuk_years

ggsave("R_Code_and_Analysis/varpart/microeukaryotes_16_17_taxa_plot.tiff", plot = venn_microeuk_years, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")


#############################
############ Inverts ########
#############################

###load phyloseq object 
inverts_data_notrarefied <- readRDS("Data/macro_eukaryotes/invert_phyloseq.rds")

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
inverts.shared <- inverts_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(inverts.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(inverts.shared)[otu >= 1] <- 1 

### Select seawater data ###
inverts.shared_16_17 <- inverts.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
inverts_group16 = subset_samples(inverts.shared_16_17, year == "2016")
inverts_group17 = subset_samples(inverts.shared_16_17, year == "2017")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
inverts_taxa_sums_16 <- as.data.frame(filter_taxa(inverts_group16, 
                                                   function(x) sum(x))) 
# do sums from presence absence OTU table
inverts_taxa_sums_17 <- as.data.frame(filter_taxa(inverts_group17, 
                                                   function(x) sum(x))) 
### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_inverts_16 <- row.names(
  inverts_taxa_sums_16)[which(inverts_taxa_sums_16[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_inverts_17 <- row.names(
  inverts_taxa_sums_17)[which(inverts_taxa_sums_17[,1] >= 2)] 

### Venn Diagram ###
inverts_years <- list("2016"=keep_inverts_16,
                       "2017"=keep_inverts_17)
venn_inverts_years <- ggvenn(
  inverts_years,
  fill_color = c("#e41a1c","#377eb8"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_inverts_years <- venn_inverts_years +
  annotate("text", x = 0, y = 1.6, label = "Inverts_16_17", size = 8)

venn_inverts_years

ggsave("R_Code_and_Analysis/varpart/Inverts_16_17_taxa_plot.tiff", plot = venn_inverts_years, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")


#################################
############ Macroeuk18S ########
#################################

###load phyloseq object 
macro18S_data_notrarefied <- readRDS("Data/micro_eukaryotes/macro18S_filtered.rds")

#### subset taxa based on pres/abs threshold per sample group ####
# make presence absence table
macro18S.shared <- macro18S_data_notrarefied # duplicate raw counts phyloseq object
otu <- as.data.frame(otu_table(macro18S.shared)) #get OTU table
#set all positive values in OTU table of project_data.pres_abs to '1'
otu_table(macro18S.shared)[otu >= 1] <- 1 

### Select seawater data ###
macro18S.shared_16_17 <- macro18S.shared %>% subset_samples(year == "2016" | year == "2017" ) 

### subset based on groups you're interested in (regions) ###
macro18S_group16 = subset_samples(macro18S.shared_16_17, year == "2016")
macro18S_group17 = subset_samples(macro18S.shared_16_17, year == "2017")

### remove all OTUs not found at threshold (N samples) ###
# do sums from presence absence OTU table
macro18S_taxa_sums_16 <- as.data.frame(filter_taxa(macro18S_group16, 
                                                  function(x) sum(x))) 
# do sums from presence absence OTU table
macro18S_taxa_sums_17 <- as.data.frame(filter_taxa(macro18S_group17, 
                                                  function(x) sum(x))) 
### select OTUs present in at least 2 samples ###
#select OTUs with sample count over your threshold
keep_macro18S_16 <- row.names(
  macro18S_taxa_sums_16)[which(macro18S_taxa_sums_16[,1] >= 2)] 
#select OTUs with sample count over your threshold
keep_macro18S_17 <- row.names(
  macro18S_taxa_sums_17)[which(macro18S_taxa_sums_17[,1] >= 2)] 

### Venn Diagram ###
macro18S_years <- list("2016"=keep_macro18S_16,
                      "2017"=keep_macro18S_17)
venn_macro18S_years <- ggvenn(
  macro18S_years,
  fill_color = c("#e41a1c","#377eb8"),
  fill_alpha = 0.6,
  stroke_color = "black",
  stroke_alpha = 0.7,
  stroke_size = 1,
  stroke_linetype = "solid",
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 5
)
venn_macro18S_years <- venn_macro18S_years +
  annotate("text", x = 0, y = 1.6, label = "Macro18S_16_17", size = 8)

venn_macro18S_years

ggsave("R_Code_and_Analysis/varpart/macro18S_16_17_taxa_plot.tiff", plot = venn_macro18S_years, width=250, height=200, units="mm",dpi=300, compression = "lzw", type = "cairo")
