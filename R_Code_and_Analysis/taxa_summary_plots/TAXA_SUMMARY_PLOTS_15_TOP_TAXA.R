### TAXA SUMMARY PLOTS ###
### Author: Bianca Trevizan Segovia ###
### Date created: December 08th, 2020 ###
### *uses the nonrarefied phyloseq object

# ============================================================
# Tutorial on drawing a taxa plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
library(phyloseq)
library(tidyverse)
library(ggplot2)

#########################################
############ 16S prokaryotes ############
#########################################
# 
# ###load phyloseq object 
# prokary_data_notrarefied <- readRDS("Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")
# 
# ### asv table and metadata
# prokary_abund_table <- otu_table(prokary_data_notrarefied)
# prokary_meta_table<-sample_data(prokary_data_notrarefied)
# 
# #Apply proportion normalisation (relative abundance)
# x<-prokary_abund_table/rowSums(prokary_abund_table)
# 
# #Sort by most abundant to less abundant
# x<-x[,order(colSums(x),decreasing=TRUE)]
# rowSums(x)
# 
# #Extract list of top N Taxa
# N<-15
# taxa_list<-colnames(x)[1:N]
# #remove "__Unknown__" and add it to others
# taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
# N<-length(taxa_list)
# 
# #Generate a new table with everything added to Others
# x <- t(x)
# colSums(x)
# new_x<-data.frame(x[1:N,]) #get most abundant taxa
# new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
# new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset
# 
# #fix column labels
# #pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
# alltaxa <- as.data.frame(unclass(tax_table(prokary_data_notrarefied)))
# taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
# colnames(new_x) 
# # colnames(new_x) = "ASV25"  "ASV28"  "ASV6"   "ASV54"  "ASV24"  "ASV14"  "ASV27"  "ASV22"  "ASV111" "ASV4" "ASV167" "ASV61"  "ASV75"  "ASV18"  "ASV37"  "Others"
# ###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
# 
# mylabels <- c("ASV25_Gammaproteobacteria_uncultured",
#               "ASV28_Gammaproteobacteria_Granulosicoccus",
#               "ASV6_Gammaproteobacteria_Methylophilus_sp._EM8",
#               "ASV54_Acidimicrobiia_uncultured",
#               "ASV24_Bacteroidia_Saprospiraceae",
#               "ASV14_Gammaproteobacteria_Methylotenera",
#               "ASV27_Acidimicrobiia_uncultured",
#               "ASV22_Alphaproteobacteria_Sulfitobacter",
#               "ASV111_Alphaproteobacteria_Sulfitobacter2",
#               "ASV4_Gammaproteobacteria_Methylophagaceae", 
#               "ASV167_Alphaproteobacteria_Rickettsia",
#               "ASV61_Alphaproteobacteria_Sulfitobacter3",
#               "ASV75_Alphaproteobacteria_Donghicola",
#               "ASV18_Gammaproteobacteria_Cocleimonas", 
#               "ASV37_Gammaproteobacteria_Granulosicoccus2" ,
#               "Others")
# 
# colour_taxa <- c("yellow", 
#                  "tan1",
#                  "darkred",
#                  "plum",
#                  "blue",
#                  "firebrick", 
#                  "plum2",
#                  "olivedrab",
#                  "olivedrab3", 
#                  "firebrick4",
#                  "aquamarine",  
#                  "olivedrab2", 
#                  "lightgoldenrod1" , 
#                  "darkslateblue",
#                  "tan3", 
#                  "grey36")
# 
# #Add metadata for plotting to new_x
# new_x$region <- sample_data(prokary_data_notrarefied)$region[which(row.names(sample_data(prokary_data_notrarefied)) %in% row.names(new_x))]
# new_x$year <- sample_data(prokary_data_notrarefied)$year[which(row.names(sample_data(prokary_data_notrarefied)) %in% row.names(new_x))]
# new_x$SampleID <- row.names(new_x)
# library(reshape)
# new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)
# 
# new_x_melt <- new_x_melt %>% 
#   arrange(region, year)
# 
# # re-order the factor levels before the plot
# new_x_melt <- new_x_melt %>% 
#   dplyr::mutate(region=recode(region,
#                               "mcmullin" = "mcmullins"))
# 
# new_x_melt <- new_x_melt %>% 
#   dplyr::mutate(region_year = paste(region, year, sep = "_"))
# 
# new_x_melt$region <- factor(new_x_melt$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# 
# new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("choked_2015", "choked_2016","choked_2017","choked_2018","pruth_2016", "pruth_2017", "pruth_2018", "triquet_2015", "triquet_2016","triquet_2017","triquet_2018","goose_2015","goose_2016","mcmullins_2015", "mcmullins_2016"))
# 
# ### Taxa plot
# p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
#   geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
#   facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
# p <- p + ggtitle("Prokaryotes")
# p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)
# 
# p <- p + theme_bw() + ylab("Relative abundance")
# p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
# p <- p + theme (axis.title.x=element_blank(),
#                 #axis.text.x=element_blank(),
#                 axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
#                 axis.ticks.x=element_blank(),
#                 axis.text.y=element_text(size = 14), #change font size of numbers
#                 axis.title.y=element_text(size = 18), #change font size of y title
#                 panel.grid.major.x=element_blank(),
#                 panel.grid.minor.x=element_blank(),
#                 legend.text=element_text(size=11),
#                 legend.text.align=0,
#                 plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
# p
# 
# ggsave("R_Code_and_Analysis/taxa_summary_plots/Prokaryotes_taxa_plot.tiff", plot = p, width=400, height=150, units="mm",dpi=300, compression = "lzw", type = "cairo")
# 
# #########################################
# ############ 18S microeukaryotes ########
# #########################################
# 
# ###load phyloseq object 
# microeuk_data_notrarefied <- readRDS("Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds")
# 
# ### asv table and metadata
# microeuk_abund_table <- otu_table(microeuk_data_notrarefied)
# microeuk_meta_table<-sample_data(microeuk_data_notrarefied)
# 
# #Apply proportion normalisation (relative abundance)
# x<-microeuk_abund_table/rowSums(microeuk_abund_table)
# 
# #Sort by most abundant to less abundant
# x<-x[,order(colSums(x),decreasing=TRUE)]
# rowSums(x)
# 
# #Extract list of top N Taxa
# N<-15
# taxa_list<-colnames(x)[1:N]
# #remove "__Unknown__" and add it to others
# taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
# N<-length(taxa_list)
# 
# #Generate a new table with everything added to Others
# x <- t(x)
# colSums(x)
# new_x<-data.frame(x[1:N,]) #get most abundant taxa
# new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
# new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset
# 
# 
# #fix column labels
# #pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
# alltaxa <- as.data.frame(unclass(tax_table(microeuk_data_notrarefied)))
# taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
# colnames(new_x) 
# # colnames(new_x) = "ASV2"   "ASV3"   "ASV16"  "ASV7"   "ASV6"   "ASV10"  "ASV18"  "ASV21"   "ASV19" "ASV31" "ASV77"  "ASV80"  "ASV35"  "ASV46"  "ASV39" "Others"
# ###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
# taxa_N
# 
# mylabels <- c("ASV2_Dinoflagellata_Prorocentrum_foraminosum",
#               "ASV3_Ochrophyta_uncultured1",
#               "ASV16_Ochrophyta_Bacillariophyceae1",
#               "ASV7_Dinoflagellata_uncultured",
#               "ASV6_Ciliophora_Ephelota_sp._QD-5",
#               "ASV10_Ochrophyta_Bacillariophyceae2",
#               "ASV18_Dinoflagellata_Amphidinium_steinii",
#               "ASV21_Ochrophyta_Bacillariophyceae3",
#               "ASV19_Protalveolata_Syndiniales", 
#               "ASV31_Ochrophyta_Bacillariophyceae4",
#               "ASV77_Ochrophyta_uncultured2",
#               "ASV80_Peronosporomycetes_Olpidiopsis_sp._4439",
#               "ASV35_Protalveolata_Perkinsidae(parasites_of_molluscs)", 
#               "ASV46_Coscinodiscophytina_Cylindrotheca_closterium" ,
#               "ASV39_Alveolata",
#               "Others")
# 
# colour_taxa <- c("turquoise3", 
#                  "slateblue4",
#                  "plum2",
#                  "turquoise",
#                  "yellow",
#                  "purple", 
#                  "turquoise4",
#                  "purple1",
#                  "plum4", 
#                  "firebrick4",
#                  "aquamarine",  
#                  "olivedrab2", 
#                  "lightgoldenrod1" , 
#                  "darkslateblue",
#                  "tan3", 
#                  "grey36")
# 
# #Add metadata for plotting to new_x
# new_x$region <- sample_data(microeuk_data_notrarefied)$region[which(row.names(sample_data(microeuk_data_notrarefied)) %in% row.names(new_x))]
# new_x$year <- sample_data(microeuk_data_notrarefied)$year[which(row.names(sample_data(microeuk_data_notrarefied)) %in% row.names(new_x))]
# new_x$SampleID <- row.names(new_x)
# library(reshape)
# new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)
# 
# # re-order the factor levels before the plot
# new_x_melt <- new_x_melt %>% 
#   dplyr::mutate(region=recode(region,
#                               "mcmullin" = "mcmullins"))
# 
# new_x_melt <- new_x_melt %>% 
#   dplyr::mutate(region_year = paste(region, year, sep = "_"))
# 
# new_x_melt$region <- factor(new_x_melt$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# 
# new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("choked_2015", "choked_2016","choked_2017","choked_2018","pruth_2016", "pruth_2017", "pruth_2018", "triquet_2015", "triquet_2016","triquet_2017","triquet_2018","goose_2015","goose_2016","mcmullins_2015"))
# 
# ### Taxa plot
# p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
#   geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
#   facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
# p <- p + ggtitle("Microeukaryotes")
# p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)
# 
# p <- p + theme_bw() + ylab("Relative abundance")
# p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
# p <- p + theme (axis.title.x=element_blank(),
#                 #axis.text.x=element_blank(),
#                 axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
#                 axis.ticks.x=element_blank(),
#                 axis.text.y=element_text(size = 14), #change font size of numbers
#                 axis.title.y=element_text(size = 18), #change font size of y title
#                 panel.grid.major.x=element_blank(),
#                 panel.grid.minor.x=element_blank(),
#                 legend.text=element_text(size=10),
#                 legend.text.align=0,
#                 plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
# p
# 
# ggsave("R_Code_and_Analysis/taxa_summary_plots/Microeukaryotes_taxa_plot.tiff", plot = p, width=420, height=150, units="mm",dpi=300, compression = "lzw", type = "cairo")
# 
#################################
############ Macroeuk18S ########
#################################

###load phyloseq object
macro18S_data_notrarefied <- readRDS("Data/micro_eukaryotes/macro18S_filtered.rds")

### asv table and metadata
macro18S_abund_table <- otu_table(macro18S_data_notrarefied)
macro18S_meta_table<-sample_data(macro18S_data_notrarefied)

#Apply proportion normalisation (relative abundance)
x<-macro18S_abund_table/rowSums(macro18S_abund_table)

#Sort by most abundant to less abundant
x<-x[,order(colSums(x),decreasing=TRUE)]
rowSums(x)

#Extract list of top N Taxa
N<-15
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
N<-length(taxa_list)

#Generate a new table with everything added to Others
x <- t(x)
colSums(x)
new_x<-data.frame(x[1:N,]) #get most abundant taxa
new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset


#fix column labels
#pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
alltaxa <- as.data.frame(unclass(tax_table(macro18S_data_notrarefied)))
taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
colnames(new_x)
# colnames(new_x) =  [1] "ASV20"  "ASV9"   "ASV34"  "ASV8"   "ASV43"  "ASV17" "ASV93"  "ASV70"  "ASV24"  "ASV171" "ASV85"  "ASV59" "ASV117" "ASV51"  "ASV97"  "Others"
###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
taxa_N

mylabels <- c("ASV20_Chromadorea(nematoda)",
              "ASV9_Hydroidolina_Clytia_gracilis",
              "ASV34_Monogononta_Ploimida(rotifer)",
              "ASV8_Copepoda_Harpacticus_sp._France_RJH_2007",
              "ASV43_Hydroidolina_Clytia_gracilis2",
              "ASV17_Podocopa_Podocopida(ostracod)",
              "ASV93_Monogononta_Ploimida(rotifer)_Proales_reinhardti",
              "ASV70_Copepoda_Porcellidium_ofunatense",
              "ASV24_Podocopa_Podocopida(ostracod)2",
              "ASV171_Terebellida_Nicolea_zostericola",
              "ASV85_Rhabditophora_Microstomum_papillosum",
              "ASV59_Heteroconchia_Fulvia_mutica",
              "ASV117_Chromadorea(nematoda)2",
              "ASV51_Podocopa_Aurila_disparata(ostracod)",
              "ASV97_Copepoda_Amonardia_coreana" ,
              "Others")

colour_taxa <- c("turquoise3",
                 "slateblue4",
                 "springgreen3",
                 "lightyellow",
                 "yellow",
                 "purple",
                 "seagreen3",
                 "navajowhite",
                 "plum4",
                 "firebrick4",
                 "lightgoldenrod1",
                 "olivedrab2",
                 "aquamarine" ,
                 "darkslateblue",
                 "papayawhip",
                 "grey36")

#Add metadata for plotting to new_x
new_x$region <- sample_data(macro18S_data_notrarefied)$region[which(row.names(sample_data(macro18S_data_notrarefied)) %in% row.names(new_x))]
new_x$year <- sample_data(macro18S_data_notrarefied)$year[which(row.names(sample_data(macro18S_data_notrarefied)) %in% row.names(new_x))]
new_x$SampleID <- row.names(new_x)
library(reshape)
new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)

# re-order the factor levels before the plot
new_x_melt <- new_x_melt %>%
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

new_x_melt <- new_x_melt %>%
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

new_x_melt$region <- factor(new_x_melt$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))

new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("choked_2015", "choked_2016","choked_2017","choked_2018","pruth_2016", "pruth_2017", "pruth_2018", "triquet_2015", "triquet_2016","triquet_2017","triquet_2018","goose_2015","goose_2016","mcmullins_2015"))

### Taxa plot
p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
  geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
  facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
p <- p + ggtitle("Macro18S")
p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)

p <- p + theme_bw() + ylab("Relative abundance")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
p <- p + theme (axis.title.x=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size = 14), #change font size of numbers
                axis.title.y=element_text(size = 18), #change font size of y title
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                legend.text=element_text(size=10),
                legend.text.align=0,
                plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
p

ggsave("R_Code_and_Analysis/taxa_summary_plots/Macro18S_taxa_plot.tiff", plot = p, width=420, height=150, units="mm",dpi=300, compression = "lzw", type = "cairo")

ggsave("R_Code_and_Analysis/taxa_summary_plots/Macro18S_taxa_plot.png", plot = p, width=420, height=150, units="mm",dpi=300)

#############################
############ Inverts ########
#############################

###load phyloseq object 
inverts_data_notrarefied <- readRDS("Data/macro_eukaryotes/invert_phyloseq.rds")

### taxa table and metadata
#Transpose the data to have sample names on rows
inverts_abund_table <- otu_table(inverts_data_notrarefied)
inverts_abund_table<-t(inverts_abund_table)
inverts_meta_table <- sample_data(inverts_data_notrarefied)

#Apply proportion normalisation (relative abundance)
x<-inverts_abund_table/rowSums(inverts_abund_table)

#Sort by most abundant to less abundant
x<-x[,order(colSums(x),decreasing=TRUE)]
rowSums(x)

#Extract list of top N Taxa
N<-15
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
N<-length(taxa_list)

#Generate a new table with everything added to Others
x <- t(x)
colSums(x)
new_x<-data.frame(x[1:N,]) #get most abundant taxa
new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset


#fix column labels
#pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
alltaxa <- as.data.frame(unclass(tax_table(inverts_data_notrarefied)))
taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
colnames(new_x) 
# colnames(new_x) =  "Caprella.laeviuscula"       "Caprella.californica"       "Gammaridean" "Lacuna.vincta"              "Porcellidium"               "Leptochelia"                "Harpacticoida"              "Aoroides"                   "Platynereis.bicanaliculata" "Lacuna"                     "Pontogeneia"                "Leptochelia.dubia"          "Pentidotea"                 "Caprella"                   "Photis"                     "Others"  
###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
taxa_N

mylabels <- c("Caprella.laeviuscula",
              "Caprella.californica",
              "Gammaridean",
              "Lacuna.vincta",
              "Porcellidium",
              "Leptochelia",
              "Harpacticoida",
              "Aoroides",
              "Platynereis.bicanaliculata",
              "Lacuna",
              "Pontogeneia",
              "Leptochelia.dubia",
              "Pentidotea",
              "Caprella",
              "Photis",
              "Others")

colour_taxa <- c("tan1", 
                 "tan2",
                 "springgreen3",
                 "purple", 
                 "darkred",
                 "yellow",
                 "red3",
                 "navajowhite",
                 "plum4", 
                 "purple4",
                 "olivedrab2",
                 "lightgoldenrod1",  
                 "aquamarine" , 
                 "sandybrown",
                 "papayawhip", 
                 "grey36")

#Add metadata for plotting to new_x
new_x$region <- sample_data(inverts_data_notrarefied)$region[which(row.names(sample_data(inverts_data_notrarefied)) %in% row.names(new_x))]
new_x$year <- sample_data(inverts_data_notrarefied)$year[which(row.names(sample_data(inverts_data_notrarefied)) %in% row.names(new_x))]
new_x$SampleID <- row.names(new_x)
library(reshape)
new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)

# re-order the factor levels before the plot
new_x_melt <- new_x_melt %>% 
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

new_x_melt <- new_x_melt %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

new_x_melt$region <- factor(new_x_melt$region, levels=c("pruth", "choked", "triquet","goose","mcmullins"))
new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("pruth_2016", "pruth_2017", "choked_2014", "choked_2015","choked_2016","choked_2017","triquet_2014", "triquet_2015","triquet_2016","triquet_2017","goose_2014","goose_2015","goose_2016","mcmullins_2014", "mcmullins_2015"))

### Taxa plot
p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
  geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
  facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
p <- p + ggtitle("inverts")
p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)

p <- p + theme_bw() + ylab("Relative abundance")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
p <- p + theme (axis.title.x=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size = 14), #change font size of numbers
                axis.title.y=element_text(size = 18), #change font size of y title
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                legend.text=element_text(size=10),
                legend.text.align=0,
                plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
p

ggsave("R_Code_and_Analysis/taxa_summary_plots/inverts_taxa_plot.tiff", plot = p, width=420, height=150, units="mm",dpi=300, compression = "lzw", type = "cairo")
ggsave("R_Code_and_Analysis/taxa_summary_plots/inverts_taxa_plot.png", plot = p, width=420, height=150, units="mm",dpi=300)

#############################
############ Inverts ########
#############################
year_16_17 <- c("2016", "2017")
new_x_melt_16_17 <- new_x_melt %>% 
  filter(year %in% year_16_17)


### Taxa plot
p <- ggplot(new_x_melt_16_17,aes(SampleID,value,fill=Taxa))+
  geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
  facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
p <- p + ggtitle("inverts")
p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)

p <- p + theme_bw() + ylab("Relative abundance")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
p <- p + theme (axis.title.x=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size = 14), #change font size of numbers
                axis.title.y=element_text(size = 18), #change font size of y title
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                legend.text=element_text(size=10),
                legend.text.align=0,
                plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
p

p <- p + theme (legend.position = "bottom")

ggsave("R_Code_and_Analysis/taxa_summary_plots/inverts_taxa_plot_16_17.tiff", plot = p, width=420, height=150, units="mm",dpi=300, compression = "lzw", type = "cairo")
ggsave("R_Code_and_Analysis/taxa_summary_plots/inverts_taxa_plot_16_17.png", plot = p, width=420, height=150, units="mm",dpi=300)

######################################################
###################   EXTRA    #######################
######################################################
 ### 16S with 30 most abundant taxa ###

# #########################################
# ############ 16S prokaryotes ############
# #########################################
# 
# ###load phyloseq object 
# prokary_data_notrarefied <- readRDS("Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")
# 
# ### asv table and metadata
# prokary_abund_table <- otu_table(prokary_data_notrarefied)
# prokary_meta_table<-sample_data(prokary_data_notrarefied)
# 
# #Apply proportion normalisation (relative abundance)
# x<-prokary_abund_table/rowSums(prokary_abund_table)
# 
# #Sort by most abundant to less abundant
# x<-x[,order(colSums(x),decreasing=TRUE)]
# rowSums(x)
# 
# #Extract list of top N Taxa
# N<-30
# taxa_list<-colnames(x)[1:N]
# #remove "__Unknown__" and add it to others
# taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
# N<-length(taxa_list)
# 
# #Generate a new table with everything added to Others
# x <- t(x)
# colSums(x)
# new_x<-data.frame(x[1:N,]) #get most abundant taxa
# new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
# new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset
# 
# #fix column labels
# #pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
# alltaxa <- as.data.frame(unclass(tax_table(prokary_data_notrarefied)))
# taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
# colnames(new_x) 
# # colnames(new_x) = "ASV25"  "ASV28"  "ASV6"   "ASV54"  "ASV24"  "ASV14"  "ASV27"  "ASV22"  "ASV111" "ASV4" "ASV167" "ASV61"  "ASV75"  "ASV18"  "ASV37"  "ASV68"  "ASV21"  "ASV58"  "ASV43"  "ASV82" "ASV259" "ASV19"  "ASV69"  "ASV64"  "ASV298" "ASV59"  "ASV94"  "ASV78"  "ASV85"  "ASV102" "Others"
# ###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
# taxa_N
# mylabels <- c("ASV25_Gammaproteobacteria_uncultured",
#               "ASV28_Gammaproteobacteria_Granulosicoccus",
#               "ASV6_Gammaproteobacteria_Methylophilus_sp._EM8",
#               "ASV54_Acidimicrobiia_uncultured",
#               "ASV24_Saprospiraceae_uncultured",
#               "ASV14_Gammaproteobacteria_Methylotenera",
#               "ASV27_Acidimicrobiia_uncultured",
#               "ASV22_Alphaproteobacteria_Sulfitobacter",
#               "ASV111_Alphaproteobacteria_Sulfitobacter2",
#               "ASV4_Gammaproteobacteria_Methylophagaceae", 
#               "ASV167_Alphaproteobacteria_Rickettsia",
#               "ASV61_Alphaproteobacteria_Sulfitobacter3",
#               "ASV75_Alphaproteobacteria_Donghicola",
#               "ASV18_Gammaproteobacteria_Cocleimonas", 
#               "ASV37_Gammaproteobacteria_Granulosicoccus2" ,
#               "ASV68_Gammaproteobacteria_Granulosicoccus3",
#               "ASV21_Flavobacteriaceae_Maribacter_sp._MGE_SAT_274",
#               "ASV58_Bacteroidia_uncultured",
#               "ASV43_Hyphomonadaceae_Hellea",
#               "ASV82_Rhodobacteraceae_uncultured",
#               "ASV259_Rhodobacteraceae_Litoreibacter_halocynthiae",
#               "ASV19_Flavobacteriaceae_Polaribacter_sp._KJF12-1",
#               "ASV69_Saprospiraceae_uncultured2",
#               "ASV64_Alphaproteobacteria_Sulfitobacter4",
#               "ASV298_Alphaproteobacteria_Loktanella_vestfoldensis",
#               "ASV59_Saprospiraceae_uncultured3" ,
#               "ASV94_Gammaproteobacteria_Granulosicoccus4" ,
#               "ASV78_Saprospiraceae_uncultured4" ,
#               "ASV85_Saprospiraceae_uncultured5",
#               "ASV102_Flavobacteriaceae_Winogradskyella",
#               "Others")
# 
# colour_taxa <- c("yellow", 
#                  "tan1",
#                  "darkred",
#                  "plum",
#                  "blue",
#                  "firebrick", 
#                  "plum2",
#                  "olivedrab",
#                  "olivedrab3", 
#                  "firebrick4",
#                  "aquamarine",  
#                  "olivedrab2", 
#                  "lightgoldenrod1" , 
#                  "darkslateblue",
#                  "tan3",
#                  "tan4", 
#                  "tan2",
#                  "springgreen3",
#                  "purple", 
#                  "firebrick2",
#                  "gold1",
#                  "red3",
#                  "navajowhite",
#                  "plum4", 
#                  "purple4",
#                  "darkseagreen4",
#                  "lightgoldenrod1",  
#                  "darkslategray4" , 
#                  "sandybrown",
#                  "papayawhip",
#                  "grey36")
# 
# #Add metadata for plotting to new_x
# new_x$region <- sample_data(prokary_data_notrarefied)$region[which(row.names(sample_data(prokary_data_notrarefied)) %in% row.names(new_x))]
# new_x$year <- sample_data(prokary_data_notrarefied)$year[which(row.names(sample_data(prokary_data_notrarefied)) %in% row.names(new_x))]
# new_x$SampleID <- row.names(new_x)
# library(reshape)
# new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)
# 
# new_x_melt <- new_x_melt %>% 
#   arrange(region, year)
# 
# # re-order the factor levels before the plot
# new_x_melt <- new_x_melt %>% 
#   dplyr::mutate(region=recode(region,
#                               "mcmullin" = "mcmullins"))
# 
# new_x_melt <- new_x_melt %>% 
#   dplyr::mutate(region_year = paste(region, year, sep = "_"))
# 
# new_x_melt$region <- factor(new_x_melt$region, levels=c("choked", "pruth", "triquet","goose","mcmullins"))
# 
# new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("choked_2015", "choked_2016","choked_2017","choked_2018","pruth_2016", "pruth_2017", "pruth_2018", "triquet_2015", "triquet_2016","triquet_2017","triquet_2018","goose_2015","goose_2016","mcmullins_2015", "mcmullins_2016"))
# 
# ### Taxa plot
# p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
#   geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
#   facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
# p <- p + ggtitle("Prokaryotes")
# p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)
# 
# p <- p + theme_bw() + ylab("Relative abundance")
# p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
# p <- p + theme (axis.title.x=element_blank(),
#                 #axis.text.x=element_blank(),
#                 axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
#                 axis.ticks.x=element_blank(),
#                 axis.text.y=element_text(size = 14), #change font size of numbers
#                 axis.title.y=element_text(size = 18), #change font size of y title
#                 panel.grid.major.x=element_blank(),
#                 panel.grid.minor.x=element_blank(),
#                 legend.text=element_text(size=11),
#                 legend.text.align=0,
#                 plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
# 
# p <- p + theme (legend.position = "bottom")
# 
# p
# 
# ggsave("R_Code_and_Analysis/taxa_summary_plots/Prokaryotes_taxa_plot.tiff", plot = p, width=510, height=250, units="mm",dpi=300, compression = "lzw", type = "cairo")


######################################
###########   EXTRA ##################
######################################

#####################################################
############ 16S prokaryotes GENUS LEVEL ############
#####################################################

###load phyloseq object 
prokary_data_notrarefied <- readRDS("Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds")
#View(as.data.frame(tax_table(prokary_data_notrarefied)))

prokary_genus_notrarefied <- prokary_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

### asv table and metadata
prokary_genus_abund_table <- otu_table(prokary_genus_notrarefied)
prokary_genus_meta_table<-sample_data(prokary_genus_notrarefied)

#Apply proportion normalisation (relative abundance)
x<-prokary_genus_abund_table/rowSums(prokary_genus_abund_table)

#Sort by most abundant to less abundant
x<-x[,order(colSums(x),decreasing=TRUE)]
rowSums(x)

#Extract list of top N Taxa
N<-30
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
N<-length(taxa_list)

#Generate a new table with everything added to Others
x <- t(x)
colSums(x)
new_x<-data.frame(x[1:N,]) #get most abundant taxa
new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset

#fix column labels
#pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
alltaxa <- as.data.frame(unclass(tax_table(prokary_genus_notrarefied)))
taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
colnames(new_x) 
# colnames(new_x) = "ASV24"  "ASV28"  "ASV22"  "ASV54"  "ASV14"  "ASV6"   "ASV167" "ASV21"  "ASV211" "ASV91" "ASV18"  "ASV83"  "ASV227" "ASV4"   "ASV312" "ASV102" "ASV104" "ASV75"  "ASV19"  "ASV531" "ASV177" "ASV43"  "ASV282" "ASV33"  "ASV259" "ASV219" "ASV326" "ASV524" "ASV246" "ASV602" "Others"
###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
taxa_N

mylabels <- c("ASV24_Saprospiraceae_uncultured",
              "ASV28_Thiohalorhabdaceae_Granulosicoccus",
              "ASV22_Rhodobacteraceae_Sulfitobacter",
              "ASV54_Microtrichaceae_uncultured",
              "ASV14_Methylophilaceae_Methylotenera",
              "ASV6_Methylophilaceae_Methylophilus",
              "ASV167_Rickettsiaceae_Rickettsia",
              "ASV21_Flavobacteriaceae_Maribacter",
              "ASV211_Verrucomicrobiales_DEV007_uncultured",
              "ASV91_Hyphomonadaceae_Litorimonas", 
              "ASV18_Thiotrichaceae_Cocleimonas",
              "ASV83_Rhodobacteraceae_Loktanella",
              "ASV227_Rubritaleaceae_Roseibacillus",
              "ASV4_Methylophagaceae_uncultured", 
              "ASV312_Pirellulaceae_Blastopirellula" ,
              "ASV102_Flavobacteriaceae_Winogradskyella",
              "ASV104_Sphingomonadaceae_Altererythrobacter",
              "ASV75_Rhodobacteraceae_Donghicola",
              "ASV19_Flavobacteriaceae_Polaribacter",
              "ASV531_Phormidesmiaceae_Phormidesmis_ANT.LACV5.1",
              "ASV177_Rhodobacteraceae_Pacificibacter",
              "ASV43_Hyphomonadaceae_Hellea",
              "ASV282_Saprospiraceae_Lewinella",
              "ASV33_Thiotrichaceae_Leucothrix",
              "ASV259_Rhodobacteraceae_Litoreibacter",
              "ASV219_Flavobacteriaceae_Kordia",
              "ASV326_Rhizobiaceae_Pseudahrensia",
              "ASV524_Cellvibrionaceae_uncultured",
              "ASV246_Arenicellacea_Arenicella",
              "ASV602_ Rubritaleaceae_Persicirhabdus",
              "Others")

colour_taxa <- c("yellow", 
                   "tan1",
                   "darkred",
                   "plum",
                   "blue",
                   "firebrick", 
                   "plum2",
                   "olivedrab",
                   "olivedrab3", 
                   "firebrick4",
                   "aquamarine",  
                   "olivedrab2", 
                   "lightgoldenrod1" , 
                   "darkslateblue",
                   "tan3",
                   "tan4", 
                   "tan2",
                   "springgreen3",
                   "purple", 
                   "firebrick2",
                   "gold1",
                   "red3",
                   "navajowhite",
                   "plum4", 
                   "purple4",
                   "darkseagreen4",
                   "lightgoldenrod1",  
                   "darkslategray4" , 
                   "sandybrown",
                   "papayawhip",
                   "grey36")

#Add metadata for plotting to new_x
new_x$region <- sample_data(prokary_genus_notrarefied)$region[which(row.names(sample_data(prokary_genus_notrarefied)) %in% row.names(new_x))]
new_x$year <- sample_data(prokary_genus_notrarefied)$year[which(row.names(sample_data(prokary_genus_notrarefied)) %in% row.names(new_x))]
new_x$SampleID <- row.names(new_x)
library(reshape)
new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)

new_x_melt <- new_x_melt %>% 
  arrange(region, year)

# re-order the factor levels before the plot
new_x_melt <- new_x_melt %>% 
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

new_x_melt <- new_x_melt %>% 
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

new_x_melt$region <- factor(new_x_melt$region, levels=c("pruth", "choked", "triquet","goose","mcmullins"))

new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("pruth_2016", "pruth_2017", "pruth_2018","choked_2015", "choked_2016","choked_2017","choked_2018", "triquet_2015", "triquet_2016","triquet_2017","triquet_2018","goose_2015","goose_2016","mcmullins_2015", "mcmullins_2016"))

### Taxa plot
p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
  geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
  facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
p <- p + ggtitle("Prokaryotes_GENUS")
p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)

p <- p + theme_bw() + ylab("Relative abundance")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
p <- p + theme (axis.title.x=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size = 14), #change font size of numbers
                axis.title.y=element_text(size = 18), #change font size of y title
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                legend.text=element_text(size=11),
                legend.text.align=0,
                plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
p

p <- p + theme (legend.position = "bottom")

p

ggsave("R_Code_and_Analysis/taxa_summary_plots/Prokaryotes_GENUS_taxa_plot.tiff", plot = p, width=510, height=250, units="mm",dpi=300, compression = "lzw", type = "cairo")

ggsave("R_Code_and_Analysis/taxa_summary_plots/Prokaryotes_GENUS_taxa_plot.png", plot = p, width=510, height=250, units="mm",dpi=300)

#####################################################
############ 18S microeukaryotes GENUS LEVEL ########
#####################################################

###load phyloseq object
microeuk_data_notrarefied <- readRDS("Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds")

microeuk_data_notrarefied <- microeuk_data_notrarefied %>% 
  tax_glom(taxrank = "Rank6") # agglomerate at "genus" level
#View(as.data.frame(tax_table(prokary_genus_notrarefied)))

### asv table and metadata
microeuk_abund_table <- otu_table(microeuk_data_notrarefied)
microeuk_meta_table<-sample_data(microeuk_data_notrarefied)

#Apply proportion normalisation (relative abundance)
x<-microeuk_abund_table/rowSums(microeuk_abund_table)

#Sort by most abundant to less abundant
x<-x[,order(colSums(x),decreasing=TRUE)]
rowSums(x)

#Extract list of top N Taxa
N<-15
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
N<-length(taxa_list)

#Generate a new table with everything added to Others
x <- t(x)
colSums(x)
new_x<-data.frame(x[1:N,]) #get most abundant taxa
new_x <- as.data.frame(t(new_x)) #transpose and cast as data frame
new_x$Others <- colSums(x[(N+1):length(row.names(x)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset


#fix column labels
#pull taxa table subset by 15 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
alltaxa <- as.data.frame(unclass(tax_table(microeuk_data_notrarefied)))
taxa_N <- alltaxa[which(row.names(alltaxa) %in% taxa_list),]
colnames(new_x)
# colnames(new_x) = "ASV2"   "ASV18"  "ASV6"   "ASV34"  "ASV80"  "ASV46"  "ASV73"  "ASV88" "ASV105" "ASV147" "ASV111" "ASV229" "ASV133" "ASV172" "ASV320" "Others"
###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
taxa_N

mylabels <- c("ASV2_Dinoflagellata_Prorocentrum",
              "ASV18_Dinoflagellata_Amphidinium",
              "ASV6_Ciliophora_Ephelota",
              "ASV34_Rotifer_Ploimida",
              "ASV80_Peronosporomycetes_Olpidiopsis(pathogens algae)",
              "ASV46_Diatom_Fragilariales",
              "ASV73_Cercozoa_Vampyrellidae",
              "ASV88_Peronosporomycetes_uncultured",
              "ASV105_Diatom_Navicula",
              "ASV147_Cercozoa_Novel_Clade_Gran-1",
              "ASV111_Ciliophora_Eufolliculina",
              "ASV229_Labyrinthulomycetes_Aplanochytrium",
              "ASV133_Diatom_Cymbella",
              "ASV172_Labyrinthulomycetes_Labyrinthula" ,
              "ASV320_Fungi_Lobulomycetaceae",
              "Others")

colour_taxa <- c("turquoise4",
                 "turquoise",
                 "chocolate4",
                 "plum2",
                 "yellow",
                 "purple",
                 "red",
                 "plum4",
                 "purple1",
                 "red3",
                 "salmon3",
                 "olivedrab2",
                 "lightgoldenrod1" ,
                 "darkslateblue",
                 "tan3",
                 "grey36")

#Add metadata for plotting to new_x
new_x$region <- sample_data(microeuk_data_notrarefied)$region[which(row.names(sample_data(microeuk_data_notrarefied)) %in% row.names(new_x))]
new_x$year <- sample_data(microeuk_data_notrarefied)$year[which(row.names(sample_data(microeuk_data_notrarefied)) %in% row.names(new_x))]
new_x$SampleID <- row.names(new_x)
library(reshape)
new_x_melt <- melt.data.frame(new_x, id.vars = c("SampleID","region", "year"), variable_name = "Taxa", na.rm=TRUE)

# re-order the factor levels before the plot
new_x_melt <- new_x_melt %>%
  dplyr::mutate(region=recode(region,
                              "mcmullin" = "mcmullins"))

new_x_melt <- new_x_melt %>%
  dplyr::mutate(region_year = paste(region, year, sep = "_"))

new_x_melt$region <- factor(new_x_melt$region, levels=c("pruth", "choked", "triquet","goose","mcmullins"))

new_x_melt$region_year <- factor(new_x_melt$region_year, levels=c("pruth_2016", "pruth_2017", "pruth_2018", "choked_2015", "choked_2016","choked_2017","choked_2018","triquet_2015", "triquet_2016","triquet_2017","triquet_2018","goose_2015","goose_2016","mcmullins_2015"))

### Taxa plot
p <- ggplot(new_x_melt,aes(SampleID,value,fill=Taxa))+
  geom_bar(stat="identity", width=1, color="black", position = "stack") + #height of the column equal the value, stacked
  facet_grid(. ~ region + year, drop=TRUE,scale="free",space="free_x")
p <- p + ggtitle("Microeukaryotes")
p <- p + scale_fill_manual(values=colour_taxa, labels = mylabels)

p <- p + theme_bw() + ylab("Relative abundance")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="snow2"))+theme(panel.spacing = unit(0.2, "lines"))
p <- p + theme (axis.title.x=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.x=element_text(angle = 90, hjust = 1, size = 5),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size = 14), #change font size of numbers
                axis.title.y=element_text(size = 18), #change font size of y title
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                legend.text=element_text(size=10),
                legend.text.align=0,
                plot.title = element_text(hjust = 0.5, size = 18, face="bold")) #center plot title and set font size
p

p <- p + theme (legend.position = "bottom")

ggsave("R_Code_and_Analysis/taxa_summary_plots/Microeukaryotes_GENUS_taxa_plot.png", plot = p, width=420, height=150, units="mm",dpi=300)

