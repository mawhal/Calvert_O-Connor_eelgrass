# started oct 13, 2024 by D. Loughnan

# starting to explore the community data:
# 1. how many sites and what sites do we have community data from
# 2. Calculations of the alpha and beta and omega diversity 

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) 

#setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/ranges")
if(length(grep("deirdre", getwd())>0)) {  setwd("~/Documents/github/Calvert_O-Connor_eelgrass")
} else {   setwd("~/Documents/git/temp")
}

# grazers
graz<- read.csv( "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv")
graz <- graz[,c(2:7)]

unique(graz$year) # 2014 2015 2016 2017
length(unique(graz$taxon)) # 194
sort(unique(graz$taxon))
# includes: 
#"Unknown.sea.star..I.think.only.1.sp."   
#"Unk.Fish.Eggs."

# some cleaning still needed:
# "Alia carinata"  "Alia.carinata" 
#"Amphissa"   merge with "Amphissa.columbiana" ?                                      
# "Amphitoidae" same as "Amphiuridae"?
# "Ampithoe.spp." and  "Ampithoidae" 
# "Aoridae" and "Aoroides" and "Aoroides.spp."
# "Caprella.californica" "Caprella californica"
#.....
# "Fluff"   "tin foil"  "foil"  "microfiber" --- like garbage?
# egg vs larvae
temp <- subset(graz, taxon == "") # 8 rows with no taxon

length(unique(graz$site)) # 10
# "choked_sandspit"
# "choked_inner" 
# "goose_south_west" 
# "goose_south_east"
# "goose_north"
# "mcmullins_south" 
# "mcmullins_north"  
# "triquet_north"    
# "triquet_south"          
# "pruth_pocket"    

hist(graz$size) # continuous data, 0.2 to 50.0

# Inverts 
invert <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_inverts_finest_1000_COVERAGE_RAREF.csv")

unique(invert$year) # 2014 2015 2016 2017
length(unique(invert$site)) #10
#"choked_inner"     "choked_sandspit"  "goose_north"      "goose_south_east" "goose_south_west"
# "mcmullins_north"  "mcmullins_south"  "pruth_pocket"     "triquet_north"    "triquet_south" 
# same 10 sites as grazers

temp <- invert[, c(8:112)]
sort(unique(colnames(temp))) # 105 unique species names---seem clean

# macro_eurkaryotes 18S
macroEuk <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv")

unique(macroEuk$year) # 2015 2016 2017 2018
length(unique(macroEuk$site)) #10

#"choked_inner"     "choked_sandspit"  "goose_south_east" "goose_south_west" "mcmullins_north"  "mcmullins_south" 
# "pruth_bay"        "pruth_pocket"     "triquet_north"    "triquet_south"
# missing goose_north; but added "pruth_bay" 

temp <- macroEuk[, c(10:214)]
sort(unique(colnames(temp))) # 205 unique ASV


microEuk <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv")

# some overlap in the ASV numbers eg 30,34
unique(microEuk$year) # 2015 2016 2017 2018
length(unique(microEuk$site)) #10

#"choked_inner"     "choked_sandspit"  "goose_south_east" "goose_south_west" "mcmullins_north"  "mcmullins_south" 
# "pruth_bay"        "pruth_pocket"     "triquet_north"    "triquet_south"
# again missing goose_north; but added "pruth_bay" 

temp <- microEuk[, c(10:1152)]
length(unique(colnames(temp))) # 1143 unique ASV

prokay <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv")

unique(microEuk$year) # 2015 2016 2017 2018
length(unique(microEuk$site)) #10

#"choked_inner"     "choked_sandspit"  "goose_south_east" "goose_south_west" "mcmullins_north"  "mcmullins_south" 
# "pruth_bay"        "pruth_pocket"     "triquet_north"    "triquet_south"
# again missing goose_north; but added "pruth_bay" 

temp <- microEuk[, c(10:1152)]
length(unique(colnames(temp))) # 1143 unique ASV


# bring in the updated data---updated data?
mtaxa_update <- read.csv( "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv" )

# data files referenced in the code but not in the repo:
microbes_18S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level.csv", header=T)
microbes_16S_ASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_ALL_RAREFIED.csv", header=T)

# Questions to look into:
# do we know what species the ASV are?
# Are the macro and the grazers the same?
prokASV <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv", header=T)

prokGenus <- read.csv("Data/R_Code_for_Data_Prep/master_data/MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv", header=T)

invertSite <- invert[, c("ID", "year", "region", "site")]
colnames(invertSite)[colnames(invertSite) == "ID"] <- "SampleID"
invertSite$data <- "invert"

macroEukSite <- macroEuk[, c("SampleID", "year", "region", "site")]
macroEukSite$data <- "macroEuk"

microEukSite <- microEuk[, c("SampleID", "year", "region", "site")]
microEukSite$data <- "microEuk"

prokaySite <- prokay[, c("SampleID", "year", "region", "site")]
prokaySite$data <- "prokay"

full <- rbind(invertSite, macroEukSite, microEukSite, prokaySite)
full <- unique(full)

species.study <- aggregate(dat.nodups ["doy"],
                           dat.nodups[c("studyid", "species")],
                           FUN = length)species.study <- aggregate(dat.nodups ["doy"],
                                                                   dat.nodups[c("studyid", "species")],
                                                                   FUN = length)

