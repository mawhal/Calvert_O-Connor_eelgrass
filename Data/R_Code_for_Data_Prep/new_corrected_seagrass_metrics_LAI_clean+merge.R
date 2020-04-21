# PROJECT: HAKAI SYNTHESIS PAPER
# PURPOSE: HAKAI EELGRASS MORPHOMETRICS (2015 - 2018 DATA)
# EELGRASS BIOMASS, DENSITY, DETRITUS, LAI, MICROEPIPHYTES, AND MACROALGAE AT THE QUADRAT LEVEL FOR EACH YEAR
# AUTHOR: EMILY ADAMCZYK
# DATE CREATED: 12 DECEMBER 2019
# DATE UPDATED: 5 FEBRUARY 2020
# DATE MODIFIED: 13 FEBRUARY 2020 *Bia multiplied 2017 and 2018 microepiphytes by 10 because we filtered 10% of the volume in those years


# libraries
library(tidyverse)

#########################################################
####################### 2015 ############################
#########################################################

# reading in hakai 2015 quadrat data
hakai2015_quad <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2015_eelgrass_quadrat_numbered.csv")

#View(hakai2015_quad)

# selecting columns I need
names(hakai2015_quad)

hakai2015_quad1 <- hakai2015_quad[-c(3, 8:10, 14)]

# making all algae types one column = macroalgae
hakai2015_quad1.2 <- hakai2015_quad1 %>% 
  group_by(site, sample.ID) %>% 
  mutate(quadrat_macroalgae_g = sum(macroepiphyte,
                                    drift.seaweed,
                                    rooted.seaweed,
                                    na.rm = T))


# renaming columns
names(hakai2015_quad1.2)

hakai2015_quad2 <- hakai2015_quad1.2 %>% 
  dplyr::rename(quadrat = sample.ID,
         quadrat_shoot_density = X..of.shoots.in.quadrat,
         live_foil_wt_g = foil.weight,
         live_wet_wt_g = foil..wet.weight,
         live_dry_wt_g = foil..dry.weight)

#View(hakai2015_quad2)

# selecting Hakai sites only
hakai2015_quad2$site

hakai2015_quad2_site <- hakai2015_quad2 %>% 
    filter(site == "sandspit" | 
           site == "choked_south_pigu" | 
           site == "mcmullins_south" | 
           site == "mcmullins north" | 
           site == "mcmullins nroth" |
           site == "goose north" | 
           site == "pruth bay south" |
           site == "pruth bay north" |
           site == "triquet south" |
           site == "goose west" |
           site == "goose east" |
           site == "triquet north" |
           site == "wolf inner ASU control")

#View(hakai2015_quad2_site)

# renaming site names
hakai2015_quad2_site$site <- as.character(hakai2015_quad2_site$site)

pp <- which(hakai2015_quad2_site$site == "pruth bay north") 
hakai2015_quad2_site$site[pp] <- "pruth_pocket"

pb <- which(hakai2015_quad2_site$site == "pruth bay south") 
hakai2015_quad2_site$site[pb] <- "pruth_bay"

ss <- which(hakai2015_quad2_site$site == "sandspit") 
hakai2015_quad2_site$site[ss] <- "choked_sandspit"

wo <- which(hakai2015_quad2_site$site == "wolf inner ASU control") 
hakai2015_quad2_site$site[wo] <- "choked_wolf"

gn <- which(hakai2015_quad2_site$site == "goose north") 
hakai2015_quad2_site$site[gn] <- "goose_north"

ge <- which(hakai2015_quad2_site$site == "goose east") 
hakai2015_quad2_site$site[ge] <- "goose_south_east"

gw <- which(hakai2015_quad2_site$site == "goose west") 
hakai2015_quad2_site$site[gw] <- "goose_south_west"

tn <- which(hakai2015_quad2_site$site == "triquet north") 
hakai2015_quad2_site$site[tn] <- "triquet_north"

ts <- which(hakai2015_quad2_site$site == "triquet south") 
hakai2015_quad2_site$site[ts] <- "triquet_south"

ic <- which(hakai2015_quad2_site$site == "choked_south_pigu") 
hakai2015_quad2_site$site[ic] <- "choked_inner"

mn <- which(hakai2015_quad2_site$site == "mcmullins north") 
hakai2015_quad2_site$site[mn] <- "mcmullins_north"

mn2 <- which(hakai2015_quad2_site$site == "mcmullins nroth") 
hakai2015_quad2_site$site[mn2] <- "mcmullins_north"



hakai2015_quad2_site$site

#############################################

# now checking out hakai 2015 epiphyte data
hakai2015_epi <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2015_eelgrass_single_shoots_numbered.csv")

#View(hakai2015_epi)

# selecting columns I need and renaming columns
names(hakai2015_epi)

hakai2015_epi1 <- hakai2015_epi[-c(1, 4, 8:12, 14, 16)]

names(hakai2015_epi1)

hakai2015_epi2 <- hakai2015_epi1 %>% 
  dplyr::rename(quadrat = sample.ID,
         epi_shoot_length_cm = shoot.length,
         epi_shoot_width_cm = shoot.width,
         epi_blade_num = blades,
         epi_microepiphyte_mg = microepiphyte,
         epi_macroalgae_g = smithora)

names(hakai2015_epi2)

# renaming sites
hakai2015_epi2$site <- as.character(hakai2015_epi2$site)

ss <- which(hakai2015_epi2$site == "choked_sandspit") 
hakai2015_epi2$site[ss] <- "choked_sandspit"

gse <- which(hakai2015_epi2$site == "goose_east") 
hakai2015_epi2$site[gse] <- "goose_south_east"

gsw <- which(hakai2015_epi2$site == "goose_west") 
hakai2015_epi2$site[gsw] <- "goose_south_west"

tb <- which(hakai2015_epi2$site == "triquet_south") 
hakai2015_epi2$site[tb] <- "triquet_south"

ic <- which(hakai2015_epi2$site == "choked_south_pigu") 
hakai2015_epi2$site[ic] <- "choked_inner"


#View(hakai2015_epi2)


################################################

# joining 2 dataframes together

hakai2015_join <- full_join(hakai2015_quad2_site, hakai2015_epi2, by = c("site", "quadrat"))

#View(hakai2015_join)


################################################

###################################################
## calculating LAI (directions from Keila Stark) ##
###################################################

# calculate length x width x # of blades x 2 (for both sides of blades) for each of the epiphyte shoots. This gives leaf area in cm2 for the singular epiphyte blade

names(hakai2015_join)

hakai2015_lai <- hakai2015_join %>% 
  group_by(site, quadrat) %>%
  mutate(epi_shoot_sa_cm = sum(epi_shoot_length_cm*epi_shoot_width_cm*epi_blade_num*2)) 

#View(hakai2015_lai)

# now that we have the surface area of the singular epiphyte blade as cm2, we can multiply the area with the total number of regular shoots in the quadrat to get the estimated leaf area of a 25cm x 25cm quadrat

hakai2015_lai <- hakai2015_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(epi_sa_quad_cm = quadrat_shoot_density*epi_shoot_sa_cm)


# Lastly, we divide by quadrat area (25cm x 25cm), or 625cm, to get LAI
hakai2015_lai <- hakai2015_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_lai = epi_sa_quad_cm/625)


# if we wanted to calculate shoot density per m^2, we would do this:
# hakai2015_density <- hakai2015_lai %>% group_by(site, quadrat) %>%mutate(density_m2 = shoot_num*16)

# calculating aboveground dry biomass per cm^2

hakai2015_biomass <- hakai2015_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_biomass_g = live_dry_wt_g - live_foil_wt_g)

#View(hakai2015_biomass)

# calculating quadrat level microepiphytes by multiplying single shoot epiphyte biomass by quadrat shoot density

hakai2015_microepi <- hakai2015_biomass %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_microepiphyte_mg = epi_microepiphyte_mg*quadrat_shoot_density)

#View(hakai2015_microepi)


# adding year 2015 to dataframe

hakai2015_microepi$year <- "2015"

# selecting relevant columns only
names(hakai2015_microepi)

hakai2015_data <- hakai2015_microepi[-c(4:9, 11:17)]

#View(hakai2015_data)




#########################################################
####################### 2016 ############################
#########################################################


# reading in hakai 2016 quadrat data
hakai2016_eg <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2016_eelgrass_quadrat_biomass.csv")

#View(hakai2016_eg1)

# removing row 49: Triquet 1 quadrat because there are two Triquet 1 quads and this particular one has incomplete data

hakai2016_eg1 <- hakai2016_eg[-c(59),]

# reading in the quadrat algae data
hakai2016_algae <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2016_eelgrass_quadrat_driftseaweed_na.csv")

#View(hakai2016_algae)

# Replace trace with NA
#hakai2016_algae[hakai2016_algae=="trace"] <- NA
#which(hakai2016_algae == "trace")

# omitting unnecessary columns
hakai2016_algae1 <- hakai2016_algae[-c(4, 6, 8:9)]

names(hakai2016_algae1)

# combining all algae biomass into one quadrat regardless of species type
hakai2016_algae2 <- hakai2016_algae1 %>% 
  group_by(site, sample.ID) %>% 
  summarise(sum_tin_weight = sum(as.numeric(drift.algae.foil.wt.), na.rm = T),
            sum_dry_bag_weight = sum(as.numeric(drift.algae.dry.wt...bag.wt.), na.rm= T)
  )

#View(hakai2016_algae2)

# joining drift algae to eelgrass quadrat dataset
hakai2016_quad <- full_join(hakai2016_eg1, hakai2016_algae2, by = c("site" = "site", "sample.ID" = "sample.ID"))

#View(hakai2016_quad)

######################################################

names(hakai2016_quad)

# adding together all macroalgae variables
hakai2016_quad1 <- hakai2016_quad %>% 
  group_by(site, sample.ID) %>% 
  mutate(quad_macroalgae_dry_g = sum(sum_dry_bag_weight,
                                epiphyte.bag..dry.wt.,
                                na.rm = TRUE)) %>% 
  mutate(quad_macroalgae_tin_g = sum(sum_tin_weight,
                                epiphyte.bag.wt.,
                                na.rm = TRUE))

names(hakai2016_quad1)

# omitting unncessary columns
hakai2016_quad1.2 <- hakai2016_quad1[-c(4, 7, 10, 12:20)]

# renaming columns
names(hakai2016_quad1.2)

hakai2016_quad2 <- hakai2016_quad1.2 %>% 
  dplyr::rename(quadrat = sample.ID,
         quadrat_shoot_density = X..of.shoots.in.quadrat,
         live_foil_wt_g = eelgrass.foil.weight,
         live_dry_wt_g = foil..dry.weight,
         quadrat_detritus_foil_wt_g = detritus.foil.weight,
         quadrat_detritus_dry_wt_g = detritus.foil..dry.weight)

#View(hakai2016_quad2)

# adding year 2016 to dataframe

hakai2016_quad2$year <- "2016"

# renaming sites
hakai2016_quad2$site

hakai2016_quad2$site <- as.character(hakai2016_quad2$site)

pp <- which(hakai2016_quad2$site == "pruth pocket") 
hakai2016_quad2$site[pp] <- "pruth_pocket"

pb <- which(hakai2016_quad2$site == "Pruth Bay") 
hakai2016_quad2$site[pb] <- "pruth_bay"

ss <- which(hakai2016_quad2$site == "Sandspit") 
hakai2016_quad2$site[ss] <- "choked_sandspit"

wo <- which(hakai2016_quad2$site == "wolf choked pass") 
hakai2016_quad2$site[wo] <- "choked_wolf"

fi <- which(hakai2016_quad2$site == "Flat Island") 
hakai2016_quad2$site[fi] <- "choked_inner"

gse <- which(hakai2016_quad2$site == "goose south east") 
hakai2016_quad2$site[gse] <- "goose_south_east"

gsw <- which(hakai2016_quad2$site == "Goose South West") 
hakai2016_quad2$site[gsw] <- "goose_south_west"

tn <- which(hakai2016_quad2$site == "Triquet North") 
hakai2016_quad2$site[tn] <- "triquet_north"

tb <- which(hakai2016_quad2$site == "Triquet Bay") 
hakai2016_quad2$site[tb] <- "triquet_south"

hakai2016_quad2$site




#############################################

# now checking out hakai 2016 epiphyte data
hakai2016_epi <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2016_eelgrass_single.shoots.csv")

#View(hakai2016_epi)

# selecting columns I need and renaming columns
names(hakai2016_epi)
hakai2016_epi1 <- hakai2016_epi[-c(1, 4, 6, 14:16, 20:24)]

names(hakai2016_epi1)

hakai2016_epi2 <- hakai2016_epi1 %>% 
  dplyr::rename(quadrat = sample.ID,
         epi_blade_num = blades,
         epi_shoot_length_cm = shoot.length,
         epi_shoot_width_cm = shoot.width,
         epi_GFC_filter_mg = GF.C.filter.weight,
         epi_GFC_filter_sample_mg = GF.C.with.dry.sample,
         epi_live_foil_wt_g = foil.weight,
         epi_live_wet_wt_g = foilwet.weight,
         epi_live_dry_wt_g = foil.dry.weight,
         epi_macroalgae_foil_g = smithora.foil.weight,
         epi_macroalgae_wet_g = smithora.wet.weight,
         epi_macroalgae_dry_g = smithora.dry.weight)

names(hakai2016_epi2)

# renaming sandspit
hakai2016_epi2$site <- as.character(hakai2016_epi2$site)

pp <- which(hakai2016_epi2$site == "pruth pocket") 
hakai2016_epi2$site[pp] <- "pruth_pocket"

pb <- which(hakai2016_epi2$site == "Pruth Bay") 
hakai2016_epi2$site[pb] <- "pruth_bay"

ss <- which(hakai2016_epi2$site == "Sandspit") 
hakai2016_epi2$site[ss] <- "choked_sandspit"

wo <- which(hakai2016_epi2$site == "wolf choked pass") 
hakai2016_epi2$site[wo] <- "choked_wolf"

fi <- which(hakai2016_epi2$site == "Flat Island") 
hakai2016_epi2$site[fi] <- "choked_inner"

gse <- which(hakai2016_epi2$site == "goose south east") 
hakai2016_epi2$site[gse] <- "goose_south_east"

gsw <- which(hakai2016_epi2$site == "Goose South West") 
hakai2016_epi2$site[gsw] <- "goose_south_west"

tn <- which(hakai2016_epi2$site == "Triquet North") 
hakai2016_epi2$site[tn] <- "triquet_north"

tb <- which(hakai2016_epi2$site == "Triquet Bay") 
hakai2016_epi2$site[tb] <- "triquet_south"

hakai2016_epi2$site

#View(hakai2016_quad2)

################################################

# joining 2 dataframes together

hakai2016_join <- full_join(hakai2016_quad2, hakai2016_epi2, by = c("site", "quadrat"))

#View(hakai2016_join)



###################################################

names(hakai2016_join)

# calculating aboveground dry biomass, detritus biomass, macroalgae biomass, and epi microepiphyte biomass per 625cm2
hakai2016_biomass <- hakai2016_join %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_biomass_g = sum(live_dry_wt_g - live_foil_wt_g),
         quadrat_detritus_g = sum(quadrat_detritus_dry_wt_g - quadrat_detritus_foil_wt_g),
         quadrat_macroalgae_g = sum(quad_macroalgae_dry_g - quad_macroalgae_tin_g),
         quadrat_microepiphyte_mg = sum((epi_GFC_filter_sample_mg - epi_GFC_filter_mg)*quadrat_shoot_density))

###################################################
## calculating LAI (directions from Keila Stark) ##
###################################################

# calculate length x width x # of blades x 2 (for both sides of blades) for each the epiphyte shoots

#View(hakai2016_biomass)

names(hakai2016_biomass)

hakai2016_lai <- hakai2016_biomass %>% 
  group_by(site, quadrat) %>%
  mutate(epi_shoot_sa_cm = sum(epi_shoot_length_cm*epi_shoot_width_cm*epi_blade_num*2)) 

#View(hakai2016_lai)

# now that we have the surface area of the singular epiphyte blade as cm2, we can multiply the area with the total number of regular shoots in the quadrat to get the estimated leaf area of a 25cm x 25cm quadrat
hakai2016_lai <- hakai2016_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(epi_sa_quad_cm = quadrat_shoot_density*epi_shoot_sa_cm)


# Lastly, we divide by quadrat area (25cm x 25cm), or 625cm, to get LAI
hakai2016_lai <- hakai2016_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_lai = epi_sa_quad_cm/625)


#View(hakai2016_lai)

# selecting relevant columns only
names(hakai2016_lai)

hakai2016_data <- hakai2016_lai[-c(1, 5:10, 12:22, 27:28)]

#View(hakai2016_data)




#########################################################
####################### 2017 ############################
#########################################################

# reading in hakai 2017 quadrat data
eg<- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2017_Eelgrass_Morphometrics_corrected.csv")

#View(eg)

# renaming column names
eg1 <- eg %>% 
  dplyr::rename(quadrat = quadrat..,
                quadrat_shoot_density = regular.shoot..,
                flowering_shoot = flowering.shoot..,
                live_paper_bag_wt_g = live.paper.bag.weight..g.,
                detritus_tin_wt_g = detritus.tin.weight..g.,
                epiphyte_bag_wt_g = epiphyte.bag.weight..g.,
                live_wet_wt_g = live.wet.weight..g.,
                detritus_wet_wt_g = detritus.wet.weight..g.,
                live_dry_wt_g = live.dry.weight..g.,
                detritus_dry_wt_g = detritus.dry.weight..g.,
                epiphyte_wet_wt_g = epi.algae.wet.weight..g.,
                epiphyte_dry_wt_g = epi.algae.dry.weight..g.,
                shoot_1_total_length_cm = shoot.1.length,
                shoot_1_sheath_length_cm = shoot.1.sheath.length,
                shoot_1_width_cm = shoot.1.width..cm.,
                shoot_1_blade_num = shoot.1...blades,
                shoot_2_total_length_cm = shoot.2.length,
                shoot_2_sheath_length_cm = shoot.2.sheath.length,
                shoot_2_width_cm = shoot.2.width..cm.,
                shoot_2_blade_num = shoot.2...blades,
                shoot_3_total_length_cm = shoot.3.length,
                shoot_3_sheath_length_cm = shoot.3.sheath.length,
                shoot_3_width_cm = shoot.3.width..cm.,
                shoot_3_blade_num = shoot.3...blades,
                shoot_4_total_length_cm = shoot.4.length,
                shoot_4_sheath_length_cm = shoot.4.sheath.length,
                shoot_4_width_cm = shoot.4.width..cm.,
                shoot_4_blade_num = shoot.4...blades,
                shoot_5_total_length_cm = shoot.5.length,
                shoot_5_sheath_length_cm = shoot.5.sheath.length,
                shoot_5_width_cm = shoot.5.width..cm.,
                shoot_5_blade_num = shoot.5...blades,
                macroalgae_spp_1 = macroalgae.species.1,
                macroalgae_spp_1_wet_wt_g = macroalgae.1.wet.weight..g.,
                macroalgae_spp_1_tin_wt_g = macroalgae.1.foil..g.,
                macroalgae_spp_1_dry_wt_g = macroalgae.1.dry..g.,
                macroalgae_spp_2 = macroalgae.species.2,
                macroalgae_spp_2_wet_wt_g = macroalgae.2.wet.weight..g.,
                macroalgae_spp_2_tin_wt_g = macroalgae.2.foil..g.,
                macroalgae_spp_2_dry_wt_g = macroalgae.2.dry..g.,
                macroalgae_spp_3 = macroalgae.species.3,
                macroalgae_spp_3_wet_wt_g = macroalgae.3.wet.weight..g.,
                macroalgae_spp_3_tin_wt_g = macroalgae.3.foil..g.,
                macroalgae_spp_3_dry_wt_g = macroalgae.3.dry..g.,
                macroalgae_spp_4 = macroalgae.species.4,
                macroalgae_spp_4_wet_wt_g = macroalgae.4.wet.weight..g.,
                macroalgae_spp_4_tin_wt_g = macroalgae.4.foil..g.,
                macroalgae_spp_4_dry_wt_g = macroalgae.4.dry..g.
  )

# renaming site names
eg1$site <- as.character(eg1$site)

pp <- which(eg1$site == "pruth pocket") 
eg1$site[pp] <- "pruth_pocket"

pb <- which(eg1$site == "pruth bay") 
eg1$site[pb] <- "pruth_bay"

ss <- which(eg1$site == "sand spit") 
eg1$site[ss] <- "choked_sandspit"

i5 <- which(eg1$site == "inner choked i5") 
eg1$site[i5] <- "choked_inner"

tn <- which(eg1$site == "triquet north") 
eg1$site[tn] <- "triquet_north"

tb <- which(eg1$site == "triquet bay") 
eg1$site[tb] <- "triquet_south"

#View(eg1)

# subtracting off dry weight from bag/tin weight to get above ground biomass weight in grams
eg2 <- eg1 %>%
  group_by(date, site, quadrat) %>% 
  mutate(quadrat_biomass_g = sum(live_dry_wt_g - live_paper_bag_wt_g),
         quadrat_detritus_g = sum(detritus_dry_wt_g - detritus_tin_wt_g),
         macroalgae_dry_wt_g = sum(macroalgae_spp_1_dry_wt_g, 
                                   macroalgae_spp_2_dry_wt_g, 
                                   macroalgae_spp_3_dry_wt_g,
                                   macroalgae_spp_4_dry_wt_g,
                                   epiphyte_dry_wt_g,
                                   na.rm = T),
         macroalgae_tin_wt_g = sum(macroalgae_spp_1_tin_wt_g, 
                                   macroalgae_spp_2_tin_wt_g, 
                                   macroalgae_spp_3_tin_wt_g,
                                   macroalgae_spp_4_tin_wt_g,
                                   epiphyte_bag_wt_g,
                                   na.rm = T)
           
           )

eg2.2 <- eg2 %>% 
  group_by(date, site, quadrat) %>% 
  mutate(quadrat_macroalgae_g = sum(macroalgae_dry_wt_g - macroalgae_tin_wt_g))

# getting rid of unncessary columns
names(eg2.2)

eg3 <- eg2.2[-c(1, 5:18, 20, 24, 27:52, 55:56)]

# adding year 2017 to dataframe
eg3$year <- "2017"

#View(eg3)

#############################################################################################

# now checking out hakai 2017 epiphyte data
hakai2017_epi <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2017_Epiphyte_Data_corrected.csv")

# removing unnecessary columns
names(hakai2017_epi)

hakai2017_epi2 <- hakai2017_epi[-c(1, 7:10, 13:50)]

# renaming columns
names(hakai2017_epi2)

hakai2017_epi3 <- hakai2017_epi2 %>% 
  dplyr::rename(quadrat = quadrat..,
         epi_blade_num = epi.blade..,
         epi_shoot_length_cm = epi.blade.length..cm.,
         epi_shoot_width_cm = epi.blade.width..cm.,
         GFC_wt_g = GF.C.weight..g.,
         GFC_wt_sample_g = GF.C.with.dry.sample..g.)

#View(hakai2017_epi3)

# renaming site names
hakai2017_epi3$site <- as.character(hakai2017_epi3$site)

pp <- which(hakai2017_epi3$site == "pruth pocket") 
hakai2017_epi3$site[pp] <- "pruth_pocket"

pb <- which(hakai2017_epi3$site == "pruth bay") 
hakai2017_epi3$site[pb] <- "pruth_bay"

ss <- which(hakai2017_epi3$site == "sand spit") 
hakai2017_epi3$site[ss] <- "choked_sandspit"

i5 <- which(hakai2017_epi3$site == "inner choked i5") 
hakai2017_epi3$site[i5] <- "choked_inner"

tn <- which(hakai2017_epi3$site == "triquet north") 
hakai2017_epi3$site[tn] <- "triquet_north"

tb <- which(hakai2017_epi3$site == "triquet bay") 
hakai2017_epi3$site[tb] <- "triquet_south"


################################################

# joining 2 dataframes together

hakai2017_join <- full_join(eg3, hakai2017_epi3, by = c("site", "quadrat"))

#View(hakai2017_join)

#############################################################################################

###################################################
## calculating LAI (directions from Keila Stark) ##
###################################################

# calculate length x width x # of blades x 2 (for both sides of blades) for each the epiphyte shoots

#View(hakai2017_join)

names(hakai2017_join)

hakai2017_lai <- hakai2017_join %>% 
  group_by(site, quadrat) %>%
  mutate(epi_shoot_sa_cm = sum(epi_shoot_length_cm*epi_shoot_width_cm*epi_blade_num*2)) 

#View(hakai2017_lai)

# now that we have the surface area of the singular epiphyte blade as cm2, we can multiply the area with the total number of regular shoots in the quadrat to get the estimated leaf area of a 25cm x 25cm quadrat
hakai2017_lai <- hakai2017_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(epi_sa_quad_cm = quadrat_shoot_density*epi_shoot_sa_cm)


# Lastly, we divide by quadrat area (25cm x 25cm), or 625cm, to get LAI
hakai2017_lai <- hakai2017_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_lai = epi_sa_quad_cm/625)


#View(hakai2017_lai)

#########

# pruth bay, triquet north, and triquet bay do not have LAI because the blade number was not counted for the epiphyte shoots. I will use shoot 2 from the quadrat data as a substitute for those sites
other_lai <- hakai2017_lai %>% 
  filter(site %in% c("pruth_bay", "triquet_south", "triquet_north"))

# calculate length x width x # of blades x 2 (for both sides of blades) for shoot 2
names(other_lai)
other_lai2 <- other_lai %>% 
  group_by(site, quadrat) %>%
  mutate(shoot_2_sa_cm = sum(shoot_2_total_length_cm*shoot_2_width_cm*shoot_2_blade_num*2)) 

#View(other_lai2)

# now that we have the surface area of the shoot 3 as cm2, we can multiply the area with the total number of regular shoots in the quadrat to get the estimated leaf area of a 25cm x 25cm quadrat
other_lai3 <- other_lai2 %>% 
  group_by(site, quadrat) %>% 
  mutate(other_sa_quad_cm = quadrat_shoot_density*shoot_2_sa_cm)


# Lastly, we divide by quadrat area (25cm x 25cm), or 625cm, to get LAI
other_lai4 <- other_lai3 %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_lai = other_sa_quad_cm/625)

# Now, using shoot number 3 for Triquet North quadrat 2 LAI because shoot number 2 did not have a blade count number
names(other_lai4)

tn <- other_lai4 %>% 
  filter(site == "triquet_north",
         quadrat == "2")

#calculate length x width x # of blades x 2 (for both sides of blades) for shoot 2
names(tn)
tn2 <- tn %>% 
  group_by(site, quadrat) %>%
  mutate(shoot_3_sa_cm = sum(shoot_3_total_length_cm*shoot_3_width_cm*shoot_3_blade_num*2)) 

#View(tn2)

# now that we have the surface area of the shoot 3 as cm2, we can multiply the area with the total number of regular shoots in the quadrat to get the estimated leaf area of a 25cm x 25cm quadrat
tn3 <- tn2%>% 
  group_by(site, quadrat) %>% 
  mutate(other_sa_quad_cm = quadrat_shoot_density*shoot_3_sa_cm)


# Lastly, we divide by quadrat area (25cm x 25cm), or 625cm, to get LAI
tn4 <- tn3 %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_lai = other_sa_quad_cm/625)


# Joining back with other data
hakai2017_lai2 <- rbind (hakai2017_lai, tn4, other_lai4)

# removing unnecessary rows
#View(hakai2017_lai2)

hakai2017_lai3 <- hakai2017_lai2[-c(19:27, 37:54, 66), ]

#View(hakai2017_lai3)

# Calculating microepiphytes at the quadrat level
hakai2017_microepi <- hakai2017_lai3 %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_microepiphyte_g = sum(GFC_wt_sample_g - GFC_wt_g)*quadrat_shoot_density)

# converting microepiphyes from grams to mg
hakai2017_microepi2 <- hakai2017_microepi %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_microepiphyte_mg = sum(quadrat_microepiphyte_g*1000))

# transforming the values to 100% (we only filtered 10% of the volume)
hakai2017_microepi2 <- hakai2017_microepi2 %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_microepiphyte_mg = sum(quadrat_microepiphyte_mg*10))

# selecting relevant columns only
names(hakai2017_microepi2)

hakai2017_data <- hakai2017_microepi2 [-c(4:9, 14:20, 22:24, 25)]

#View(hakai2017_data)

######################################################



#########################################################
####################### 2018 ############################
#########################################################

# reading in hakai 2018 quadrat data
eg18 <- read.csv("Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2018_Eelgrass_Morphometrics_Macroalgae_Epiphytes_corrected.csv")

#View(eg18)

# selecting relevant columns and renaming column names
names(eg18)

eg18_1 <- eg18[-c(2, 5, 8:9, 12:32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 55:58, 61:103)]

names(eg18_1)

# renaming some columns

eg18_1.1 <- eg18_1 %>% 
  dplyr::rename(site = `site`,
         quadrat = quadrat..,
         quadrat_shoot_density = regular.shoot..,
         epi_blade_num = epi.blade..,
         epi_shoot_length_cm = epi.blade.length..cm.,
         epi_shoot_width_cm = epi.blade.width..cm.)

# combining all macroalgae columns into one
# also subtracting off dry weight from bag/tin weight to get above ground biomass weight in grams
eg18_2 <- eg18_1.1 %>%
  group_by(site, quadrat) %>% 
  mutate(quadrat_biomass_g = sum(live.dry.weight..g. - live.paper.bag.weight..g.),
         quadrat_detritus_g = sum(detritus.dry.weight..g. - detritus.tin.weight..g.),
         epi_microepiphyte_g = sum(GF.C.with.dry.sample..g. - GF.C.weight..g.),
         macroalgae_dry_wt_g = as.numeric(sum(macroalgae.dry.1..g., 
                                   macroalgae.dry.2..g.,
                                   macroalgae.dry.3..g.,
                                   macroalgae.dry.4..g.,
                                   macroalgae.dry.5..g.,
                                   na.rm = T)),
         macroalgae_tin_wt_g = sum(macroalgae.foil1..g.,
                                   macroalgae.foil.2..g.,
                                   macroalgae.foil.3..g.,
                                   macroalgae.foil.4..g.,
                                   macroalgae.foil.5..g.,
                                   na.rm = T)
  )

#View(eg18_2)

# Calculating quadrat macroalgae biomass and calculating microepiphytes to reflect quadrat level microepi biomass
eg18_2.2 <- eg18_2 %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_macroalgae_g = sum(macroalgae_dry_wt_g - macroalgae_tin_wt_g),
         quadrat_microepiphyte_mg = sum(epi_microepiphyte_g*1000))

# transforming the values to 100% (we only filtered 10% of the volume)
eg18_2.2 <- eg18_2.2 %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_macroalgae_g = sum(macroalgae_dry_wt_g - macroalgae_tin_wt_g),
         quadrat_microepiphyte_mg = sum(quadrat_microepiphyte_mg*10))


# Omitting unnecessary columns
names(eg18_2.2)

eg18_3 <- eg18_2.2[-c(4:17, 21:22)]

names(eg18_3)

# renaming site names
eg18_3$site <- as.character(eg18_3$site)

pp <- which(eg18_3$site == "PruthPocket") 
eg18_3$site[pp] <- "pruth_pocket"

pb <- which(eg18_3$site == "Pruth Bay") 
eg18_3$site[pb] <- "pruth_bay"

ss <- which(eg18_3$site == "Sandspit") 
eg18_3$site[ss] <- "choked_sandspit"

i5 <- which(eg18_3$site == "Choked I5") 
eg18_3$site[i5] <- "choked_inner"

tb <- which(eg18_3$site == "TriquetBay") 
eg18_3$site[tb] <- "triquet_south"

tn <- which(eg18_3$site == "TriquetNorth") 
eg18_3$site[tn] <- "triquet_north"

eg18_3$site

#############################################################################################

###################################################
## calculating LAI (directions from Keila Stark) ##
###################################################

# calculate length x width x # of blades x 2 (for both sides of blades) for each the epiphyte shoots

#View(eg18_3)

names(eg18_3)

hakai2018_lai <- eg18_3 %>% 
  group_by(site, quadrat) %>%
  mutate(epi_shoot_sa_cm = sum(epi_shoot_length_cm*epi_shoot_width_cm*epi_blade_num*2)) 

#View(hakai2018_lai)

# now that we have the surface area of the singular epiphyte blade as cm2, we can multiply the area with the total number of regular shoots in the quadrat to get the estimated leaf area of a 25cm x 25cm quadrat
hakai2018_lai <- hakai2018_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(epi_sa_quad_cm = quadrat_shoot_density*epi_shoot_sa_cm)


# Lastly, we divide by quadrat area (25cm x 25cm), or 625cm, to get LAI
hakai2018_lai <- hakai2018_lai %>% 
  group_by(site, quadrat) %>% 
  mutate(quadrat_lai = epi_sa_quad_cm/625)

# adding year 2018 to dataframe
hakai2018_lai$year <- "2018"

# selecting relevant columns only
names(hakai2018_lai)

hakai2018_data <- hakai2018_lai [-c(4:6, 9:11, 14:15)]

#View(hakai2018_data)


##########################################################################################
##########################################################################################
##########################################################################################


#########################################
### combining all dataframes together ###
#########################################

hakai_synthesis <- rbind(hakai2015_data, hakai2016_data, hakai2017_data, hakai2018_data) #nrow = 261

#View(hakai_synthesis)

# renaming quadrat to 

hakai_synthesis_1 <-  hakai_synthesis %>% 
  rename(quadrat_id = quadrat)

View(hakai_synthesis_1)

write.csv(hakai_synthesis_1, "Data/R_Code_for_Data_Prep/master_data/MASTER_seagrass_metrics_20200214.csv")

# double checking that all rows were adding together
nrow(hakai2015_data) #71
nrow(hakai2016_data) #81
nrow(hakai2017_data) #54
nrow(hakai2018_data) #54



######################################################################################

## YEET! Let's make some dank plots

#### Environment Settings ####
# Set global ggplot2 theme; this is most common ggplot format for publication
theme_set(theme_bw())

###########################

# BIOMASS
names(hakai_synthesis)

# BIOMASS SITE
biomass_site <- ggplot(hakai_synthesis, aes(x = year, y = quadrat_biomass_g, color = year)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.4) +
  facet_wrap(~site, scales = "free") +
  theme (#axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 12), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         #legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 12), #font size of legend
         plot.title = element_text(hjust = 0.5, size = 16)) + #center plot title and set font size
  xlab("") +
  ylab("aboveground dry biomass per 625cm^2 (g)")


biomass_site
ggsave("hakai_eelgrass_biomass_site.png", biomass_site, width = 8, height = 6)


# BIOMASS YEAR
biomass_year <- ggplot(hakai_synthesis, aes(x = year, y = quadrat_biomass_g, color = site)) +
  geom_boxplot() +
  facet_wrap(~year, scales = "free") +
  theme (#axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 12), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         #legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 12), #font size of legend
         plot.title = element_text(hjust = 0.5, size = 16)) + #center plot title and set font size
  xlab("") +
  ylab("aboveground dry biomass per 625cm^2 (g)")


biomass_year

ggsave("hakai_eelgrass_biomass_year.png", biomass_year, width = 8, height = 6)



##################

# DENSITY

# DENSITY SITE
density_site <- ggplot(hakai_synthesis, aes(x = year, y = quadrat_shoot_density, color = year)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.4) +
  facet_wrap(~site, scales = "free") +
  theme (#axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 12), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         #legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 12), #font size of legend
         plot.title = element_text(hjust = 0.5, size = 16)) + #center plot title and set font size
  xlab("") +
  ylab("Eelgrass density (number of shoots per 625cm^2)")


density_site

ggsave("hakai_eelgrass_density_site.png", density_site, width = 8, height = 6)

# DENSITY YEAR

density_year <- ggplot(hakai_synthesis, aes(x = year, y = quadrat_shoot_density, color = site)) +
  geom_boxplot() +
  facet_wrap(~year, scales = "free") +
  theme (#axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.text.x = element_blank(),
         axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 12), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         #legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 12), #font size of legend
         plot.title = element_text(hjust = 0.5, size = 16)) + #center plot title and set font size
  xlab("") +
  ylab("Eelgrass density (number of shoots per 625cm^2)")


density_year

ggsave("hakai_eelgrass_density_year.png", density_year, width = 8, height = 6)



##################

# LAI

# LAI SITE
lai_site <- ggplot(hakai_synthesis, aes(x = year, y = quadrat_lai, color = year)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.4) +
  facet_wrap(~site, scales = "free") +
  theme (axis.text.x = element_blank(),
         #axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 12), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         #legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 12), #font size of legend
         plot.title = element_text(hjust = 0.5, size = 16)) + #center plot title and set font size
  xlab("") +
  ylab("Eelgrass LAI 625cm^2")


lai_site


ggsave("hakai_eelgrass_lai_site.png", lai_site, width = 8, height = 6)


# LAI YEAR

lai_year <- ggplot(hakai_synthesis, aes(x = year, y = quadrat_lai, color = site)) +
  geom_boxplot() +
  facet_wrap(~year, scales = "free") +
  theme (axis.text.x = element_blank(),
         #axis.title.x = element_text(size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
         axis.title.y = element_text(size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
         axis.text = element_text(size = 12), #font size of numbers in axis
         panel.grid.major = element_blank(), #remove major grid
         panel.grid.minor = element_blank(), #remove minor grid
         axis.line = element_line(colour = "black"), #draw line in the axis
         panel.border = element_blank(), #remove lines outside the graph
         #legend.title=element_blank(), #remove legend title
         legend.direction = "vertical", #direction
         legend.justification = c(1, 1), legend.position = "right", #legend is top right
         legend.key.size = unit(2.0, 'lines'), #spacing between legends
         legend.text = element_text(size = 12), #font size of legend
         plot.title = element_text(hjust = 0.5, size = 16)) + #center plot title and set font size
  xlab("") +
  ylab("Eelgrass LAI 625cm^2")


lai_year


ggsave("hakai_eelgrass_lai_year.png", lai_year, width = 8, height = 6)
