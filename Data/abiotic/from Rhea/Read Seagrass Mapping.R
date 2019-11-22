# load libraries
library( tidyverse )

# read data
d <- read.delim( "seagrass.mapping.all_18s.txt" )


# select relevant columns
dsel <- d %>% 
  select( sample_id = X.SampleID, swab_id, project_name, date, region, site, host_species, host_type,
          survey_type, host_tissue = growth, quadrat_id = meso_quadrat_id, shoot_id = meso_shoot_id,
          shoot_length = meso_shoot_length, shoot_width = meso_shoot_width, number_of_blades = meso_number_of_blades,
          total_dry_wt = meso_total_dry_wt, blade_dry_wt = meso_blade_dry_wt, microepiphyte_wt = meso_microepiphyte_wt,
          macroepiphyte_wt = meso_macroepiphyte_wt, smithora_wt = meso_smithora_wt, bryozoans = meso_bryozoans,
          quadrat_density = meso_quadrat_density, quadrat_biomass = meso_quadrat_biomass, broken, smithora_pa = smithora,
          depth_mean = mean_depth, depth_ctd = ctd_depth, depth_collection = collection_depth, conductivity, temperature, 
          pressure = Pressure, par, chl = fluorometry_chlorophyll, turbidity, do = dissolved_oxygen, 
          do_concentration = dissolved_oxygen_concentration, salinity, speed_of_sound )


## merge abiotic data from 2016 and 2017 --  2018? 2015 data are from CTDs from Rhea
# read 2016 data
d16 <- read.csv( "Hakai YSI data 2016 (compiled) - Sheet1.csv" )
# select columns
d16sel <- d16 %>%
  select( date, site, coordinates, depth_collection = depth..m., do = do..mg.l., do_concentration = do...., 
          conductivity = conductivity..us.cm., salinity = salinity..ppt., temperature = temperature..c.,
          depth_mean = max.depth.at.site..m. ) # NOTE WE RENAMED THE DEPTH TO BE MEAN INSTEAD OF MAX

# read 2017 data
d17 <- read.csv( "Hakai abiotic data 2017 - Sheet1.csv" )
# select columns
d17sel <- d17 %>%
  select( site, date, depth_collection = depth.m., temperature = temp.c., pressure = pressure.mmhg.,
          conductivity = cond.uscm., salinity = salinity.ppt., ph )

# what will the 2018 data look like?

# merge data 