
## Data files in the Data/R_Code_for_Data_Prep/master_data/ folder

## Data files
#### MASTER_merged_explanatory_20200214.csv
* Using data from MASTER_abiotic_20200214.csv, coordinates.csv, MASTER_seagrass_metrics_20200214.csv
* Includes data on site location, conditions (bed area, depth), environmental conditions (temp, salinity, sparse data on pH), quadrat level seagrass metrics (LAI, shoot density, biomass).
* What is X---looks like row numbers mistakenly included?
* Meta-data with units?

#### MASTER_abiotic_20200214.csv
* Created from data files: abiotic/hakai_abiotic_2015.csv, hakai_abiotic_2016.csv, corrected_hakai_abiotic_2017.csv, corrected_hakai_abiotic_2018.csv
* Code selects the deepest temp
* Data on environmental conditions, with salinity and temperature being best represented, but most parameters being patchy. Data merged across years.
* How is the data summarized, is the temp at this max depth? What is the diff between site depth and depth?
* No data from 2014?

#### MASTER_grazers.csv
* Combination of data from Data/macro_eukaryotes/: hakai_grazers_2014.csv, hakai_grazers_2015.csv (NA's converted to zeros), hakai_grazers_2016.csv (sieves no longer used, start of body size), hakai_grazers_2017.csv 
* There was an issue with the quadrat numbers not lining up, these were fixed based on the data in 2015_quadrat_sample_match_old_names.csv
* Data from 2014-2017
* Data of sample (quadrat), sieve size (sieve size, mm), taxon for different sites across years, size = body size
* Size values appear to be highly duplicated. Are these true values or an artifact of the data cleaning process?

#### MASTER_inverts_finest_1000_COVERAGE_RAREF.csv
* Created from MASTER_grazers_finest_to_phyloseq.csv, 
* Data from 2014-2017
* Rows are unique samples, columns are mostly species names. Assume the values in the matrix are abundances
* Are any of the columns/rows just zeros?

#### MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv
* Created from macro_18S/18S_allyears_unfiltered.RDS
* Data from 2015-2018
* Rows are unique samples, all from old leaves
* What do the column names mean? ASV?

#### MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv
* Data from 2015-2018
* Rows are unique samples, all from old leaves

#### MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv
* Data from 2015-2018
* Rows are unique samples, all from old leaves

#### MASTER_seagrass_metrics_20200214.csv
* Seems to be summary stats of seagrass metrics at the quadrat level
* How were these values generated?

#### MASTER_shoots.csv
* Data about seagrass shoots and blades (length, width, etc)
* Collected from 2015-2018

#### MASTER_quadrats.csv
* Unclear what this data is. Includes information about the seagrass (biomass, shoot count), but also columns that I don't recognize (drift, anchored, genus names, etc)
* Data has lots of NA's but also not counts, often to many decimal places

## R files:
#### merge_abiotic_seagrass_metrics_coordinates.R 
* Creates MASTER_merged_explanatory_20200214.csv 

##### abiotic_clean+merge.R
* Code creating MASTER_abiotic_20200214.csv

#### macroeuk_clean+merge.R 
* File creates MASTER_grazers.csv
* Based on header sounds like it was updated more recently than macroeuk_clean+merge_v2.R

#### inverts_pipeline_COVERAGE_BASED_RAREF.R
* Creates MASTER_inverts_finest_1000_COVERAGE_RAREF.csv
* Creates macro_eukaryotes/inverts_family_average_otu_tables_1000.csv
* Creates MASTER_inverts_family_1000_COVERAGE_RAREF.csv

#### FINAL_macroeuk18S_pipeline_ONLY_COVERAGE_BASED_RAREF_1000.R
* Creates MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv
* Creates macro_18S_average_otu_tables_1000.csv
* Creates macro_18S_average_otu_tables_3.csv
* Code outlines all the assumptions made, the removal of plants etc

### Cleaned names
"choked south pigu" <- "choked_inner"
"triquet bay south" <- "triquet_south"
"flat island" <- "choked_inner"
"sandspit" <- "choked_sandspit"
"triquet bay" <- "triquet_south"
"choked interior i5" <- "choked_inner"
"sand spit" <- "choked_sandspit"
"choked I5" <- "choked_inner"
"sandspit" <- "choked_sandspit"
"Lower" <- "Lower Choked"

