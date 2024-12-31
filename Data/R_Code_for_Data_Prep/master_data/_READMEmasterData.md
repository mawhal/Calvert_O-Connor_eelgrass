
## Data files in the Data/R_Code_for_Data_Prep/master_data/ folder. This folder contains 20 data file

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
* There was an issue with the quadrat numbers not lineing up, these were fixed based on the data in 2015_quadrat_sample_match_old_names.csv
* Not a matrix--long format of indiv sizes
* Data from 2014-2017
* Data of sample (quadrat), taxon for different sites across years, size = sieve size (sieve size, mm) or body size
* Size values appear to be highly duplicated---reflects the binning by sieve size

#### MASTER_grazers_family_to_phyloseq.csv
* Created in Data/macro_eukaryotes/save_data_frame_for_phyloseq.R as the family-level dataset
* Created using "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv" and "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv"
* Read into Data/macro_eukaryotes/inverts_pipeline_COVERAGE_BASED_RAREF.R as "inverts_family"

#### MASTER_grazers_finest_to_phyloseq.csv
* Created in Data/macro_eukaryotes/save_data_frame_for_phyloseq.R as the community dataset
* Created using "Data/R_Code_for_Data_Prep/master_data/MASTER_grazers.csv" and "R_Code_and_Analysis/output_data/O'Connor_hakai_seagrass_taxa_edit.csv"
* Read into Data/macro_eukaryotes/inverts_pipeline_COVERAGE_BASED_RAREF.R as "inverts"

#### MASTER_inverts_family_1000_COVERAGE_RAREF.csv
* Created in inverts_pipeline_COVERAGE_BASED_RAREF.R and Data/macro_eukaryotes/Entamoeba_inverts_pipeline_COVERAGE_BASED_RAREF.R
* Read into R_Code_and_Analysis/betadiversity/NMDS_plots_PERMANOVA_coverage_based_raref.R as "inverts_family"

#### MASTER_inverts_finest_1000_COVERAGE_RAREF.csv
* Created from MASTER_grazers_finest_to_phyloseq.csv, 
* Data from 2014-2017
* Rows are unique samples, columns are mostly species names. Assume the values in the matrix are abundances
* Are any of the columns/rows just zeros?

#### MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv
* Created from macro_18S/18S_allyears_unfiltered.RDS
* Data from 2015-2018
* Rows are unique samples, all from old leaves
* Amplicon Sequence Variants (ASV)--like OTU's but not grouped---each unique sequence kept separate
* What do the number mean? 1000 Corrected for singletons

#### MASTER_macroeuk18S_ASV_level_1000_NOT_COR_SING_COVERAGE_RAREF.csv
* Read into in R_Code_and_Analysis/betadiversity/tests_macroeuk18S/NMDS_plots_tests_macroeuk18S.R
* Unclear what code creates the file

#### MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv
* Data from 2015-2018
* Rows are unique samples, all from old leaves
* Created in Data/micro_eukaryotes/tests/microeuk_pipeline_ONLY_COVERAGE_BASED_RAREF.R
* Read into:
  - R_Code_and_Analysis/varpart/varpart_pRDA_2016_2017.R
  - Data/micro_eukaryotes/FINAL_microeuk_pipeline_ONLY_COVERAGE_BASED_RAREF_1000.R
  - R_Code_and_Analysis/alphadiversity/alpha_coverage_based_raref.R
  - R_Code_and_Analysis/betadiversity/Jaccard_NMDS_plots_cover_based_raref.R
  - R_Code_and_Analysis/gammadiversity/gamma_partitioning_coverage_based_raref.R
  - R_Code_and_Analysis/gammadiversity/choked_triquet_iNEXT_richness_estimate_ONLY_COVERAGE_BASED_RAREF.R
  - R_Code_and_Analysis/gammadiversity/iNEXT_richness_estimate_ONLY_COVERAGE_BASED_RAREF.R
  - R_Code_and_Analysis/betadiversity/NMDS_plots_PERMANOVA_coverage_based_raref.R

#### MASTER_microeuk_family_level_1000_COVERAGE_RAREF.csv
* Created in Data/micro_eukaryotes/FINAL_microeuk_pipeline_ONLY_COVERAGE_BASED_RAREF_1000.R from "Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds" (which is created in the same file) and "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv"
* Read into:
  - R_Code_and_Analysis/gammadiversity/gamma_partitioning_coverage_based_raref.R

#### MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv
* Created in Data/micro_eukaryotes/FINAL_microeuk_pipeline_ONLY_COVERAGE_BASED_RAREF_1000.R from "Data/micro_eukaryotes/all_years_18S_filtered_meso_Zos_ASV.rds" (which is created in the same file) and "Data/R_Code_for_Data_Prep/master_data/MASTER_microeuk_ASV_level_1000_COVERAGE_RAREF.csv"
* Read into:
  - R_Code_and_Analysis/betadiversity/LCBD_all.R
  - R_Code_and_Analysis/mantel/18S_dist_script_genus.R
  - R_Code_and_Analysis/varpart/varpart_pRDA_all_years.R
  - R_Code_and_Analysis/gammadiversity/gamma_partitioning_coverage_based_raref.R
  - R_Code_and_Analysis/distance_decay/dd_2016_all_together_REMOVE_ENV_EFFECT.R
  - R_Code_and_Analysis/distance_decay/dd_2016_all_together_WITHOUT_PRUTH.R
  - R_Code_and_Analysis/alphadiversity/alpha_coverage_based_raref.R
  - R_Code_and_Analysis/distance_decay/dd_2017_all_together.R
  - R_Code_and_Analysis/distance_decay/dd_2016_all_together.R
  - R_Code_and_Analysis/varpart/varpart_pRDA_2016_2017.R
  - R_Code_and_Analysis/gammadiversity/iNEXT_richness_estimate_ONLY_COVERAGE_BASED_RAREF.R

#### MASTER_prokary_ASV_level_1000_COVERAGE_RAREF.csv
* Data from 2015-2018 of prokaryotes
* Rows are unique samples, all from old leaves
* Created in: 
  - Data/prokaryotes/tests/prokary_pipeline_ONLY_COVERAGE_BASED_RAREF.R
  - Data/prokaryotes/entamoeba_prokary_pipeline_ONLY_COVERAGE_BASED_RAREF.R
  - Data/prokaryotes/FINAL_prokary_pipeline_COVERAGE_BASED_RAREF_1000.R
* Read into:
  - R_Code_and_Analysis/varpart/varpart_pRDA_2016_2017.R
  - R_Code_and_Analysis/betadiversity/Jaccard_NMDS_plots_cover_based_raref.R
  - R_Code_and_Analysis/alphadiversity/alpha_coverage_based_raref.R
  - R_Code_and_Analysis/gammadiversity/gamma_partitioning_coverage_based_raref.R
  - R_Code_and_Analysis/gammadiversity/choked_triquet_iNEXT_richness_estimate_ONLY_COVERAGE_BASED_RAREF.R
  - R_Code_and_Analysis/gammadiversity/iNEXT_richness_estimate_ONLY_COVERAGE_BASED_RAREF.R
  - R_Code_and_Analysis/betadiversity/NMDS_plots_PERMANOVA_coverage_based_raref.R

#### MASTER_prokary_family_level_1000_COVERAGE_RAREF.csv
* Created in Data/prokaryotes/FINAL_prokary_pipeline_COVERAGE_BASED_RAREF_1000.R from "Data/prokaryotes/all_years_16S_filtered_meso_Zos_ASV.rds"
* Read into R_Code_and_Analysis/gammadiversity/gamma_partitioning_coverage_based_raref.R, but commented out

#### MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv
* Data from 2015-2018
* Rows are unique samples, all from old leaves
* File used to calculate alpha div
* Unclear what code creates this file
* Read into:
  - R_Code_and_Analysis/betadiversity/LCBD_all.R
  - R_Code_and_Analysis/mantel/16S_dist_script_genus.R
  - R_Code_and_Analysis/gammadiversity/YEARS_gamma_partitioning_coverage_based_raref.R
  - R_Code_and_Analysis/varpart/varpart_pRDA_all_years.R
  - R_Code_and_Analysis/varpart/varpart_pRDA_2016_2017.R
  - R_Code_and_Analysis/distance_decay/dd_2016_all_together.R
  - R_Code_and_Analysis/distance_decay/dd_2017_all_together.R
  - R_Code_and_Analysis/alphadiversity/alpha_coverage_based_raref.R
  - R_Code_and_Analysis/distance_decay/dd_2016_all_together_WITHOUT_PRUTH.R
  - R_Code_and_Analysis/distance_decay/dd_2016_all_together_REMOVE_ENV_EFFECT.R
  - R_Code_and_Analysis/gammadiversity/gamma_partitioning_coverage_based_raref.R
  - R_Code_and_Analysis/gammadiversity/iNEXT_richness_estimate_ONLY_COVERAGE_BASED_RAREF.R

#### MASTER_seagrass_metrics_20200214.csv
* Seems to be summary stats of seagrass metrics at the quadrat level
* Code to calculate LAI and other metrics
* Created in Data/R_Code_for_Data_Prep/new_corrected_seagrass_metrics_LAI_clean+merge.R
* Combined data from:
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2015_eelgrass_quadrat_numbered.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2015_eelgrass_single_shoots_numbered.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2016_eelgrass_quadrat_biomass.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2016_eelgrass_quadrat_driftseaweed_na.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2016_eelgrass_single.shoots.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2017_Epiphyte_Data_corrected.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2017_Eelgrass_Morphometrics_corrected.csv"
  - "Data/seagrass_metrics/seagrass_corrected_files_used_final_merge_code/Hakai_2018_Eelgrass_Morphometrics_Macroalgae_Epiphytes_corrected.csv"
* Read into:
  - Data/R_Code_for_Data_Prep/merge_abiotic_seagrass_metrics_coordinates.R as eelgrass_metrics

#### MASTER_shoots.csv
* Data about seagrass shoots and blades (length, width, etc)
* Collected from 2015-2018
* Created in: Data/R_Code_for_Data_Prep/old/old_seagrass_metrics_clean+merge.R
* Read into:
  - R_Code_and_Analysis/sampling effort.R
  - R_Code_and_Analysis/eelgrass_summary.R

#### MASTER_quadrats.csv
* Unclear what this data is. Includes information about the seagrass (biomass, shoot count), but also columns that I don't recognize (drift, anchored, genus names, etc)
* Data has lots of NA's but also not counts, often to many decimal places
* Created in: Data/R_Code_for_Data_Prep/old/old_seagrass_metrics_clean+merge.R
* Read into:
  - R_Code_and_Analysis/sampling effort.R
  - R_Code_and_Analysis/macroeuk_sites_samplesize.R as quadrate level data
  - R_Code_and_Analysis/macroeuk_initial_patterns.R
  - R_Code_and_Analysis/macroeuk_mobr.R

#### unique_taxa_raw.csv
* Code to write this file in 3 files:  
  - Data/macro_eukaryotes/macroeuk_clean+merge.R
  - Data/macro_eukaryotes/macroeuk_clean+merge_ver2.R
  - Data/macro_eukaryotes/macroeuk_data_fromRhea/Hakai_seagrass_annual_mesograzer_clean+merge_ver2.R
* Combined data of annual hakai_grazers data (Data/macro_eukaryotes/hakai_grazers_2014.csv etc.)
* printed list of species names from taxon column 
  
## R files used to generate data files: 
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
* Code outlines all the assumptions made, the removal of plants etc.

#### alphadiversity/alpha_chao1_shannon_pielou.R
* The older of the two alpha R code files
* Uses data from: MASTER_prokary_genus_level_1000_COVERAGE_RAREF.csv, MASTER_microeuk_genus_level_1000_COVERAGE_RAREF.csv
* calculates chao1, shannon, pielou diversity for each sample

#### alphadiversity/alpha_coverage_based_raref.R

### Cleaned names
* "choked south pigu" <- "choked_inner"
* "triquet bay south" <- "triquet_south"
* "flat island" <- "choked_inner"
* "sandspit" <- "choked_sandspit"
* "triquet bay" <- "triquet_south"
* "choked interior i5" <- "choked_inner"
* "sand spit" <- "choked_sandspit"
* "choked I5" <- "choked_inner"
* "sandspit" <- "choked_sandspit"
* "Lower" <- "Lower Choked"

