---
title: "_READMEmasterData"
output: html_document
---

## Data files in the Data/R_Code_for_Data_Prep/master_data/ folder

#### MASTER_merged_explanatory_20200214.csv
* Includes data on site location, conditions (bed area, depth), environmental conditions (temp, salinity, sparse data on pH), quadrat level seagrass metrics (LAI, shoot density, biomass).
* What is X---looks like row numbers mistakenly included?
* Meta-data with units?

#### MASTER_abiotic_20200214.csv.csv
* Data on environmental conditions, with salinity and temperature being best represented, but most parameters being patchy. Data merged across years.
* How is the data summarized, is the depth max depth and the temp at this max depth? What is the diff between site depth and depth?
* Data from 2014?

#### MASTER_grazers.csv
* Data from 2014-2017
* Data of sample, size, taxon for different sites across years
* What does sample and size mean? Organism size? Are samples different quadrats? Can they be pooled?
* Size values appear to be highly duplicated. Are these true values or an artifact of the data cleaning process?

#### MASTER_inverts_finest_1000_COVERAGE_RAREF.csv
* Data from 2014-2017
* Rows are unique samples, columns are mostly species names. Assume the values in the matrix are abundances
* Are any of the columns/rows just zeros?

#### MASTER_macroeuk18S_ASV_level_1000_COR_SING_COVERAGE_RAREF.csv
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