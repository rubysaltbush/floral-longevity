# script to match taxonomy between symmetry and longevity data, and produce 
# clean taxonomically matched data set

#### GET DATA ####

# source symmetry and longevity data from respective scripts
source("scripts/symmetry_data.R")
source("scripts/longevity_data.R")

### MATCH TAXONOMY ####

# match taxonomy to World Flora Online using TNRS R package
# first reduce to taxa
longev_taxa <- longevity_all %>% 
  dplyr::select(og_species_patch) %>%
  dplyr::distinct()

sym_taxa <- sym_data %>%
  dplyr::select(og_species) %>%
  dplyr::distinct()

# resolve to World Flora Online using TNRS
#longev_tnrs <- TNRS(longev_taxa)
#sym_tnrs <- TNRS(sym_taxa)
# Problem with the API: HTTP Status 400 :(

# export taxa to csv, paste into TNRS webtool, then read in results
#readr::write_csv(longev_taxa, "data_output/longev_taxa.csv")
#readr::write_csv(sym_taxa, "data_output/sym_taxa.csv")
# more than 5000 sym taxa, have to paste in 5k at a time

# used online TNRS version 5.1, https://tnrs.biendata.org/ , 
# and downloaded csv of best matches
# read these back in
longev_taxa <- readr::read_csv("data_input/longevityall_tnrs_result_best.csv")
sym_taxa <- readr::read_csv("data_input/symtaxa_tnrs_result1best.csv")



# NEXXT - match both to accepted taxa, then match symmetry data 
# to longevity data, then export for filling in 

