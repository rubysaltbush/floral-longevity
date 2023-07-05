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
# NOT WORKING! Error message:
# Problem with the API: HTTP Status 400 :(

# export taxa to csv, paste into TNRS webtool, then read in results
#readr::write_csv(longev_taxa, "data_output/longev_taxa.csv") # 1544 longevity taxa
#readr::write_csv(sym_taxa, "data_output/sym_taxa.csv")
# more than 5000 sym taxa, have to paste in 5k at a time

# used online TNRS version 5.1, https://tnrs.biendata.org/ , 
# and downloaded csv of best matches
# read these back in
longev_taxa <- readr::read_csv("data_input/longevityall_tnrs_result_best.csv")
sym_taxa <- readr::read_csv("data_input/symtaxa_tnrs_result1best.csv")

# okay looking through longevity taxa matching, many with 1 to 1 match
# or slight orthographic variations, only a few that appear incorrectly matched
# for some reason. Many taxa with supspecies and variants have only
# been matched at species level which should be fine, only mystery 
# is why Solanum nudum doesn't match when it's accepted in POWO
# also mysterious - why did it return only 1530 rows when I submitted
# 1544 distinct names?

# reduce taxa matching to few columns of interest then match this to 
# longevity & symmetry data

longev_taxa <- longev_taxa %>%
  dplyr::select(og_species_patch = Name_submitted, Overall_score, 
                Taxonomic_status, Accepted_name, Accepted_name_rank,
                Accepted_family) %>%
  dplyr::distinct()
longev_taxa$Accepted_genus <- gsub(" .*", "", longev_taxa$Accepted_name)

sym_taxa <- sym_taxa %>%
  dplyr::select(og_species_patch = Name_submitted, Accepted_name, 
                Accepted_name_rank, Accepted_family) %>%
  dplyr::distinct()
sym_taxa$Accepted_genus <- gsub(" .*", "", sym_taxa$Accepted_name)

# join taxonomy info back onto data
longevity_all <- longevity_all %>%
  dplyr::left_join(longev_taxa, by = "og_species_patch")

sym_data$og_species_patch <- sym_data$og_species

sym_data <- sym_data %>%
  dplyr::left_join(sym_taxa, by = "og_species_patch")

rm(longev_taxa, sym_taxa)

#### MATCH SYM AND LONGEV ####

# need to rename symmetry data columns to be unique from longevity data columns
# and simplify symmetry so that single value per accepted taxon
# FOR NOW JUST GET RID OF SOURCE COLUMN
sym_accepted <- sym_data %>%
  dplyr::group_by(Accepted_name, Accepted_name_rank, Accepted_genus, Accepted_family) %>%
  dplyr::summarise(sym_all = stringr::str_flatten(sort(unique(sym_all)), collapse = " "))
sym_accepted$sym_all <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_accepted$sym_all)
sym_accepted$sym_all <- gsub("zygomorphic zygomorphic", "zygomorphic", sym_accepted$sym_all)

# try just a straight join first
sym_long <- longevity_all %>%
  dplyr::left_join(sym_accepted, by = c("Accepted_name", "Accepted_name_rank",
                                    "Accepted_family", "Accepted_genus"))
# decent, would be good to match on genus as well

# longevity data at different ranks, can do multiple matches
table(longevity_all$Accepted_name_rank)


# first separate symmetry data into species, genus and family level data
sym_species <- sym_data %>%
  dplyr::filter(Accepted_name_rank %in% c("species", "subspecies", "variety"))
