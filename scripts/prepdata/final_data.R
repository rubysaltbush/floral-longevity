# script to match taxonomy between symmetry and longevity data, and produce 
# clean taxonomically matched data set

sym_long <- cache_csv("data_output/longevity_symmetry_all.csv", function() {

#### GET DATA ####

# source symmetry and longevity data from respective scripts
source("scripts/prepdata/symmetry_data.R")
source("scripts/prepdata/longevity_data.R")

### MATCH TAXONOMY ####

# match taxonomy to World Flora Online using TNRS R package

# longevity first
longev_tnrs <- cache_csv("data_output/longev_taxa.csv", function() {
  # first reduce to taxa
  longev_taxa <- longevity_all %>% 
    dplyr::select(og_species_patch) %>%
    dplyr::distinct()
  # resolve to World Flora Online using TNRS
  longev_tnrs <- TNRS::TNRS(longev_taxa$og_species_patch)
  # export matches to csv to cache
  readr::write_csv(longev_tnrs, "data_output/longev_taxa.csv")
})

# then symmetry (will take longer)
sym_tnrs <- cache_csv("data_output/sym_taxa.csv", function() {
 # first reduce to taxa
 sym_taxa <- sym_data %>%
   dplyr::select(og_species) %>%
   dplyr::distinct()
 # resolve to World Flora Online using TNRS
 sym_tnrs <- TNRS::TNRS(sym_taxa$og_species)
 # export matches to csv to cache
 readr::write_csv(sym_tnrs, "data_output/sym_taxa.csv")
})

# have checked through all matches <1 in longevity data, agree with TNRS for these
# for symmetry will filter out matches <0.5 as some of these iffy
sym_tnrs <- sym_tnrs %>%
  dplyr::filter(Overall_score >= 0.5)

# reduce taxa matching to few columns of interest then match this to 
# longevity & symmetry data

longev_tnrs <- longev_tnrs %>%
  dplyr::select(og_species_patch = Name_submitted, Accepted_name, 
                Accepted_name_rank, Accepted_family) %>%
  dplyr::distinct()
longev_tnrs$Accepted_genus <- gsub(" .*", "", longev_tnrs$Accepted_name)

sym_tnrs <- sym_tnrs %>%
  dplyr::select(og_species_patch = Name_submitted, Accepted_name, 
                Accepted_name_rank, Accepted_family) %>%
  dplyr::distinct()
sym_tnrs$Accepted_genus <- gsub(" .*", "", sym_tnrs$Accepted_name)

# join taxonomy info back onto longevity data
longevity_all <- longevity_all %>%
  dplyr::left_join(longev_tnrs, by = "og_species_patch")
# patch Solanum nudum which is accepted according to POWO, not sure why it didn't match
longevity_all$Accepted_name[longevity_all$og_species_patch == "Solanum nudum"] <- "Solanum nudum"
longevity_all$Accepted_name_rank[longevity_all$og_species_patch == "Solanum nudum"] <- "species"

# join taxonomy info back onto symmetry data
sym_data$og_species_patch <- sym_data$og_species
sym_data <- sym_data %>%
  dplyr::left_join(sym_tnrs, by = "og_species_patch")

rm(longev_tnrs, sym_tnrs)

#### MATCH SYM AND LONGEV ####

# need to rename symmetry data columns to be unique from longevity data columns
# and simplify symmetry so that single value per accepted taxon
# for now will lose sources of symmetry data but can trace this back later
sym_accepted <- sym_data %>%
  dplyr::group_by(Accepted_name, Accepted_name_rank, Accepted_genus, Accepted_family) %>%
  dplyr::summarise(sym_all = stringr::str_flatten(sort(unique(sym_all)), collapse = " ")) %>%
  dplyr::ungroup()
sym_accepted$sym_all <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_accepted$sym_all)
sym_accepted$sym_all <- gsub("zygomorphic zygomorphic", "zygomorphic", sym_accepted$sym_all)
table(sym_accepted$sym_all)

# longevity data at different ranks, can do multiple matches
table(longevity_all$Accepted_name_rank)

# first separate symmetry data into different taxonomic ranks
# variety and subspecies
sym_var_ssp <- sym_accepted %>%
  dplyr::filter(Accepted_name_rank %in% c("variety", "subspecies")) %>%
  dplyr::select(Accepted_name, sym_sspvar = sym_all)
# species
sym_species <- sym_accepted %>%
  dplyr::filter(Accepted_name_rank %in% c("species", "subspecies", "variety"))
sym_species$Accepted_name <- gsub(" var\\..*| subsp\\..*", "", sym_species$Accepted_name)
sym_species <- sym_species %>%
  dplyr::group_by(Accepted_name, Accepted_genus, Accepted_family) %>%
  dplyr::summarise(sym_species = stringr::str_flatten(sort(unique(sym_all)), collapse = " "))
sym_species$sym_species <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_species$sym_species)
sym_species$sym_species <- gsub("zygomorphic zygomorphic", "zygomorphic", sym_species$sym_species)
table(sym_species$sym_species)
# genus
sym_genus <- sym_species %>%
  dplyr::group_by(Accepted_genus, Accepted_family) %>%
  dplyr::summarise(sym_genus = stringr::str_flatten(sort(unique(sym_species)), collapse = " "))
sym_genus$sym_genus <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_genus$sym_genus)
sym_genus$sym_genus <- gsub("zygomorphic zygomorphic", "zygomorphic", sym_genus$sym_genus)
table(sym_genus$sym_genus)
# family
sym_family <- sym_genus %>%
  dplyr::group_by(Accepted_family) %>%
  dplyr::summarise(sym_family = stringr::str_flatten(sort(unique(sym_genus)), collapse = " "))
sym_family$sym_family <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_family$sym_family)
sym_family$sym_family <- gsub("zygomorphic zygomorphic", "zygomorphic", sym_family$sym_family)
table(sym_family$sym_family)

# now join these various instances of symmetry data to longevity data
sym_long <- longevity_all %>%
  dplyr::left_join(sym_var_ssp, by = "Accepted_name") %>%
  dplyr::left_join(sym_species, by = c("Accepted_name", "Accepted_genus", "Accepted_family")) %>%
  dplyr::left_join(sym_genus, by = c("Accepted_genus", "Accepted_family")) %>%
  dplyr::left_join(sym_family, by = c("Accepted_family"))
rm(sym_family, sym_genus, sym_species, sym_var_ssp, sym_accepted, sym_data, longevity_all)

# paste var and subsp symmetry into sym_species column to simplify
sym_long$sym_species <- ifelse(is.na(sym_long$sym_species), 
                               sym_long$sym_sspvar, 
                               sym_long$sym_species)
sym_long <- dplyr::select(sym_long, -sym_sspvar)
# and paste RS scored symmetry into sym_species column also
sym_long$sym_species <- ifelse(is.na(sym_long$sym_species), 
                               sym_long$sym_RS_scored, 
                               sym_long$sym_species)
sum(!is.na(sym_long$sym_species))
# symmetry available for 761 of 2030 observations, 1269 to score

# exported to csv to score symmetry for remaining taxa
# read in scored symmetry data thus far, and join onto sym_long
symscored <- readxl::read_xlsx("data_input/longevity_symmetry_all_dataentry.xlsx",
                               sheet = 1, na = c("", "NA"), guess_max = 2000)
symscored <- symscored %>%
  dplyr::select(Accepted_name, sym_species) %>%
  dplyr::filter(!is.na(Accepted_name)) %>%
  dplyr::distinct()
sum(is.na(symscored$sym_species))
# 25 taxa left to score symmetry for, all difficult remainders

sym_long <- sym_long %>%
  dplyr::select(-sym_species) %>%
  dplyr::left_join(symscored, by = "Accepted_name")
rm(symscored)

# export final version to csv 
readr::write_csv(sym_long, "data_output/longevity_symmetry_all.csv")
})
