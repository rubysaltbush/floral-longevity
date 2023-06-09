# script to assemble available symmetry data and align taxa to 
# World Checklist of Vascular Plants (WCVP)

### GET DATA ####
# read in data from different sources
yoderetal <- readr::read_csv("data_input/Yoder_et_al._(2020)_plant_degree-sharing.csv")
jolyschoen <- readr::read_csv("data_input/Joly&Schoen_SupplementaryDataset.csv")

# cut down to variables of interest, species and symmetry
# Yoder et al. data just has symmetry and species
yoderetal <- yoderetal %>%
  dplyr::select(yea_name = plant, symmetry) %>%
  dplyr::distinct()
table(yoderetal$symmetry)
# 2736 separate species, 496 zygomorphic

# Joly and Schoen data itself assembled from various sources and taxonomically matched
jolyschoen <- jolyschoen %>%
  dplyr::select(js_name = Species, symmetry = Symmetry) %>%
  dplyr::distinct()
table(jolyschoen$symmetry)
# 3081 species, 1019 zygomorphic

# reformat jolyschoen to match Yoder names and symmetry format
jolyschoen$js_name <- gsub("_", " ", jolyschoen$js_name)
jolyschoen$symmetry <- gsub("Act", "actinomorphic", jolyschoen$symmetry)
jolyschoen$symmetry <- gsub("Zygo", "zygomorphic", jolyschoen$symmetry)

# how many species overlap in these two lists?
paste(sum(jolyschoen$js_name %in% yoderetal$yea_name), "of", nrow(jolyschoen), "names from Joly and Schoen in Yoder et al data" )
paste(sum(yoderetal$yea_name %in% jolyschoen$js_name), "of", nrow(yoderetal), "names from Yoder et al in Joly and Schoen data" )

jolyschoen$og_species <- jolyschoen$js_name
yoderetal$og_species <- yoderetal$yea_name

# try appending them, see if these 299 agree on symmetry
sym_data <- jolyschoen %>%
  dplyr::full_join(yoderetal, by = c("og_species", "symmetry"))
# hmm 38 species have disagreement on symmetry, interesting
# check these
disagreements <- sym_data %>%
  dplyr::group_by(og_species) %>%
  dplyr::filter(n()>1)
# will ultimately have to resolve these manually with checking of original sources
rm(disagreements, jolyschoen, yoderetal)

# add in data for Australian taxa that I scored myself
aus_sym <- readr::read_csv("data_input/aus_symmetry_20210423.csv")
aus_sym <- aus_sym %>%
  dplyr::select(auspc_name = `Taxon name`, symmetry_full = Symmetry) %>%
  dplyr::distinct()
table(aus_sym$symmetry_full)
# categorising disymmetric as zygomorphic as per Joly & Schoen in below
aus_sym$symmetry <- ifelse(str_detect(aus_sym$symmetry_full, ".actinomorphic."), "actinomorphic", "zygomorphic")
table(aus_sym$symmetry)
# need to flatten to one obs per species
aus_sym <- aus_sym %>%
  dplyr::group_by(auspc_name) %>%
  dplyr::summarise(symmetry = stringr::str_flatten(symmetry, collapse = " "))
aus_sym$symmetry <- gsub("zygomorphic actinomorphic", "actinomorphic zygomorphic", aus_sym$symmetry)
table(aus_sym$symmetry)
# then prep to join onto other symmetry data
aus_sym$og_species <- aus_sym$auspc_name
# how many species overlap from this and other data?
paste(sum(aus_sym$auspc_name %in% sym_data$og_species), "of", nrow(aus_sym), "names from my symmetry scoring in Joly and Schoen and Yoder et al data" )

# try appending them, see if these 100 agree on symmetry - will be 8010 records if so
sym_data <- sym_data %>%
  dplyr::full_join(aus_sym, by = c("og_species", "symmetry"))
# 15 species disagree on symmetry, check these
disagreements <- sym_data %>%
  dplyr::group_by(og_species) %>%
  dplyr::filter(n() > 1)
# looks like the disagreement is for Asteraceae and other taxa with pseudanthia
# will ultimately have to resolve these manually by checking original sources
# and strict definition of what I mean by symmetry
rm(disagreements, aus_sym)

# TO DO - PROTEUS data? Have checked and not much extra in Schonenberger et al. (2020) data

# add in and check out TRY data
try_sym <- rtry::rtry_import("data_input/27350_06062023071220_TRY20230606.txt")
# filter out lat/long etc associated data
try_sym <- try_sym %>%
  dplyr::filter(TraitName == "Flower symmetry type (flower shape)")
table(try_sym$Dataset)
# all TRY symmetry data are for African plant species, from two studies - Renske Onstein's Cape Region and Marcos Schmidt's Trait data for African Plants
table(try_sym$OrigValueStr)
# need to convert "yes" and "no" to corresponding symmetry values
test <- try_sym %>%
  dplyr::mutate(symmetry = ifelse(OrigValueStr == "yes", OriglName, OrigValueStr)) %>%
  dplyr::filter(symmetry != "no") %>%
  dplyr::select(try_name = AccSpeciesName, symmetry) %>%
  dplyr::mutate(og_species = try_name) %>%
  dplyr::distinct()
table(test$symmetry)


# Standardization of species names (AccSpeciesName) in TRY version 6: The Plant List 
# has been static since 2013 and is assumed to be outdated. We therefore used the 
# following reference floras for the standardization of species names: World Flora 
# Online (WFO, http://www.worldfloraonline.org/), Leipzig Catalogue of Vascular 
# Plant Names (LCVP, https://idiv-biodiversity.github.io/lcvplants), Tropicos 
# (https://www.tropicos.org), and Index Fungorum (http://www.indexfungorum.org). 
# If a fuzzy match with acceptable levenshtein distance (in general: <2 in the 
# genus, <3 in the epithet) could be established to an accepted name in WFO, this 
# name was used for standardization, else the accepted name in one of the other 
# reference floras was used, if available.

### MATCH TAXONOMY ####

