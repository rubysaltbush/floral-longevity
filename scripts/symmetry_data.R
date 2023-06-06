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

# TO DO - PROTEUS data? Have checked and not much extra in Schonenberger et al. (2020) data
# add in and check out TRY data (could get via R package??)

### MATCH TAXONOMY ####

