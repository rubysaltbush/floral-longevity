# assemble available symmetry data from multiple sources

### symmetry data ####

sym_data <- cache_csv("data_output/sym_data.csv", function() {

# read in data from different sources
#* published symmetry studies ----
yoderetal <- readr::read_csv("data_input/Yoder_et_al._(2020)_plant_degree-sharing.csv")
jolyschoen <- readr::read_csv("data_input/Joly&Schoen_SupplementaryDataset.csv")

# cut down to variables of interest, species and symmetry
# Yoder et al. data just has symmetry and species
yoderetal <- yoderetal %>%
  dplyr::select(yea_name = plant, symmetry) %>%
  dplyr::filter(!is.na(symmetry)) %>%
  dplyr::distinct()
table(yoderetal$symmetry)
# 2736 separate species, 496 zygomorphic

# Joly and Schoen data itself assembled from various sources and taxonomically matched
jolyschoen <- jolyschoen %>%
  dplyr::select(js_name = Species, symmetry = Symmetry) %>%
  dplyr::filter(!is.na(symmetry)) %>%
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

jolyschoen <- jolyschoen %>%
  dplyr::mutate(source = "jolyandschoen2021") %>%
  dplyr::select(og_species = js_name, symmetry, source)
yoderetal <- yoderetal %>%
  dplyr::mutate(source = "yoderetal2020") %>%
  dplyr::select(og_species = yea_name, symmetry, source)

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

#* RS Australian symmetry data ----
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
aus_sym <- aus_sym %>%
  dplyr::mutate(source = "stephens2021") %>%
  dplyr::select(og_species = auspc_name, symmetry, source)
# how many species overlap from this and other data?
paste(sum(aus_sym$og_species %in% sym_data$og_species), "of", nrow(aus_sym), "names from my symmetry scoring in Joly and Schoen and Yoder et al data" )

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

#* TRY database symmetry data ----
# add in and check out TRY data
try_sym <- rtry::rtry_import("data_input/27350_06062023071220_TRY20230606.txt")
# filter out lat/long etc associated data
try_sym <- try_sym %>%
  dplyr::filter(TraitName == "Flower symmetry type (flower shape)")
table(try_sym$Dataset)
# all TRY symmetry data are for African plant species, from two studies - Renske Onstein's Cape Region and Marco Schmidt's Trait data for African Plants
table(try_sym$OrigValueStr)
# need to convert "yes" and "no" to corresponding symmetry values
try_sym <- try_sym %>%
  dplyr::mutate(symmetry = ifelse(OrigValueStr == "yes", OriglName, OrigValueStr)) %>%
  dplyr::filter(symmetry != "no") %>%
  dplyr::select(og_species = AccSpeciesName, symmetry, source = DatasetID) %>%
  dplyr::mutate(source = paste("2023TRY", source, sep = "_")) %>%
  dplyr::distinct()
table(try_sym$symmetry)
try_sym$symmetry <- gsub("bilateral", "zygomorphic", try_sym$symmetry)
try_sym$symmetry <- gsub("radial", "actinomorphic", try_sym$symmetry)
try_sym$symmetry <- gsub("Zygomorphic", "zygomorphic", try_sym$symmetry)
try_sym$symmetry <- gsub("Actinomorphic", "actinomorphic", try_sym$symmetry)
table(try_sym$symmetry)
# how many species overlap from this and other data?
paste(sum(try_sym$og_species %in% sym_data$og_species), "of", nrow(try_sym), "names from TRY symmetry data in Joly and Schoen, Yoder et al and my Aus scored data" )

# try appending them, see if these 340 agree on symmetry - will be 13660 records if so
sym_data <- sym_data %>%
  dplyr::full_join(try_sym, by = c("og_species", "symmetry"))
# check species that disagree on symmetry
disagreements <- sym_data %>%
  dplyr::group_by(og_species) %>%
  dplyr::filter(n() > 1) # 741 rows woah!!
# will ultimately have to resolve these manually by checking original sources
# and strict definition of what I mean by symmetry
rm(disagreements, try_sym)

# 741 records in sym_data now have multiple records, paste these together 
# so I can easily filter and check them
# first consolidate source columns
sym_data <- sym_data %>%
  dplyr::mutate(source = paste(source.x, source.y, source.x.x, source.y.y)) %>%
  dplyr::select(og_species, symmetry, source)
sym_data$source <- gsub(" NA", "", sym_data$source) # get rid of pasted NA values
sym_data$source <- gsub("NA ", "", sym_data$source) # get rid of pasted NA values

# now summarise down to combine symmetry types so one row per taxon
sym_data_sources <- sym_data %>%
  dplyr::group_by(og_species) %>%
  dplyr::summarise(sources = stringr::str_flatten(source, collapse = " "))
sym_data_sources$sources <- gsub("2023TRY_319 2023TRY_319", "2023TRY_319", sym_data_sources$sources)
sym_data_sources$sources <- gsub("2023TRY_350 2023TRY_350", "2023TRY_350", sym_data_sources$sources)
table(sym_data_sources$sources)
sym_data <- sym_data %>%
  dplyr::group_by(og_species) %>%
  dplyr::summarise(sym_all = stringr::str_flatten(symmetry, collapse = " "))
sym_data$sym_all <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_data$sym_all)
sym_data$sym_all <- gsub("zygomorphic zygomorphic", "zygomorphic", sym_data$sym_all)
sym_data$sym_all <- gsub("zygomorphic actinomorphic", "actinomorphic zygomorphic", sym_data$sym_all)
sym_data$sym_all <- gsub("actinomorphic actinomorphic", "actinomorphic", sym_data$sym_all)
table(sym_data$sym_all)
sym_data <- sym_data %>%
  dplyr::left_join(sym_data_sources, by = "og_species")
rm(sym_data_sources)

# export csv of assembled symmetry data for now
readr::write_csv(sym_data, "data_output/sym_data.csv")
})

