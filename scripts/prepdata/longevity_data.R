#### longevity data ####

# assemble data on floral longevity from multiple sources

longevity_all <- cache_csv("data_output/longevity_data_all.csv", function() {

#* my field data ----
# longevity and symmetry data for 34 native Australian plant species
# collected on fieldwork in Ku-ring-gai National Park and Springwood, NSW
source("scripts/prepdata/field_data.R")

# summarise down to species, mean and SE longevity (days), symmetry, lat + longs
fieldlong <- fieldlong %>%
  dplyr::mutate(SE_long = sd_long/sqrt(n)) %>%
  dplyr::mutate(sym_RS_scored = str_replace(symmetry, "actinomorphic.*", "actinomorphic")) %>%
  dplyr::select(og_species = species, mean_long_days, SE_long, sym_RS_scored, 
                Site = site, Lat = latitude, Lon = longitude) %>%
  dplyr::mutate(sources = "Ruby fieldwork", og_species_patch = og_species)

#* Marcos Méndez's data ----

marcoslong <- cache_csv("data_input/marcoslong_speciespatched.csv", function() {
                          
  # read in .xls sheet with Marcos' highest quality data
  longevhighq <- readxl::read_xls("data_input/Floral_longevity_20230623.xls", sheet = 1)
  longevhighq$sources <- "marcoshighquality"
  # and .xls sheet with lower quality data as well
  longevlowq <- readxl::read_xls("data_input/Floral_longevity_20230623.xls", sheet = 2)
  longevlowq$sources <- "marcoslowquality"
  longevlowq$Exotic <- as.character(longevlowq$Exotic)
  longevlowq$`Altitude (m)` <- as.character(longevlowq$`Altitude (m)`)
  # and .xls sheet with community based data
  longevitycomm <- readxl::read_xls("data_input/Floral_longevity_20230623.xls", sheet = 3)
  longevitycomm$sources <- "marcoscommunity"
  longevitycomm$Exotic <- as.character(longevitycomm$Exotic)
  
  # bind all Marcos' data into one df
  # for now will include all quality levels (0-3) and exotics
  # SHOULD I FILTER OUT ALL GREENHOUSE HABITAT VALUES??? LEAVING IN FOR NOW
  marcoslong <- longevhighq %>%
    dplyr::bind_rows(longevlowq) %>%
    dplyr::bind_rows(longevitycomm) %>%
    dplyr::select(og_species = Species, og_family = Family, 
                  og_longevity = `Floral.longevity (days)`, 
                  SE_long = `SE...8`, Quality, Site, Lat, Lon, 
                  alt_m = `Altitude (m)`, Habitat, reference = Source, 
                  Exotic, Pseudanthium, self_incom = SI, Nectar, notes = ...30,
                  longev_univisit = `Floral.longevity.un/visited (days)`, 
                  long_unv_SE = SE...11, Duration, sources) %>% # keeping everything except flower size for now
    tidyr::fill(og_species, og_family, reference) %>% #fill blank values from row above 
    # Site, Lat, Lon often also blank but hard to autofill as many genuine NAs in these columns
    # will have to fix this manually if important
    dplyr::filter(is.na(Pseudanthium) | Pseudanthium != 1) %>% # exclude Pseudanthium longevity
    dplyr::filter(!is.na(og_longevity)) %>% # exclude taxa missing longevity data (24)
    dplyr::select(-Pseudanthium) %>% # column now blank, don't need
    dplyr::distinct()
  rm(longevhighq, longevlowq, longevitycomm)
  
  # many rows are duplicates from community and other sheets - remove dupes!
  # first paste together og source for each row (high qual, low qual or comm data)
  temp <- marcoslong %>%
    dplyr::group_by(og_species, reference) %>% # group by species and reference
    dplyr::summarise(sources = paste(sources, collapse = " "))
  
  # replace source column with pasted together column to get rid of duplicate rows
  marcoslong <- marcoslong %>%
    dplyr::select(-sources) %>%
    dplyr::left_join(temp, by = c("og_species", "reference")) %>%
    dplyr::distinct()
  rm(temp)
  
  table(marcoslong$og_longevity)  
  
  # now need to process longevity into numeric column
  # need to convert:
  # <12 h -> 0.5
  marcoslong$mean_long_days <- gsub("< 12 h", "0.5", marcoslong$og_longevity)
  # <24 h -> 1
  marcoslong$mean_long_days <- gsub("< 24 h", "1", marcoslong$mean_long_days)
  # 1 (12 h) -> 0.5? or 1? 0.5 for now, though for many studies 1 would be the minimum possible value
  marcoslong$mean_long_days <- gsub("1 \\(12 h\\)", "0.5", marcoslong$mean_long_days)
  # dawn to early evening should be roughly 0.5 by current definition
  # all values starting "dawn" roughly half a day, this will do for now
  marcoslong$mean_long_days <- gsub("[D|d]awn to.*", "0.5", marcoslong$mean_long_days)
  # get rid of notes in longevity column
  marcoslong$mean_long_days <- gsub("Outcrossed.*?: ", "", marcoslong$mean_long_days)
  # and this weird one, take average of midpoints
  marcoslong$mean_long_days <- gsub("\\(9\\)-15-16-\\(>19\\)", "15.5", marcoslong$mean_long_days)
  # finally up to 3 days can be just 3 as no minimum given
  marcoslong$mean_long_days <- gsub("up to 3 d", "3", marcoslong$mean_long_days)
  # and 5 to 6 (10) can be 5.5
  marcoslong$mean_long_days <- gsub("5 to 6 \\(10\\)", "5.5", marcoslong$mean_long_days)
  
  # loop through remaining to fix
  # if 1 or 2, 1 to 2 -> 1.5; take average of min and max e.g. 7 to 10 = 8.5
  # same for hours but /24 to get longevity measure in days
  mean_long_days <- c()
  for(n in 1:nrow(marcoslong)){
    longdat <- marcoslong[[n, "mean_long_days"]]
    if (str_detect(longdat, "^[\\d.]* to [\\d.]* h|^[\\d.]* or [\\d.]* h|^[\\d.]*-[\\d.]* hours")){ # if the longevity is in the form "2 to 3 h" or "1 or 2 h" or "8-9 hours"
      minmax <- as.data.frame(str_match(longdat, "(?<min>[\\d.]*)[ tor-]*(?<max>[\\d.]*)")) # then extract the min and max number
      meanl <- mean(c(as.numeric(minmax$min), as.numeric(minmax$max)))/24 # and return their mean /24 to give longevity in days
    } else if (str_detect(longdat, "^[\\d.]* to [\\d. days]*$|^[\\d.]* or [\\d.]*$|^[\\d.]*-[\\d. d]*$")){ # if the longevity is in the form "2 to 3", "1 or 2" or "1-2" with or without d or days on the end 
      minmax <- as.data.frame(str_match(longdat, "(?<min>[\\d.]*)[ tor-]*(?<max>[\\d.]*)")) # then extract the min and max number
      meanl <- mean(c(as.numeric(minmax$min), as.numeric(minmax$max))) # and return their mean
    } else if (str_detect(longdat, "^\\d*[+]|^\\d* or more|^\\d* to many")){
      meanl <- str_match(longdat, "^\\d*") # with 1+, 2 or more, or 3 to many, return the min number as don't know max
    } else {
      meanl <- longdat
    }
    mean_long_days <- c(mean_long_days, meanl)
  }
  marcoslong$mean_long_days <- mean_long_days
  rm(n, meanl, minmax, longdat, mean_long_days)
  
  marcoslong$mean_long_days <- as.numeric(marcoslong$mean_long_days)
  hist(marcoslong$mean_long_days)
  # woot finally!!!
  
  # now: leaving multiple values per species for now (atomise data)
  # next: patch species names
  # then: match taxonomy between longevity and symmetry to see how many 
  #       taxa I need to score symmetry for
  
  # patch species names in longevity data
  # original names (left) matched by Marcos to World Flora Online (in brackets)
  marcoslong$og_species_patch <- gsub("^.*\\(=|\\)", "", marcoslong$og_species) # remove og names, just keep names manually matched by Marcos
  marcoslong$og_species_patch <- gsub("\\.", " ", marcoslong$og_species_patch) # replace full stops with spaces
  
  # for species patching will need to manually fix a few values, export csv
  readr::write_csv(marcoslong, "data_output/marcoslong_topatchspecies.csv")
  # mostly this manual patching involved converting genus initials to
  # full genus names, plus removing a few notes in species name column
  # and read back in manually patched results
  marcoslong <- readr::read_csv("data_input/marcoslong_speciespatched.csv")
})

#* Song et al (2022) data ----

# read in data downloaded from Supp Mat of Song, B., Sun, L., Barrett, S. C. H., 
# Moles, A. T., Luo, Y.-H., Armbruster, W. S., Gao, Y.-Q., Zhang, S., Zhang, Z.-Q., 
# & Sun, H. (2022). Global analysis of floral longevity reveals latitudinal 
# gradients and biotic and abiotic correlates. New Phytologist, 235. 
# https://doi.org/10.1111/nph.18271

songlong <- readxl::read_xlsx("data_input/Data_for_Dryad.xlsx")

# simplify to columns that will match other longevity data
songlong <- songlong %>%
  dplyr::select(og_species = Species, og_family = Family, Lat = Latitude,
                alt_m = Elevation, mean_long_days = Floral.longevity, 
                pollination = Pollination.mode,
                self_incom = Self.compatibility, reference = referrence) %>%
  dplyr::mutate(sources = "songetal20223", og_species_patch = og_species)

# some species in Song et al data appended "_male" or "_female", transfer
# this info to the notes column
songlong$notes <- NA
for (i in 1:nrow(songlong)) {
  if (str_detect(songlong$og_species_patch[[i]], "[_-]male|[_-]female|[_-]perfect")) { # if the original species name has _male or -female or _perfect in it
    songlong$notes[[i]] <- str_match(songlong$og_species_patch[[i]], "[_-]male|[_-]female|[_-]perfect") # put this in the notes column
    songlong$og_species_patch[[i]] <- gsub("[_-]male|[_-]female|[_-]perfect", "", songlong$og_species_patch[[i]]) # and delete it from the species name column
  } else if (str_detect(songlong$og_species_patch[[i]], " [sx]$")) { # or if it has a " s" (for self pollen) or " x" (for cross pollen)
    songlong$notes[[i]] <- str_match(songlong$og_species_patch[[i]], " [sx]$") # put this in the notes column
    songlong$og_species_patch[[i]] <- gsub(" [sx]$", "", songlong$og_species_patch[[i]]) # and delete it from the species name column
  }
}

# check how many Song et al taxa in Marcos' data
sum(songlong$og_species_patch %in% marcoslong$og_species_patch)
# 260 names directly match, will be good to compare longevity values for these

#** all longev together ----
longevity_all <- marcoslong %>%
  dplyr::bind_rows(songlong) %>%
  dplyr::bind_rows(fieldlong)
rm(marcoslong, fieldlong, songlong)

# patching some erroneous names detected during taxonomic matching
longevity_all$og_species_patch <- gsub("  ", " ", longevity_all$og_species_patch) # few double spaces, convert to single in case that helps at all
longevity_all$og_species_patch[longevity_all$og_species_patch == "Typhonium angustifolium"] <- "Typhonium angustilobum" # checked against original paper, name typo
longevity_all$og_species_patch[longevity_all$og_species_patch == "Walsura trifolia"] <- "Walsura trifoliolata" # checked against distribution, wrong name in og paper

# perhaps here filter out Marcos' wind pollinated taxa??? He says:
# YES, I included some abiotically (wind) pollinated plants. They are good to 
# have as a reference of floral longevity when animals are not involved. You can 
# remove then but remember we have that info in case we feel the need of comparing 
# how long a flower can last when released from the burden of having to attract 
# a pollinator. The families that should be removed from my date base are: 
# Amaranthaceae (I believe), Casuarinaceae, Plantaginaceae (or whatever family 
# they have been synonimized, but retain Veronica species), Poaceae.

# first class Casuarinaceae, Poaceae as abiotically pollinated
longevity_all$pollination[longevity_all$og_family == "Poaceae"] <- "abiotically"
longevity_all$pollination[longevity_all$og_family == "Casuarinaceae"] <- "abiotically"
longevity_all$pollination[longevity_all$og_family == "Cyperaceae"] <- "abiotically" # checked these myself, all wind pollinated
# Sarcopoterium spinosum og paper says wind pollinated
longevity_all$pollination[longevity_all$og_species_patch == "Sarcopoterium spinosum"] <- "abiotically"

# have checked and many Plantaginaceae biotically pollinated, only Plantago lagopus may be ambophilous but will leave for now
# have checked and all Amaranthaceae appear insect pollinated
table(longevity_all$pollination)
# filter out abiotically pollinated taxa (symmetry not relevant for these)
longevity_all <- longevity_all %>%
  dplyr::filter(is.na(pollination)|pollination %in% c("vertebrate", "invertebrate", "NA"))

# export this! then match taxonomy
readr::write_csv(longevity_all, "data_output/longevity_data_all.csv")
})
