# script to taxonomically match species with longevity and symmetry data
# to closest representatives in ALLOTB and GBOTB trees from Smith and Brown
# (2018)

phylo_names_match <- cache_csv("data_output/data_phylo_names_matched.csv",
                              function() {

# read in short Smith and Brown tree (GBOTB.tre with 79,881 tips)
gbotb <- ape::read.tree("data_input/GBOTB.tre")
# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# check how many species names have direct match in tree
# first filter data down to species names
phylo_names_match <- sym_long %>%
  dplyr::select(species = Accepted_name, og_species_patch, 
                genus = Accepted_genus, family = Accepted_family) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(species))
phylo_names_match$species <- gsub(" ", "_", phylo_names_match$species)
phylo_names_match$og_species_patch <- gsub(" ", "_", phylo_names_match$og_species_patch)

# check how many accepted species directly match
sum(unique(phylo_names_match$species) %in% gbotb$tip.label)
# 808 species with direct matches! 
sum(unique(phylo_names_match$species) %in% allotb$tip.label)
# 1212 species with direct matches!

# and check original names too in case these directly match
sum(unique(phylo_names_match$og_species_patch) %in% gbotb$tip.label)
# 780 with direct matches
sum(unique(phylo_names_match$og_species_patch) %in% allotb$tip.label)
# 1170 with direct matches

# now check which names match and how to improve matching
# first build df of names from trees
tree_names <- tibble(allotb = allotb$tip.label)
gbotb_names <- tibble(allotb = gbotb$tip.label, gbotb = gbotb$tip.label)
tree_names <- dplyr::left_join(tree_names, gbotb_names, by = "allotb")
tree_names$species <- tree_names$allotb
tree_names$match_level <- "direct_accepted"
rm(gbotb_names)

# match names from trees to names in data, accepted names first
phylo_names_match <- dplyr::left_join(phylo_names_match, tree_names, by = "species")
# then filter out these matches
mismatches <- dplyr::filter(phylo_names_match, is.na(allotb))
# then match again, this time based on original species names
tree_names$match_level <- "direct_original"
mismatches <- mismatches %>%
  dplyr::select(species, og_species_patch, genus, family) %>%
  dplyr::left_join(tree_names, by = c("og_species_patch" = "species"))
# matches 59 more, patch these back in
matchedpatch <- dplyr::filter(mismatches, !is.na(allotb))
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(matchedpatch, by = c("species", "og_species_patch",
                                          "genus", "family"))
rm(matchedpatch)

# okay now 192 to go
mismatches <- dplyr::filter(phylo_names_match, is.na(allotb))

# some not matching because of subspecies and var issues
# first create tree_names version with no subsp and var etc
tree_names_nosubsp <- tree_names
tree_names_nosubsp$species_nosubsp <- gsub("_(?:subsp|var|f)\\._.+$", "", tree_names_nosubsp$species)
tree_names_nosubsp <- tree_names_nosubsp %>%
  dplyr::group_by(species_nosubsp) %>%
  dplyr::slice_head(n = 1) %>% # randomly select a variant/subspecies (manually checked, this gets all possible GBOTB matches)
  dplyr::ungroup() %>%
  dplyr::select(-species)
# above is slow!

# now create matching column in mismatches with no var or subsp etc
mismatches$species_nosubsp <- gsub("_(?:subsp|var|f)\\._.+$", "", mismatches$species)
# and get rid of odd hyphens (will fix one mismatch)
mismatches$species_nosubsp <- gsub("-", "", mismatches$species_nosubsp)
# now match based on no subsp in tree name
tree_names_nosubsp$match_level <- "direct_accepted_nosubsp"
mismatches <- mismatches %>%
  dplyr::select(species, og_species_patch, genus, family, species_nosubsp) %>%
  dplyr::left_join(tree_names_nosubsp, by = "species_nosubsp")
# matches 86 more, patch these back in
matchedpatch <- mismatches %>%
  dplyr::filter(!is.na(allotb)) %>%
  dplyr::select(-species_nosubsp)
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(matchedpatch, by = c("species", "og_species_patch",
                                          "genus", "family"))
rm(matchedpatch)

# okay now 106 to go
mismatches <- dplyr::filter(mismatches, is.na(allotb))

# tried stripping names and fuzzy matching but doesn't match many more
# read in some manual updates from checking synonyms on POWO etc
tree_names <- dplyr::select(tree_names, -match_level, -species)
fuzzymatchpatch <- readr::read_csv("data_input/fuzzy_matches.csv") %>%
  dplyr::select(species:family, allotb = manual_match_allotb, match_level) %>%
  dplyr::filter(!is.na(allotb)) %>%
  dplyr::left_join(tree_names, by = "allotb")

# patch in these 14 next, then move to genus level matching
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(fuzzymatchpatch, by = c("species", "og_species_patch",
                                             "genus", "family"))
rm(fuzzymatchpatch)

mismatches <- dplyr::filter(phylo_names_match, is.na(allotb))
# 93 to go, these remaining appear to be missing at species level (even with
# synonyms) from ALLOTB tree so will have to match at accepted genus level
# with preference for genus representatives that have GBOTB match too

# first label taxa with GBOTB presence
tree_names$gbotb_present <- ifelse(is.na(tree_names$gbotb), "0", "1")
# then make genus column
tree_names$genus <- gsub("_.*", "", tree_names$allotb)
# then randomly choose genus level match, preferencing GBOTB match
tree_names_genus <- tree_names %>%
  dplyr::group_by(genus) %>%
  dplyr::slice_max(n = 1, order_by = gbotb_present) %>% # first get all with GBOTB present, or random if no GBOTB match in genus
  dplyr::slice_sample(n = 1) %>% # then randomly choose one in each genus
  dplyr::ungroup()
# match by genus
mismatches <- mismatches %>%
  dplyr::select(species:family) %>%
  dplyr::left_join(tree_names_genus, by = "genus") %>%
  dplyr::select(-gbotb_present) %>%
  dplyr::mutate(match_level = "direct_genus_accepted")
# matches 83 more, patch these back in
matchedpatch <- dplyr::filter(mismatches, !is.na(allotb))
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(matchedpatch, by = c("species", "og_species_patch",
                                          "genus", "family"))
rm(matchedpatch)
# 10 to go, try matching on og genus
mismatches <- dplyr::filter(phylo_names_match, is.na(allotb)) %>%
  dplyr::rename(accepted_genus = genus) %>%
  dplyr::mutate(genus = gsub("_.*", "", og_species_patch))
# match by genus
mismatches <- mismatches %>%
  dplyr::select(species:family, genus) %>%
  dplyr::left_join(tree_names_genus, by = "genus") %>%
  dplyr::select(-gbotb_present) %>%
  dplyr::mutate(match_level = "direct_genus_original")
# this fixes 2 more, patch back in
matchedpatch <- mismatches %>%
  dplyr::filter(!is.na(allotb)) %>%
  dplyr::select(-genus) %>%
  dplyr::rename(genus = accepted_genus)
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(matchedpatch, by = c("species", "og_species_patch",
                                          "genus", "family"))
rm(matchedpatch)

# 8 remaining, guess these will have to be matched to nearest genus in family
mismatches <- mismatches %>%
  dplyr::select(-genus) %>%
  dplyr::rename(genus = accepted_genus) %>%
  dplyr::filter(is.na(allotb))
# write to csv for manual matching
readr::write_csv(mismatches, "data_output/last_few_mismatches.csv")

# have matched last few mismatches to nearest genus representative
# in same family
mismatches <- readr::read_csv("data_input/last_few_mismatches.csv")

# this fixes final 8, patch back in
matchedpatch <- mismatches %>%
  dplyr::filter(!is.na(allotb)) %>%
  dplyr::select(-match_justification)
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(matchedpatch, by = c("species", "og_species_patch",
                                          "genus", "family"))
rm(matchedpatch)
table(phylo_names_match$match_level)
sum(is.na(phylo_names_match$allotb))
rm(mismatches)

# done for ALLOTB!

sum(is.na(phylo_names_match$gbotb))
# now find genus-level matches for 492 without GBOTB match
# first make new column so I can keep track of where GBOTB match has come from
phylo_names_match$match_level_gbotb <- ifelse(is.na(phylo_names_match$gbotb), "NA", phylo_names_match$match_level)
gbotb_missing <- phylo_names_match %>%
  dplyr::filter(is.na(gbotb)) %>%
  dplyr::select(species:family, match_level_gbotb) %>%
  dplyr::left_join(tree_names_genus, by = "genus") %>%
  dplyr::mutate(match_level_gbotb = "direct_genus_accepted")
sum(is.na(gbotb_missing$gbotb))
# matches all but 65, this will do - create patch
gbotb_patch <- gbotb_missing %>%
  dplyr::filter(!is.na(gbotb)) %>%
  dplyr::select(species:family, gbotb, match_level_gbotb)
# and patch back in
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(gbotb_patch, by = c("species", "og_species_patch",
                                          "genus", "family")) %>%
  dplyr::rename(match_level_allotb = match_level)
# matching on og genus only adds 6, think not worth it as a bit confusing taxonomically
rm(gbotb_patch, gbotb_missing, tree_names, tree_names_genus, tree_names_nosubsp)
table(phylo_names_match$match_level_allotb)
table(phylo_names_match$match_level_gbotb)

# 5 taxa matched on original names need to be patched manually as inconsistencies
matchpatch <- readr::read_csv("data_input/phylo_match_patch.csv")
phylo_names_match <- phylo_names_match %>%
  dplyr::rows_update(matchpatch, by = c("species", "og_species_patch",
                                         "genus", "family"))
rm(matchpatch)

# write to csv for now
readr::write_csv(phylo_names_match, "data_output/data_phylo_names_matched.csv")

phylo_names_match
})


