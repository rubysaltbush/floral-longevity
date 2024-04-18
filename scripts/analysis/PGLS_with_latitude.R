# script to repeat PGLS analyses but with species mean absolute latitude added in
# to check the effect of latitudinal variation in longevity on relationship
# with floral symmetry

#### prep data and trees ####

#* prep longev data ----

# reduce longevity data to accepted species, symmetry, latitude
# and remove duplicates
sym_long_lat <- sym_long %>%
  dplyr::select(species = Accepted_name, sym_species, mean_long_days, Lat, og_species_patch) %>%
  dplyr::filter(complete.cases(.)) %>% # filter out ~125 obs without latitude 
  dplyr::distinct() %>% # remove ~121 duplicated observations (same species, longevity and latitude)
  as.data.frame()

# add underscore to match tip labels
sym_long_lat$species <- gsub(" ", "_", sym_long_lat$species)
sym_long_lat$og_species_patch <- gsub(" ", "_", sym_long_lat$og_species_patch)

# match allotb and gbotb names to this data
# read in taxonomic name matching key
source("scripts/prepdata/phylo_names_match.R")
# phylo_names_match has 3 duplicated species because different og names give
# different match levels, will remove worse matching ones
sym_long_lat <- sym_long_lat %>%
  dplyr::left_join(phylo_names_match, by = c("species", "og_species_patch"))
rm(phylo_names_match)

# define heirarchy of taxonomic matches so can sample best matches first
sym_long_lat$allotb_matchrank <- ifelse(sym_long_lat$match_level_allotb %in% c("direct_accepted", 
                                                                               "direct_accepted_nosubsp",
                                                                               "manual_misspelling"), 
                                        "1", 
                                        ifelse(sym_long_lat$match_level_allotb %in% c("direct_original",
                                                                                      "manual_synonym"), 
                                               "2",
                                               ifelse(sym_long_lat$match_level_allotb == "direct_genus_accepted", 
                                                      "3",
                                                      ifelse(sym_long_lat$match_level_allotb %in% c("direct_genus_original",
                                                                                                    "closest_genus",
                                                                                                    "genus_synonym"), 
                                                             "4", "5"))))
table(sym_long_lat$allotb_matchrank)

# remove 15 duplicate observations, choosing best allotb matchrank
sym_long_lat <- sym_long_lat %>%
  dplyr::group_by(species, sym_species, mean_long_days, Lat, genus, family, 
                  allotb, gbotb) %>% # first group by all columns that should be identical if obs identical
  dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # choose best possible taxonomic match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::select(-og_species_patch) %>% # don't need this column anymore
  dplyr::ungroup()

# redefine symmetry as 0 and 1
sym_long_lat$sym_species <- gsub("zygomorphic", "1", sym_long_lat$sym_species)
sym_long_lat$sym_species <- gsub("actinomorphic", "0", sym_long_lat$sym_species)

# PGLS can't handle multiple observations per species, take species mean
# of longevity and of (absolute) latitude per ACCEPTED species
sym_long_lat <- sym_long_lat %>%
  dplyr::group_by(species, sym_species, allotb, allotb_matchrank) %>%
  dplyr::summarise(spmean_long_days = mean(mean_long_days), 
                   spmeanlatabs = mean(abs(Lat))) %>%
  dplyr::ungroup() %>%
  # and then filter out multiple observations assigned to the same ALLOTB species
  dplyr::group_by(allotb) %>% # group by allotb species
  dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # choose best possible taxonomic match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::ungroup()

sym_long_lat %>%
  dplyr::select(sym_species) %>%
  table()
# 971 actinomorphic to 452 zygomorphic taxa

#* prep and prune trees ----

# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# prune allotb tree to 1423 matched taxa
allotbnames <- sym_long_lat %>%
  dplyr::select(allotb) %>%
  dplyr::distinct()
allotb <- ape::drop.tip(allotb, allotb$tip.label[-match(allotbnames$allotb, allotb$tip.label)])
length(allotb$tip.label)
plot(allotb, type = "fan", show.tip.label = FALSE)
rm(allotbnames)

#### PGLS analyses ####

#* ALLOTB PGLS ----

#** prep data ----

pgls_allotb <- sym_long_lat

# reorder data so species order matches order of tips in tree
pgls_allotb <- as.data.frame(allotb$tip.label) %>%
  dplyr::left_join(pgls_allotb, by = c("allotb$tip.label" = "allotb"))

# taxon names for evolutionary correlation calculation
spp <- pgls_allotb[,1]
# taxon_name to row names
rownames(pgls_allotb) <- pgls_allotb[,1]
pgls_allotb[,1] <- NULL

hist(pgls_allotb$spmean_long_days)
hist(log(pgls_allotb$spmean_long_days))
# will log transform longevity to meet assumptions of Brownian motion

#** figure ----
# plot species mean longevity by species mean latitude, coloured by symmetry
pgls_allotb %>%
  dplyr::select(spmeanlatabs, spmean_long_days, sym_species) %>% # select columns of interest
  ggplot(aes(x = spmeanlatabs, y = log(spmean_long_days), colour = sym_species)) +
  geom_point() +
  scale_colour_discrete(type = setNames(my_colours$symmetry, c(1, 0)),
                        labels = c("0" = "actinomorphic", "1" = "zygomorphic"),
                        name = "Symmetry") +
  ggpubr::theme_pubr(legend = "right") +
  xlab("Latitude (absolute)") +
  ylab("Floral longevity (log # days)")
ggsave("figures/symmetry_loglongevity_abs_latitude_allotb.pdf", width = 9, height = 5)

#** PGLS ----

# list to store model output
pgls_models <- list()

# what is the effect of symmetry and latitude on longevity ?
pgls_models$allotbsymlatlong <- nlme::gls(log(spmean_long_days) ~ 
                                            scale(as.numeric(sym_species)) 
                                          + scale(spmeanlatabs),
                                          correlation = ape::corBrownian(phy = allotb, 
                                                                         form = ~spp),
                                          data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlatlong)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ scale(as.numeric(sym_species)) + scale(spmeanlatabs) 
# Data: pgls_allotb 
# AIC      BIC    logLik
# 4491.906 4512.948 -2241.953
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)                     0.9090730 0.7370253  1.233435  0.2176
# scale(as.numeric(sym_species))  0.1367322 0.0408635  3.346076  0.0008
# scale(spmeanlatabs)            -0.2087726 0.0332832 -6.272614  0.0000
# 
# Correlation: 
#   (Intr) s(.(_)
#            scale(as.numeric(sym_species)) 0.017        
#            scale(spmeanlatabs)            0.011  0.044 
#            
#            Standardized residuals:
#              Min           Q1          Med           Q3          Max 
#            -1.007123509 -0.319472645  0.009280592  0.279001215  0.838333798 
#            
#            Residual standard error: 3.010204 
#            Degrees of freedom: 1423 total; 1420 residual

# with all predictors scaled (and therefore numeric) all come out
# as significant, comparing effect sizes latitude has strongest (-ve) effect,
# followed by symmetry (+ve)

rm(spp, pgls_allotb, allotb, sym_long_lat, pgls_models)

# tried Bayesian models as per Song et al. (2022) paper, using brms package
# and running full interaction model symmetry still has significant +ve effect,
# this different analysis method gives a +ve effect of latitude on 
# longevity as reported in the Song et al and no interaction between symmetry
# and longevity. Using a different analysis method and getting a different
# result for latitude is slightly concerning but given symmetry comes out as
# having a significant +ve effect in both analyses I will focus on this for
# now, as I was not aiming to replicate Song et al's latitude result.
# I suspect the different effects of latitude on longevity with different 
# analysis methods are due to the different ways phylogeny is handled by 
# these analyses.