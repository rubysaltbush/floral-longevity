# script to repeat PGLS analyses but with latitude added in as a random factor,
# and therefore individual rather than species mean measures of longevity

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
sym_long_lat %>%
  dplyr::filter(!duplicated(species)) %>%
  dplyr::select(sym_species) %>%
  table()
# 977 actinomorphic to 462 zygomorphic taxa

#* prep and prune trees ----

# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# prune allotb tree to 1420 matched taxa
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
# not removing duplicated species as these are multiple obs at dif latitudes etc
pgls_allotb <- sym_long_lat

# reorder data so species order matches order of tips in tree
pgls_allotb <- as.data.frame(allotb$tip.label) %>%
  dplyr::left_join(pgls_allotb, by = c("allotb$tip.label" = "allotb"))

# taxon_name to row names
#rownames(pgls_allotb) <- pgls_allotb[,1] # gah can't, dupes not allowed
#pgls_allotb[,1] <- NULL
# taxon names for evolutionary correlation calculation
spp <- pgls_allotb[,1]

hist(pgls_allotb$mean_long_days)
hist(log(pgls_allotb$mean_long_days))
# will log transform longevity to meet assumptions of Brownian motion

#** figure ----
# plot individuals longevity by latitude, coloured by symmetry
pgls_allotb %>%
  dplyr::select(Lat, mean_long_days, sym_species) %>% # select columns of interest
  ggplot(aes(x = Lat, y = mean_long_days, colour = sym_species)) +
  geom_point() +
  scale_colour_discrete(type = setNames(my_colours$symmetry, c(1, 0)),
                        labels = c("0" = "actinomorphic", "1" = "zygomorphic"),
                        name = "Symmetry") +
  ggpubr::theme_pubr(legend = "right") +
  xlab("Latitude") +
  ylab("Floral longevity (days)")
ggsave("figures/symmetry_longevity_latitude_allotb.pdf", width = 9, height = 5)

#** ttest & PGLS ----
# t-test of longevity by symmetry
ttests <- list()
ttests$allotbspecies <- t.test(pgls_allotb$mean_long_days[pgls_allotb$sym_species == "1"], 
                               pgls_allotb$mean_long_days[pgls_allotb$sym_species == "0"])
ttests$allotbspecies
# Welch Two Sample t-test
# 
# data:  pgls_allotb$mean_long_days[pgls_allotb$sym_species == "1"] and pgls_allotb$mean_long_days[pgls_allotb$sym_species == "0"]
# t = 3.05, df = 867.04, p-value = 0.002358
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.2578909 1.1889634
# sample estimates:
#   mean of x mean of y 
# 4.647234  3.923807 

# run PGLS without latitude as random factor first, to check symmetry by itself
pgls_models <- list()
pgls_models$allotbsymlong <- nlme::gls(log(mean_long_days) ~ sym_species,
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlong)
# Generalized least squares fit by maximum likelihood
# Model: log(mean_long_days) ~ sym_species 
# Data: pgls_allotb 
# AIC  BIC logLik
# -Inf -Inf    Inf
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value  Std.Error  t-value p-value
# (Intercept)  0.4054225 0.30881388 1.312838  0.1894
# sym_species1 0.1862016 0.08386685 2.220206  0.0265
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.029
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -3.3509211 -0.4682190  0.6389591  1.4358032  3.4597957 
# 
# Residual standard error: 0.865882 
# Degrees of freedom: 1702 total; 1700 residual

# why are AIC, BIC etc -Inf???

# seems to have somewhat lower power from excluded taxa, any effect of latitude on longevity?
pgls_models$allotblonglat <- nlme::gls(log(mean_long_days) ~ abs(Lat),
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotblonglat)
# Generalized least squares fit by maximum likelihood
# Model: log(mean_long_days) ~ abs(Lat) 
# Data: pgls_allotb 
# AIC  BIC logLik
# -Inf -Inf    Inf
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept) -0.004466087 0.30106984 -0.014834  0.9882
# abs(Lat)     0.021984497 0.00198582 11.070732  0.0000
# 
# Correlation: 
#   (Intr)
# abs(Lat) -0.129
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -4.1253542 -0.4649974  0.5067301  1.2321093  3.9617194 
# 
# Residual standard error: 0.8374731 
# Degrees of freedom: 1702 total; 1700 residual

# hmm AIC/BIC/loglik also inf here :/

# yes looks like latitude has stronger effect on longevity than symmetry
# can't do effect of latitude on symmetry bc binary variable
# but can dummy check the reverse??
pgls_models$allotbsymlat <- nlme::gls(abs(Lat) ~ sym_species,
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlat)
# Generalized least squares fit by maximum likelihood
# Model: abs(Lat) ~ sym_species 
# Data: pgls_allotb 
# AIC  BIC logLik
# -Inf -Inf    Inf
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  19.473488  3.647534 5.338810  0.0000
# sym_species1  0.587961  0.990587 0.593548  0.5529
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.029
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -1.8763647 -0.3526659  0.8151594  1.8114754  4.8676047 
# 
# Residual standard error: 10.22731 
# Degrees of freedom: 1702 total; 1700 residual

# latitude does not vary by symmetry, which agrees with the scatter plot

# how do symmetry and latitude interact as fixed effects on longevity
pgls_models$allotbsymlatlong <- nlme::gls(log(mean_long_days) ~ sym_species*abs(Lat),
                                        correlation = ape::corBrownian(phy = allotb, 
                                                                       form = ~spp),
                                        data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlatlong)
# Generalized least squares fit by maximum likelihood
# Model: log(mean_long_days) ~ sym_species * abs(Lat) 
# Data: pgls_allotb 
# AIC  BIC logLik
# -Inf -Inf    Inf
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept)           -0.04820067 0.30176500 -0.159729  0.8731
# sym_species1           0.29032678 0.13163966  2.205466  0.0276
# abs(Lat)               0.02344960 0.00240130  9.765379  0.0000
# sym_species1:abs(Lat) -0.00431903 0.00382919 -1.127923  0.2595
# 
# Correlation: 
#   (Intr) sym_s1 abs(L)
# sym_species1          -0.078              
# abs(Lat)              -0.150  0.437       
# sym_species1:abs(Lat)  0.078 -0.788 -0.563
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -4.1569225 -0.4459726  0.4417050  1.1686265  4.0138449 
# 
# Residual standard error: 0.836035 
# Degrees of freedom: 1702 total; 1698 residual

# looks like the effect of symmetry on floral longevity is independent of the 
# effect of latitude on floral longevity, at least in these data

rm(spp, pgls_allotb, allotb, sym_long_lat, ttests, pgls_models)

# not sure why all models coming out with infinite AIC, BIC and loglik, 
# is it because I am leaving duplicates of individuals in the analysis?
# may have to exclude these or average latitudes and longevity by species?
