# script to repeat PGLS analyses but with latitude added in as a random factor,
# and therefore individual rather than species mean measures of longevity

#### prep trees & data ####

#* prep longev data ----

# reduce longevity data to accepted species, symmetry, latitude
# and remove duplicates
sym_long_lat <- sym_long %>%
  dplyr::select(species = Accepted_name, sym_species, mean_long_days, Lat, og_species_patch) %>%
  dplyr::filter(complete.cases(.)) %>% # filter out 354 obs without latitude (all from Marcos' data, some duplicated but with latitude in Song et al data)
  # actually not happy with this many missing obs, will see if possible to score latitude quickly for these??
  dplyr::distinct() %>% # remove 134 duplicated observations (same species, longevity and latitude)
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
# 890 actinomorphic to 420 zygomorphic taxa

#* prep and prune trees ----

# might remove gbotb???

# read in short Smith and Brown tree (GBOTB.tre with 79,881 tips)
gbotb <- ape::read.tree("data_input/GBOTB.tre")
# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# prune allotb tree to 1295 matched taxa
allotbnames <- sym_long_lat %>%
  dplyr::select(allotb) %>%
  dplyr::distinct()
allotb <- ape::drop.tip(allotb, allotb$tip.label[-match(allotbnames$allotb, allotb$tip.label)])
length(allotb$tip.label)
plot(allotb, type = "fan", show.tip.label = FALSE)
rm(allotbnames)

# prune gbotb tree to 1083 matched taxa
# first filter out NAs
gbotbnames <- sym_long_lat %>%
  dplyr::select(gbotb) %>%
  dplyr::filter(!is.na(gbotb)) %>%
  dplyr::distinct()
gbotb <- ape::drop.tip(gbotb, gbotb$tip.label[-match(gbotbnames$gbotb, gbotb$tip.label)])
length(gbotb$tip.label)
plot(gbotb, type = "fan", show.tip.label = FALSE)
rm(gbotbnames)

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
# t = 3.5043, df = 735.41, p-value = 0.0004855
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.3822576 1.3561641
# sample estimates:
#   mean of x mean of y 
# 4.686797  3.817586 

# run PGLS without latitude as random factor first, to check statistical power
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
# (Intercept)  0.4611582 0.30838511 1.495397  0.1350
# sym_species1 0.1934898 0.08447752 2.290429  0.0221
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.038
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -3.4566537 -0.5390382  0.5778391  1.3861755  3.3217901 
# 
# Residual standard error: 0.8555205 
# Degrees of freedom: 1552 total; 1550 residual

# seems to have lower power from excluded taxa, any effect of latitude on longevity?
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
# (Intercept) 0.1156316 0.30021909  0.385157  0.7002
# abs(Lat)    0.0217685 0.00206119 10.561134  0.0000
# 
# Correlation: 
#   (Intr)
# abs(Lat) -0.118
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -4.3077093 -0.6109077  0.3753145  1.1044921  3.8093974 
# 
# Residual standard error: 0.8277037 
# Degrees of freedom: 1552 total; 1550 residual

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
# (Intercept)  17.04801  3.676386 4.637165  0.0000
# sym_species1  0.48822  1.007091 0.484782  0.6279
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.038
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6437546 -0.1158289  0.9950606  2.0543159  5.1189232 
# 
# Residual standard error: 10.19901 
# Degrees of freedom: 1552 total; 1550 residual

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
# (Intercept)            0.06591911 0.30088045  0.219087  0.8266
# sym_species1           0.29822614 0.13298971  2.242475  0.0251
# abs(Lat)               0.02323757 0.00248334  9.357401  0.0000
# sym_species1:abs(Lat) -0.00435313 0.00396274 -1.098517  0.2721
# 
# Correlation: 
#   (Intr) sym_s1 abs(L)
# sym_species1          -0.082              
# abs(Lat)              -0.139  0.435       
# sym_species1:abs(Lat)  0.076 -0.789 -0.559
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -4.3341063 -0.5838175  0.3145016  1.0374136  3.5248952 
# 
# Residual standard error: 0.826044 
# Degrees of freedom: 1552 total; 1548 residual

# looks like the effect of symmetry on floral longevity is independent of the 
# effect of latitude on floral longevity, at least in these data

rm(spp, pgls_allotb)

