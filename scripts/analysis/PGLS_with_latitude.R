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
  dplyr::filter(!duplicated(species)) %>%
  dplyr::select(sym_species) %>%
  table()
# 972 actinomorphic to 451 zygomorphic taxa

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

#** ttest & PGLS ----
# t-test of longevity by symmetry
ttests <- list()
ttests$allotblong <- t.test(pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "1"], 
                               pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "0"])
ttests$allotblong
# Welch Two Sample t-test
# 
# data:  pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "1"] and pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "0"]
# t = 4.3647, df = 710.24, p-value = 1.462e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.6167554 1.6252529
# sample estimates:
#   mean of x mean of y 
# 4.682592  3.561587 

# t-test of latitude by symmetry
ttests$allotblat <- t.test(pgls_allotb$spmeanlat[pgls_allotb$sym_species == "1"], 
                            pgls_allotb$spmeanlat[pgls_allotb$sym_species == "0"])
ttests$allotblat
# Welch Two Sample t-test
# 
# data:  pgls_allotb$spmeanlat[pgls_allotb$sym_species == "1"] and pgls_allotb$spmeanlat[pgls_allotb$sym_species == "0"]
# t = 1.1786, df = 921.11, p-value = 0.2389
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.5881533  2.3567204
# sample estimates:
#   mean of x mean of y 
# 26.37643  25.49215 

# run PGLS without latitude first, to check symmetry by itself
pgls_models <- list()
pgls_models$allotbsymlong <- nlme::gls(log(spmean_long_days) ~ 
                                         scale(as.numeric(sym_species)),
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlong)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ scale(as.numeric(sym_species)) 
# Data: pgls_allotb 
# AIC      BIC   logLik
# 4524.62 4540.402 -2259.31
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)                    0.9563604 0.7457652 1.282388  0.1999
# scale(as.numeric(sym_species)) 0.1426238 0.0413410 3.449937  0.0006
# 
# Correlation: 
#   (Intr)
# scale(as.numeric(sym_species)) 0.017 
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.10113651 -0.28198308 -0.01443645  0.23146808  0.81545011 
# 
# Residual standard error: 3.047146 
# Degrees of freedom: 1423 total; 1421 residual

# effect of (absolute) latitude on longevity?
pgls_models$allotblonglat <- nlme::gls(log(spmean_long_days) ~ scale(spmeanlatabs),
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotblonglat)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ scale(spmeanlatabs) 
# Data: pgls_allotb 
# AIC      BIC    logLik
# 4496.767 4512.548 -2245.383
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)          0.8665400 0.7384333  1.173484  0.2408
# scale(spmeanlatabs) -0.2110542 0.0332768 -6.342376  0.0000
# 
# Correlation: 
#   (Intr)
# scale(spmeanlatabs) 0.01  
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -1.0204262 -0.3280829  0.0318071  0.3029884  0.8935845 
# 
# Residual standard error: 3.017469 
# Degrees of freedom: 1423 total; 1421 residual


# yes looks like latitude has an effect on longevity
# can't do effect of latitude on symmetry bc binary variable (would need phylo logistic regression)
# but can dummy check the reverse??
pgls_models$allotbsymlat <- nlme::gls(spmeanlatabs ~ sym_species,
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlat)
# Generalized least squares fit by maximum likelihood
# Model: spmeanlatabs ~ sym_species 
# Data: pgls_allotb 
# AIC      BIC    logLik
# 11242 11257.78 -5617.998
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)  23.076481  7.901225  2.920620  0.0035
# sym_species1 -1.478082  0.940960 -1.570824  0.1164
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.021
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -0.7060763 -0.2233425  0.1016666  0.4622944  1.4485092 
# 
# Residual standard error: 32.28142 
# Degrees of freedom: 1423 total; 1421 residual

# latitude does not vary by symmetry, which agrees with the scatter plot

# how do symmetry and latitude interact as fixed effects on longevity
pgls_models$allotbsymlatlong <- nlme::gls(log(spmean_long_days) ~ 
                                            scale(as.numeric(sym_species))*scale(spmeanlatabs),
                                        correlation = ape::corBrownian(phy = allotb, 
                                                                       form = ~spp),
                                        data = pgls_allotb, method = "ML")
summary(pgls_models$allotbsymlatlong)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ scale(as.numeric(sym_species)) * scale(spmeanlatabs) 
# Data: pgls_allotb 
# AIC      BIC    logLik
# 4458.418 4484.721 -2224.209
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)                                         0.9536605 0.7281932  1.309626  0.1905
# scale(as.numeric(sym_species))                      0.1691634 0.0409237  4.133630  0.0000
# scale(spmeanlatabs)                                -0.1239595 0.0359223 -3.450767  0.0006
# scale(as.numeric(sym_species)):scale(spmeanlatabs) -0.1628957 0.0287153 -5.672776  0.0000
# 
# Correlation: 
#   (Intr) sc(.(_)) scl(s)
# scale(as.numeric(sym_species))                      0.019                
# scale(spmeanlatabs)                                 0.014  0.102         
# scale(as.numeric(sym_species)):scale(spmeanlatabs) -0.011 -0.160   -0.405
# 
# Standardized residuals:
#   Min           Q1          Med           Q3          Max 
# -1.115757251 -0.285272896  0.003389668  0.264009912  0.857956007 
# 
# Residual standard error: 2.972902 
# Degrees of freedom: 1423 total; 1419 residual

# with all predictors scaled (and therefore numeric) all come out
# as significant, comparing effect sizes symmetry has strongest (+ve) effect,
# closely followed by interaction (-ve) and then latitude (-ve)

rm(spp, pgls_allotb, allotb, sym_long_lat, ttests, pgls_models)

#### BAYESIAN MIXED MODEL APPROACH ####

# going to try using exact same model approach as Song et al. (2022) paper
# to be sure results are comparable, given odd difference in effect
# of latitude on longevity (positive in Song et al paper, negative here)

library(brms) # Bayesian mixed models package

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
# sym_species
# 0   1 
# 980 462 

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

# reorder data so species order matches order of tips in tree
brm_allotb <- as.data.frame(allotb$tip.label) %>%
  dplyr::left_join(sym_long_lat, by = c("allotb$tip.label" = "allotb")) %>%
  dplyr::rename(allotb = `allotb$tip.label`)

# for Bayesian modelling, want a few extra columns of variables/cofactors
brm_allotb <- brm_allotb %>%
  dplyr::mutate(abs_lat = abs(Lat))
speciesmeans <- brm_allotb %>%
  dplyr::group_by(allotb) %>%
  dplyr::summarise(spmeanlongdays = mean(mean_long_days),
                   spmeanabslat = mean(abs_lat))
brm_allotb <- brm_allotb %>%
  dplyr::left_join(speciesmeans, by = "allotb")
rm(speciesmeans)

# set up model ----

# get variance-covariance matrix for phylogeny
allotb_vcv <- ape::vcv.phylo(phy = allotb, model = "Brownian")

# longevity by latitude ----

# try running Bayesian mixed model as per Song et al. (2022), with 
# latitude as fixed effect and phylo and species as random effects
brm_model_longlat <- brms::brm(log(mean_long_days) ~ 
                                 scale(spmeanabslat) +
                                 (1|gr(allotb, cov = allotb_vcv)) +
                                 (1|species),
          data = brm_allotb,
          family = gaussian(),
          data2 = list(allotb_vcv = allotb_vcv),
          iter = 4000,
          cores = 2)

summary(brm_model_longlat)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: log(mean_long_days) ~ scale(spmeanabslat) + (1 | gr(allotb, cov = allotb_vcv)) + (1 | species) 
# Data: brm_allotb (Number of observations: 1705) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~allotb (Number of levels: 1423) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.08      0.00     0.07     0.09 1.00     1125     2236
# 
# ~species (Number of levels: 1442) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.42      0.03     0.37     0.48 1.01      771     1597
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.99      0.24     0.52     1.45 1.00     2056     3485
# scalespmeanabslat     0.32      0.03     0.27     0.38 1.00     3939     5233
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.37      0.02     0.34     0.41 1.01      909     1610
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# standardised effect size of 0.32 for latitude, 95% CI of 0.27-0.38

plot(brm_model_longlat, nvariables = 2, ask = FALSE)
# at this stage don't actually know what to interpret from these plots -
# estimates all have normal distribution, that seems good?

plot(brms::conditional_effects(brm_model_longlat), points = TRUE)
# need to develop export-worthy version of above plot

# compute phylogenetic signal, following method in vignette
# which gives estimate of intra-class correlation
hyp <- "sd_allotb__Intercept^2 / (sd_allotb__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(brm_model_longlat, hyp, class = NULL))
# Hypothesis Tests for class :
#                   Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (sd_allotb__Inter... = 0     0.05      0.01     0.04     0.06         NA        NA    *
#      ---
#      'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#    '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#    for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#    Posterior probabilities of point hypotheses assume equal prior probabilities.
plot(hyp)
rm(hyp)

# but! need to add intraspecific variability of latitude in
# make new column of latitude variability
brm_allotb$within_spec_lat <- brm_allotb$abs_lat - brm_allotb$spmeanabslat

# and then fit model again using within_spec_lat as an additional predictor.

brm_model_longlat2 <- update(
  brm_model_longlat, formula = ~ . + within_spec_lat,
  newdata = brm_allotb, cores = 2,
  iter = 4000
)

summary(brm_model_longlat2)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: log(mean_long_days) ~ scale(spmeanabslat) + (1 | gr(allotb, cov = allotb_vcv)) + (1 | species) + within_spec_lat 
# Data: brm_allotb (Number of observations: 1705) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~allotb (Number of levels: 1423) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.08      0.00     0.07     0.09 1.00     1018     2311
# 
# ~species (Number of levels: 1442) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.42      0.03     0.37     0.48 1.01      645     1648
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             0.99      0.23     0.53     1.45 1.00     2147     3304
# scalespmeanabslat     0.32      0.03     0.27     0.38 1.00     3396     4728
# within_spec_lat      -0.02      0.01    -0.04     0.00 1.00    13306     6582
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.37      0.02     0.34     0.41 1.00      776     1976
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# looks like within species latitude variation has no impact on longevity
# given within spec estimate 95% CI includes 0 (no slope)

# check no change to phylogenetic signal
hyp <- "sd_allotb__Intercept^2 / (sd_allotb__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(brm_model_longlat2, hyp, class = NULL))
# Hypothesis Tests for class :
#   Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (sd_allotb__Inter... = 0     0.05      0.01     0.04     0.06         NA        NA    *
#      ---
#      'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#    '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#    for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#    Posterior probabilities of point hypotheses assume equal prior probabilities.
rm(hyp)
# yep, phylogenetic signal unchanged

# longevity by symmetry ----

# symmetry has to be numeric to scale
brm_allotb$sym_species <- as.numeric(brm_allotb$sym_species)

# run model with phylogenetic covariance and species as random factors
brm_model_longsym <- brms::brm(log(mean_long_days) ~ 
                                 scale(sym_species) +
                                 (1|gr(allotb, cov = allotb_vcv)) +
                                 (1|species),
                               data = brm_allotb,
                               family = gaussian(),
                               data2 = list(allotb_vcv = allotb_vcv),
                               iter = 4000,
                               cores = 2)

summary(brm_model_longsym)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: log(mean_long_days) ~ scale(sym_species) + (1 | gr(allotb, cov = allotb_vcv)) + (1 | species) 
# Data: brm_allotb (Number of observations: 1705) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~allotb (Number of levels: 1423) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.10      0.00     0.09     0.11 1.00     1071     2649
# 
# ~species (Number of levels: 1442) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.39      0.03     0.33     0.45 1.01      621      746
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            0.92      0.29     0.37     1.49 1.00     1825     2469
# scalesym_species     0.09      0.04     0.01     0.17 1.00     5188     5372
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.38      0.02     0.35     0.41 1.01      644     1493
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# symmetry still has effect on longevity, standardised effect estimate of 0.09
# so weaker than latitude but 95% CI of 0.01-0.17 so non-zero positive effect


plot(brm_model_longsym, nvariables = 2, ask = FALSE)
# at this stage don't actually know what to interpret from these plots -
# estimates all have normal distribution, that seems good?

plot(brms::conditional_effects(brm_model_longsym), points = TRUE)
# TO DO - make above plot export worthy

# compute phylogenetic signal, following method in vignette
# which gives estimate of intra-class correlation
hyp <- "sd_allotb__Intercept^2 / (sd_allotb__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(brm_model_longsym, hyp, class = NULL))
# Hypothesis Tests for class :
#   Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (sd_allotb__Inter... = 0     0.06      0.01     0.05     0.08         NA        NA    *
#      ---
#      'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#    '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#    for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#    Posterior probabilities of point hypotheses assume equal prior probabilities.
plot(hyp)
rm(hyp)


#  symmetry and latitude combined ----

# try adding symmetry and latitude together!
brm_model_longsymlat <- brms::brm(log(mean_long_days) ~ 
                                 scale(sym_species) +
                                 scale(spmeanabslat) +
                                 (1|gr(allotb, cov = allotb_vcv)) +
                                 (1|species),
                               data = brm_allotb,
                               family = gaussian(),
                               data2 = list(allotb_vcv = allotb_vcv),
                               iter = 4000,
                               cores = 2)

summary(brm_model_longsymlat)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: log(mean_long_days) ~ scale(sym_species) + scale(spmeanabslat) + (1 | gr(allotb, cov = allotb_vcv)) + (1 | species) 
# Data: brm_allotb (Number of observations: 1705) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~allotb (Number of levels: 1423) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.08      0.00     0.07     0.09 1.00      950     2091
# 
# ~species (Number of levels: 1442) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.42      0.03     0.37     0.48 1.00      702     1239
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept             1.03      0.24     0.57     1.51 1.00     1691     2905
# scalesym_species      0.09      0.04     0.01     0.16 1.00     3675     5178
# scalespmeanabslat     0.33      0.03     0.27     0.38 1.00     3537     4897
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.37      0.02     0.34     0.41 1.00      897     2209
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

brm_model_longsymlatint <- brms::brm(log(mean_long_days) ~ 
                                       scale(sym_species) +
                                       scale(spmeanabslat) +
                                      scale(sym_species):scale(spmeanabslat) +
                                    (1|gr(allotb, cov = allotb_vcv)) +
                                    (1|species),
                                  data = brm_allotb,
                                  family = gaussian(),
                                  data2 = list(allotb_vcv = allotb_vcv),
                                  iter = 4000,
                                  cores = 2)

summary(brm_model_longsymlatint)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: log(mean_long_days) ~ scale(sym_species) + scale(spmeanabslat) + scale(sym_species):scale(spmeanabslat) + (1 | gr(allotb, cov = allotb_vcv)) + (1 | species) 
# Data: brm_allotb (Number of observations: 1705) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~allotb (Number of levels: 1423) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.08      0.00     0.07     0.09 1.00     1299     2565
# 
# ~species (Number of levels: 1442) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.42      0.03     0.37     0.48 1.01      680     1404
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                              1.02      0.24     0.53     1.48 1.00     1669     2881
# scalesym_species                       0.09      0.04     0.01     0.16 1.00     4383     5435
# scalespmeanabslat                      0.33      0.03     0.27     0.38 1.00     3560     5103
# scalesym_species:scalespmeanabslat     0.00      0.03    -0.05     0.05 1.00     4586     5586
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.37      0.02     0.34     0.41 1.00      731     1725
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# question to answer for myself - is it possible for estimates in this analysis to be negative? should be?

