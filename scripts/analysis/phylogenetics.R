# this script adapts Smith and Brown (2018) evolutionary trees to sample taxa,
# runs phylogenetic generalised least squares regressions (PGLS) to see if 
# flowers have evolved to last longer or shorter in taxa with zygomorphic 
# flowers, and test the phylogenetic signal and evolutionary rate of floral 
# symmetry and floral longevity evolution

#### prep trees & data ####

# read in short Smith and Brown tree (GBOTB.tre with 79,881 tips)
gbotb <- ape::read.tree("data_input/GBOTB.tre")
# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# read in taxonomic name matching key
source("scripts/prepdata/phylo_names_match.R")

#* prune trees ----

# prune allotb tree to 1433 matched taxa
allotb <- ape::drop.tip(allotb, allotb$tip.label[-match(phylo_names_match$allotb, allotb$tip.label)])
length(allotb$tip.label)
plot(allotb, type = "fan", show.tip.label = FALSE)

# prune gbotb tree to 1187 matched taxa (many duplicated as genus-only matches)
# first filter out NAs
gbotbnames <- dplyr::filter(phylo_names_match, !is.na(phylo_names_match$gbotb))
gbotbnames <- gbotbnames$gbotb
gbotb <- ape::drop.tip(gbotb, gbotb$tip.label[-match(gbotbnames, gbotb$tip.label)])
rm(gbotbnames)
length(gbotb$tip.label)
plot(gbotb, type = "fan", show.tip.label = FALSE)

#* prep longev data ----

# summarise longevity per ACCEPTED species
# take mean per species for longevity
spmean_long <- sym_long %>%
  dplyr::group_by(species = Accepted_name, sym_species) %>%
  dplyr::summarise(spmean_long_days = mean(mean_long_days)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(species)) %>%
  as.data.frame()
# add underscore to match tip labels
spmean_long$species <- gsub(" ", "_", spmean_long$species)

# match allotb and gbotb names to this data
phylo_names_match <- phylo_names_match %>%
  dplyr::select(species, genus:match_level_gbotb) %>%
  dplyr::distinct()
# above df has 3 rows more than spmean_long because of different og name matches giving
# different match levels, will leave these duplicates for now and should be
# fixed when I choose one species per genus
spmean_long <- spmean_long %>%
  dplyr::left_join(phylo_names_match, by = "species")
rm(phylo_names_match)

# redefine symmetry as 0 and 1
spmean_long$sym_species <- gsub("zygomorphic", "1", spmean_long$sym_species)
spmean_long$sym_species <- gsub("actinomorphic", "0", spmean_long$sym_species)
table(spmean_long$sym_species)
# 986 actinomorphic to 469 zygomorphic taxa (will change with subsampling)

# define heirarchy of taxonomic matches so can sample best matches first
spmean_long$allotb_matchrank <- ifelse(spmean_long$match_level_allotb %in% c("direct_accepted", 
                                                                             "direct_accepted_nosubsp",
                                                                             "manual_misspelling"), 
                                       "1", 
                                       ifelse(spmean_long$match_level_allotb %in% c("direct_original",
                                                                                    "manual_synonym"), 
                                              "2",
                                              ifelse(spmean_long$match_level_allotb == "direct_genus_accepted", 
                                                     "3",
                                                     ifelse(spmean_long$match_level_allotb %in% c("direct_genus_original",
                                                                                                  "closest_genus",
                                                                                                  "genus_synonym"), 
                                                            "4", "5"))))
table(spmean_long$allotb_matchrank)

# create column of allotb genus for subsampling
spmean_long$allotb_genus <- gsub("_.*", "", spmean_long$allotb)

# add orders and clades for labelling in phylogeny and subsampling in clades
familyorderclade <- readr::read_csv("data_input/family_order_clade.csv")
spmean_long <- spmean_long %>%
  dplyr::left_join(familyorderclade, by = "family")
rm(familyorderclade)

#### SPECIES MEAN LONGEVITY DATA DESCRIPTION ####

# what are the longest lived flowers ON AVERAGE?
spmean_long %>%
  dplyr::arrange(dplyr::desc(spmean_long_days)) %>%
  dplyr::select(species, family, sym_species, spmean_long_days) %>%
  dplyr::distinct() %>%
  head()
#                 species           family sym_species spmean_long_days
# 1  Telipogon_peruvianus      Orchidaceae           1            33.00
# 2     Siparuna_muricata     Siparunaceae           0            28.75
# 3     Encyclia_mapuerae      Orchidaceae           1            27.74
# 4 Pinguicula_longifolia Lentibulariaceae           1            27.10
# 5       Pieris_japonica        Ericaceae           0            23.35
# 6   Corybas_cheesemanii      Orchidaceae           1            23.14

# what are the shortest lived flowers ON AVERAGE?
spmean_long %>%
  dplyr::arrange(spmean_long_days) %>%
  dplyr::select(species, family, sym_species, spmean_long_days) %>%
  dplyr::distinct() %>%
  head()
#                   species        family sym_species spmean_long_days
# 1     Heteracia_szovitsii    Asteraceae           0       0.08240741
# 2    Moraea_pseudospicata     Iridaceae           0       0.10416667
# 3     Talinum_paniculatum    Talinaceae           0       0.16666667
# 4     Agalinis_auriculata Orobanchaceae           1       0.18750000
# 5 Tuberaria_bupleurifolia     Cistaceae           0       0.19791667
# 6       Tuberaria_guttata     Cistaceae           0       0.19791667

# what are the longest lived FAMILIES on average?
spmean_long %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(meanlong = mean(spmean_long_days), n = n()) %>%
  dplyr::arrange(dplyr::desc(meanlong)) %>%
  dplyr::filter(n >= 10) %>%
  head()

# family           meanlong     n
# <chr>               <dbl> <int>
# 1 Lentibulariaceae    11.1     10
# 2 Orchidaceae         11.0     56
# 3 Saxifragaceae       10.6     12
# 4 Ranunculaceae        8.50    45
# 5 Violaceae            7.55    10
# 6 Liliaceae            7.53    12

# checked and all Lentibulariaceae in one genus (Pinguicola) so not sure if
# representative of family as a whole

# what are the shortest lived FAMILIES on average?
spmean_long %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(meanlong = mean(spmean_long_days), n = n()) %>%
  dplyr::arrange(meanlong) %>%
  dplyr::filter(n >= 10) %>%
  head()

# family          meanlong     n
# <chr>              <dbl> <int>
# 1 Convolvulaceae     0.700    14
# 2 Cistaceae          0.822    12
# 3 Cactaceae          0.913    44
# 4 Bromeliaceae       1.25     32
# 5 Melastomataceae    1.37     30
# 6 Euphorbiaceae      1.57     10

# what are the longest lived ORDERS on average?
spmean_long %>%
  dplyr::group_by(order) %>%
  dplyr::summarise(meanlong = mean(spmean_long_days), n = n()) %>%
  dplyr::arrange(dplyr::desc(meanlong)) %>%
  dplyr::filter(n >= 10) %>%
  head()
# order        meanlong     n
# <chr>           <dbl> <int>
# 1 Liliales         9.34    27
# 2 Saxifragales     9.31    23
# 3 Ranunculales     7.69    63
# 4 Asparagales      7.10   121
# 5 Laurales         7.06    11
# 6 Ericales         5.20    91

# what are the shortest lived ORDERS on average?
spmean_long %>%
  dplyr::group_by(order) %>%
  dplyr::summarise(meanlong = mean(spmean_long_days), n = n()) %>%
  dplyr::arrange(meanlong) %>%
  dplyr::filter(n >= 10) %>%
  head()
# order        meanlong     n
# <chr>           <dbl> <int>
# 1 Poales           1.09    39
# 2 Boraginales      1.91    27
# 3 Solanales        1.99    43
# 4 Arecales         2.09    10
# 5 Fabales          2.26   117
# 6 Zingiberales     2.35    42

#### PGLS analyses ####

#* ALLOTB PGLS ----

#** prep data ----
# remove allotb duplicate matches from data, randomly choose best matches
pgls_allotb <- spmean_long

pgls_allotb <- pgls_allotb %>%
  dplyr::group_by(allotb) %>% # group by allotb species
  dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # choose best possible taxonomic match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::ungroup()

# reorder data so species order matches order of tips in tree
pgls_allotb <- as.data.frame(allotb$tip.label) %>%
  dplyr::left_join(pgls_allotb, by = c("allotb$tip.label" = "allotb"))

# taxon_name to row names
rownames(pgls_allotb) <- pgls_allotb[,1]
pgls_allotb[,1] <- NULL
# taxon names for evolutionary correlation calculation
spp <- rownames(pgls_allotb)

hist(pgls_allotb$spmean_long_days)
hist(log(pgls_allotb$spmean_long_days))
# will log transform longevity to meet assumptions of Brownian motion

#** figure ----
# boxplot of longevity by symmetry, allotb species
ggplot(data = pgls_allotb, aes(x = spmean_long_days, y = sym_species, fill = sym_species)) +
  geom_violin() +
  geom_boxplot(width = 0.15, colour = "grey", alpha = 0.3) +
  scale_fill_discrete(type = setNames(my_colours$symmetry, c(1, 0))) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("") +
  scale_y_discrete(labels = c("0" = "actinomorphic", "1" = "zygomorphic"))
ggsave("figures/symmetry_longevity_allotbspecies_boxplot.pdf", width = 9, height = 5)

#** phylo-signal ----
# longevity
phylosigdata <- log(pgls_allotb$spmean_long_days)
names(phylosigdata) <- rownames(pgls_allotb)

# Pagel's lambda
phytools::phylosig(allotb, x = phylosigdata, method = "lambda", test = TRUE)
# Phylogenetic signal lambda : 0.806212 
# logL(lambda) : -1709.47 
# LR(lambda=0) : 829.31 
# P-value (based on LR test) : 2.28926e-182 

# lambda > 0.8, intermediate to strong phylogenetic signal

# Blomberg et al's K
phytools::phylosig(allotb, x = phylosigdata, method = "K", test = TRUE)
# Phylogenetic signal K : 0.125343 
# P-value (based on 1000 randomizations) : 0.001

# p-value suggests some phylogenetic signal in floral longevity, though <1 so
# less phylogenetic signal than expected under Brownian motion evolution

# now symmetry, binary so use Fritz and Purvis' D
phylosigdata <- as.data.frame(pgls_allotb$sym_species)
names(phylosigdata) <- "symmetry"
phylosigdata$treenames <- rownames(pgls_allotb)

# first remove node labels which confuse caper
allotb2 <- allotb
allotb2$node.label <- NULL

caper::phylo.d(data = phylosigdata, 
               phy = allotb2, 
               names.col = treenames, 
               binvar = symmetry)
# Calculation of D statistic for the phylogenetic structure of a binary variable

# Data :  data
# Binary variable :  symmetry
# Counts of states:  0 = 976
#                    1 = 457
# Phylogeny :  phy
# Number of permutations :  1000
# 
# Estimated D :  -0.304861
# Probability of E(D) resulting from no (random) phylogenetic structure :  0
# Probability of E(D) resulting from Brownian phylogenetic structure    :  1
rm(phylosigdata, allotb2)

#** ttest & PGLS ----
# t-test of longevity by symmetry
ttests <- list()
ttests$allotbspecies <- t.test(pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "1"], 
                               pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "0"])
ttests$allotbspecies
# Welch Two Sample t-test
# 
# data:  pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "1"] and pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "0"]
# t = 4.3817, df = 720.42, p-value = 1.352e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.6162555 1.6167801
# sample estimates:
#   mean of x mean of y 
# 4.679864  3.563347 

# exact results vary with random subsampling

pgls_models <- list()
pgls_models$allotbspecies <- nlme::gls(log(spmean_long_days) ~ sym_species,
                                       correlation = ape::corBrownian(phy = allotb, 
                                                                      form = ~spp),
                                       data = pgls_allotb, method = "ML")
summary(pgls_models$allotbspecies)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ sym_species 
# Data: pgls_allotb 
# AIC      BIC    logLik
# 4595.938 4611.741 -2294.969
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  0.8689184 0.7578723 1.146523  0.2518
# sym_species1 0.3083904 0.0894011 3.449514  0.0006
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.021
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.08674314 -0.28062157 -0.01482579  0.22145390  0.80406738 
# 
# Residual standard error: 3.096406 
# Degrees of freedom: 1433 total; 1431 residual

rm(spp, pgls_allotb)

#* GBOTB PGLS ----

#** prep data ----
# remove gbotb duplicates from data, randomly choose best matches
pgls_gbotb <- spmean_long
# build column of match ranking for gbotb matches
pgls_gbotb$gbotb_matchrank <- ifelse(pgls_gbotb$match_level_gbotb %in% c("direct_accepted", 
                                                               "direct_accepted_nosubsp",
                                                               "manual_misspelling"), 
                                     "1", 
                                     ifelse(pgls_gbotb$match_level_gbotb %in% c("direct_original",
                                                                      "manual_synonym"),
                                            "2",
                                            ifelse(pgls_gbotb$match_level_gbotb == "direct_genus_accepted", 
                                              "3",
                                              ifelse(pgls_gbotb$match_level_gbotb %in% c("direct_genus_original",
                                                                                    "closest_genus",
                                                                                    "genus_synonym"),
                                                     "4", "5"))))
pgls_gbotb <- pgls_gbotb %>%
  dplyr::filter(!is.na(gbotb)) %>%
  dplyr::group_by(gbotb) %>% # group by gbotb species
  dplyr::slice_min(order_by = gbotb_matchrank, n = 1) %>% # choose best possible match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::ungroup()
# 1187 taxa remaining

# reorder data so species order matches order of tips in tree
pgls_gbotb <- as.data.frame(gbotb$tip.label) %>%
  dplyr::left_join(pgls_gbotb, by = c("gbotb$tip.label" = "gbotb"))

# taxon_name to row names
rownames(pgls_gbotb) <- pgls_gbotb[,1]
pgls_gbotb[,1] <- NULL
# taxon names for evolutionary correlation calculation
spp <- rownames(pgls_gbotb)

hist(pgls_gbotb$spmean_long_days)
hist(log(pgls_gbotb$spmean_long_days))
# will log transform longevity to meet assumptions of Brownian motion

#** figure ----
# boxplot of longevity by symmetry, allotb species
ggplot(data = pgls_gbotb, aes(x = spmean_long_days, y = sym_species, fill = sym_species)) +
  geom_violin() +
  geom_boxplot(width = 0.15, colour = "grey", alpha = 0.3) +
  scale_fill_discrete(type = setNames(my_colours$symmetry, c(1, 0))) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("") +
  scale_y_discrete(labels = c("0" = "actinomorphic", "1" = "zygomorphic"))
ggsave("figures/symmetry_longevity_gbotbspecies_boxplot.pdf", width = 9, height = 5)

#** phylo-signal ----
# first longevity
phylosigdata <- log(pgls_gbotb$spmean_long_days)
names(phylosigdata) <- rownames(pgls_gbotb)

# Pagel's lambda
phytools::phylosig(gbotb, x = phylosigdata, method = "lambda", test = TRUE)
# Phylogenetic signal lambda : 0.79123 
# logL(lambda) : -1445.82 
# LR(lambda=0) : 623.6 
# P-value (based on LR test) : 1.23246e-137 

# lambda around 0.8, intermediate phylogenetic signal

# Blomberg et al's K
phytools::phylosig(gbotb, x = phylosigdata, method = "K", test = TRUE)
# Phylogenetic signal K : 0.100588 
# P-value (based on 1000 randomizations) : 0.001 

# p-value suggests some phylogenetic signal in floral longevity, though <1 so
# less phylogenetic signal than expected under Brownian motion evolution

# now symmetry, binary so use Fritz and Purvis' D
phylosigdata <- as.data.frame(pgls_gbotb$sym_species)
names(phylosigdata) <- "symmetry"
phylosigdata$treenames <- rownames(pgls_gbotb)

# first remove node labels which confuse caper
gbotb2 <- gbotb
gbotb2$node.label <- NULL

caper::phylo.d(data = phylosigdata, 
               phy = gbotb2, 
               names.col = treenames, 
               binvar = symmetry)
# Calculation of D statistic for the phylogenetic structure of a binary variable
# 
# Data :  data
# Binary variable :  symmetry
# Counts of states:  0 = 822
# 1 = 365
# Phylogeny :  phy
# Number of permutations :  1000
# 
# Estimated D :  -0.2934379
# Probability of E(D) resulting from no (random) phylogenetic structure :  0
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.997
rm(phylosigdata, gbotb2)

#** ttest & PGLS ----
# t-test of longevity by symmetry
ttests$gbotbspecies <- t.test(pgls_gbotb$spmean_long_days[pgls_gbotb$sym_species == "1"], 
                              pgls_gbotb$spmean_long_days[pgls_gbotb$sym_species == "0"])
ttests$gbotbspecies
# Welch Two Sample t-test
# 
# data:  pgls_gbotb$spmean_long_days[pgls_gbotb$sym_species == "1"] and pgls_gbotb$spmean_long_days[pgls_gbotb$sym_species == "0"]
# t = 4.4206, df = 544.08, p-value = 1.189e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.7318375 1.9023602
# sample estimates:
#   mean of x mean of y 
# 4.942259  3.625160 

# PGLS of longevity by symmetry with phylogeny considered
pgls_models$gbotbspecies <- nlme::gls(log(spmean_long_days) ~ sym_species,
                                      correlation = ape::corBrownian(phy = gbotb,
                                                                     form = ~spp),
                                      data = pgls_gbotb, method = "ML")
summary(pgls_models$gbotbspecies)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ sym_species 
# Data: pgls_gbotb 
# AIC      BIC    logLik
# 4056.306 4071.543 -2025.153
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  0.8261584 0.8517499 0.969954  0.3323
# sym_species1 0.4121765 0.1072434 3.843373  0.0001
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.022
# 
# Standardized residuals:
#   Min           Q1          Med           Q3          Max 
# -0.964298653 -0.239797191  0.001959438  0.205377769  0.674681576 
# 
# Residual standard error: 3.445238 
# Degrees of freedom: 1187 total; 1185 residual

rm(spp, pgls_gbotb)

#* subsampling PGLS ----

# loop to subsample one species per genus for ALLOTB analyses
# loop through, subsample and run ttest and PGLS for each
for (n in 1:50){
  
  # first build data frame, randomly sampling one taxon per genus
  pgls_onepergenus <- spmean_long %>%
    dplyr::group_by(allotb_genus) %>%
    dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # first choose taxa with best taxonomic matches in tree
    dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
    dplyr::ungroup()
  # this selects 804 taxa, one in each genus
  
  # first run straight t-test to get means and sample sizes in each group
  ttests[[paste0("subsample", n)]] <- t.test(pgls_onepergenus$spmean_long_days[pgls_onepergenus$sym_species == "1"],
                                             pgls_onepergenus$spmean_long_days[pgls_onepergenus$sym_species == "0"])
  
  # prune tree to 804 selected taxa
  tree_nomissing <- ape::drop.tip(allotb, allotb$tip.label[-match(pgls_onepergenus$allotb, allotb$tip.label)])
  length(tree_nomissing$tip.label)
  
  # reorder data so species order matches order of tips in tree
  pgls_onepergenus <- as.data.frame(tree_nomissing$tip.label) %>%
    dplyr::left_join(pgls_onepergenus, by = c("tree_nomissing$tip.label" = "allotb"))
  
  # taxon_name to row names
  rownames(pgls_onepergenus) <- pgls_onepergenus[,1]
  pgls_onepergenus[,1] <- NULL
  # taxon names for evolutionary correlation calculation
  spp <- rownames(pgls_onepergenus)
  # run model!
  pgls_models[[paste0("subsample", n)]] <- nlme::gls(log(spmean_long_days) ~ 
                                                       sym_species, 
                                                     correlation = ape::corBrownian(phy = tree_nomissing,
                                                                                    form = ~spp),
                                                     data = pgls_onepergenus, method = "ML")
  
}

rm(n, spp, tree_nomissing, pgls_onepergenus)

#* export results ----
 
# from all PGLS models in simple table with subsamples and their mean

# create data frame for regression results table
pglsresults <- data.frame()

# loop through all results and build table of relevant values
for (model in names(pgls_models)) {
  t <- as.data.frame(table(pgls_models[[model]]$fitted)) # 2 values of fitted indicate n in each symemtry group
  new_row <- tibble(model = model,
                    total_n = pgls_models[[model]]$dims$N,
                    zyg_n = min(t$Freq), # assume zyg will always be smaller n
                    act_n = max(t$Freq), # assume actin always higher n
                    zygmean = ttests[[model]]$estimate[1],
                    actmean = ttests[[model]]$estimate[2],
                    pgls_logLik = pgls_models[[model]]$logLik,
                    pgls_intercept = pgls_models[[model]]$coefficients[1],
                    pgls_symspecies1 = pgls_models[[model]]$coefficients[2],
                    pgls_pvalue = summary(pgls_models[[model]])$tTable[2,4]
  )
  pglsresults <- rbind(pglsresults, new_row)
}

# amazing! 
rm(new_row, t, model)

# export results to csv
readr::write_csv(pglsresults, "results/PGLS_models_key_results.csv")

#* check assumptions ----

# export qqplots for each model
pdf("figures/qqplots_PGLS.pdf", width = 10, height = 10)
par(mfrow = c(8, 7), mar = c(1, 1, 1, 1))
for (model in names(pgls_models)) {
  
  print(qqnorm(resid(pgls_models[[model]]), main = "",
               xlab = "", ylab = "",))
  print(qqline(resid(pgls_models[[model]])))
  
}
dev.off()

# export residual vs fitted plots for each model
pdf("figures/residual_vs_fitted_PGLS.pdf", width = 10, height = 10)
par(mfrow = c(8, 7), mar = c(1, 1, 1, 1))
for (model in names(pgls_models)) {
  
  print(plot(resid(pgls_models[[model]], type = "p") ~ fitted(pgls_models[[model]]),
             col = "mediumblue", xlab = "", ylab = ""))
  abline(0, 0)
  grid()
  
}
dev.off()

# mostly look okay, 1-2 outlying low residuals in actinomorphic category but
# they don't look so outlying as to be hugely influential

rm(pgls_models, model, ttests, pglsresults)

#### plot phylogeny ####

#* prep data ----

# random subsample of one species per allotb match
spmean_long_sub <- spmean_long %>%
  dplyr::select(allotb, family, sym_species, spmean_long_days, allotb_matchrank) %>%
  dplyr::group_by(allotb) %>% # group by allotb species
  dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # choose best possible taxonomic match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::ungroup()

# add column with species' tip position in tree for plotting
treetips <- data.frame(allotb = allotb$tip.label, position = c(1:1433))
spmean_long_sub <- spmean_long_sub %>%
  dplyr::left_join(treetips, by = "allotb")
rm(treetips)

# rename some paraphyletic family labels so they display properly
spmean_long_sub$family[spmean_long_sub$position %in% 589:592] <- "Lamiaceae (a)"
spmean_long_sub$family[spmean_long_sub$position %in% 426:469] <- "Lamiaceae (b)"
spmean_long_sub$family[spmean_long_sub$position == 205] <- "Fabaceae (a)"
spmean_long_sub$family[spmean_long_sub$position %in% 95:201] <- "Fabaceae (b)"
spmean_long_sub$family[spmean_long_sub$position == 886] <- "Primulaceae (a)"
spmean_long_sub$family[spmean_long_sub$position %in% 864:876] <- "Primulaceae (b)"

# change longevity data into a named vector for phytools
spmean_long_subV <- log(spmean_long_sub$spmean_long_days)
names(spmean_long_subV) <- spmean_long_sub$allotb
# and prep symmetry data also
symV <- as.factor(spmean_long_sub$sym_species)
names(symV) <- spmean_long_sub$allotb
# make sure discrete character is in the order of tree
symV <- symV[allotb$tip.label]
# set factor colours
cols <- setNames(c(my_colours$symmetry[2], my_colours$symmetry[1]), c("0", "1"))

# simulate longevity evolution across phylogeny to visualise
contmap <- phytools::contMap(allotb, spmean_long_subV, plot = FALSE, res = 500)
# re-colour contmap with custom scale
contmap <- phytools::setMap(contmap, my_colours$longevity)

#* tall plot ----

# dummy plot to get locations to draw tip labels coloured by symmetry
plot(contmap, fsize = c(0.5, 0.7))
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# export figure
pdf(file = "figures/tall_contmap_spmeanlongevity.pdf", width = 20, height = 120)
plot(contmap, legend = 0.7*max(nodeHeights(allotb)), sig = 1, 
     lwd = 6, outline = FALSE, ftype = "off",
     xlim = c(lastPP$x.lim[1], lastPP$x.lim[2] + 20),
     leg.txt = "Floral longevity (log mean # days)")
for(i in 1:length(symV)) {
  text(lastPP$xx[i], lastPP$yy[i], 
       gsub("_", " ", contmap$tree$tip.label[i]),
       pos = 4, cex = 0.6, col = cols[symV[i]], font = 3)
}

# function to work out highest and lowest tip numbers for each family,
# then loop through these to draw segment and text labels for each family
family_labels <- function(xpos = 200){
  tip_pos <- spmean_long_sub %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(y0 = min(position), y1 = max(position))
  for (i in 1:nrow(tip_pos)){
    segments(x0 = xpos, y0 = tip_pos$y0[[i]], x1 = xpos, y1 = tip_pos$y1[[i]], lwd = 2)
    text(x = xpos+1, y = tip_pos$y0[[i]] + ((tip_pos$y1[[i]] - tip_pos$y0[[i]])/2), 
         paste(tip_pos$family[[i]], sep = ""), adj = 0, srt = 0, cex = 0.5)
  }
  rm(tip_pos)
}

family_labels(xpos = 163)

# function to work out highest and lowest tip numbers for each order,
# then loop through these to draw segment and text labels for each order
order_labels <- function(xpos = 200){
  tip_pos <- spmean_long_sub %>%
    dplyr::group_by(order) %>%
    dplyr::summarise(y0 = min(position), y1 = max(position))
  for (i in 1:nrow(tip_pos)){
    segments(x0 = xpos, y0 = tip_pos$y0[[i]], x1 = xpos, y1 = tip_pos$y1[[i]], lwd = 2)
    text(x = xpos + 1, y = tip_pos$y0[[i]] + ((tip_pos$y1[[i]] - tip_pos$y0[[i]])/2), 
         paste(tip_pos$order[[i]], sep = ""), adj = 0, srt = 0, cex = 0.5)
  }
  rm(tip_pos)
}

order_labels(xpos = 171.5)

# label clades
tall_cladelabels <- function(xpos = 185){
  # clade labelling as per Ramirez-Barahona et al. (2020)
  segments(x0 = xpos, y0 = 1432, x1 = xpos, y1 = 1433, lwd = 3, col = "#BBCDE9")
  text(x = xpos+1.5, y = 1432.5, "ANA", srt = 0, adj = 0, cex = 1, col = "#BBCDE9")
  segments(x0 = xpos, y0 = 1388, x1 = xpos, y1 = 1431, lwd = 3, col = "#47A1D1")
  text(x = xpos+2, y = 1409.5, "Magnoliids", srt = 270, cex = 1, col = "#47A1D1")
  segments(x0 = xpos, y0 = 1135, x1 = xpos, y1 = 1387, lwd = 3, col = "#59BE1C")
  text(x = xpos+2, y = 1261, "Monocots", srt = 270, cex = 1, col = "#59BE1C")
  segments(x0 = xpos-5.5, y0 = 1135, x1 = xpos-5.5, y1 = 1229, lwd = 3, col = "#0C9934")
  text(x = xpos-3.5, y = 1182, "Commelinids", srt = 270, cex = 1, col = "#0C9934")
  segments(x0 = xpos, y0 = 1, x1 = xpos, y1 = 1134, lwd = 3, col = "#F0D01B")
  text(x = xpos+2, y = 567.5, "Eudicots", srt = 270, cex = 1, col = "#F0D01B")
  segments(x0 = xpos-5.5, y0 = 1, x1 = xpos-5.5, y1 = 376, lwd = 3, col = "#F2B211")
  text(x = xpos-3.5, y = 188, "Rosids", srt = 270, cex = 1, col = "#F2B211")
  segments(x0 = xpos-5.5, y0 = 399, x1 = xpos-5.5, y1 = 956, lwd = 3, col = "#FCB98E")
  text(x = xpos-3.5, y = 677.5, "Asterids", srt = 270, cex = 1, col = "#FCB98E")
}

tall_cladelabels()

# insert legend
legend(x = "bottomright", legend = c("actinomorphic", "zygomorphic"), bg = "white",
       fill = cols, cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Floral symmetry")
dev.off()

rm(lastPP, family_labels, order_labels, tall_cladelabels, i)

#* circular plot ----

# build fan style plot for main figure
pdf(file = "figures/allotb_longevity_contMap_fan.pdf", width = 15, height = 15)

plot(contmap, type = "fan", legend = FALSE, lwd = 3, outline = FALSE, 
     ftype = "off", xlim = c(-185, 170))

# fill background in so pale colours in contMap stand out more clearly
plotrix::draw.circle(0, 0, radius = max(nodeHeights(allotb)), 
                     col = "#dadada", lty = 0)
# label end Cretaceous period at 66 mya
plotrix::draw.circle(0, 0, radius = max(nodeHeights(allotb)) - 66, 
                     border = "black", lty = 3)

# plot contMap again
par(new = TRUE) # hack to force below to plot on top of above 
plot(contmap, type = "fan", legend = FALSE, lwd = 3, outline = FALSE, 
     ftype = "off", xlim = c(-185, 170), add = TRUE)

# below adapted from http://blog.phytools.org/2016/08/vertical-legend-in-contmap-style-plots.html
# add bud size legend using phytools function
phytools::add.color.bar(leg = 100,
                        cols = contmap$cols,
                        title = "Floral longevity (log days)",
                        lims = NULL,
                        digits = 2,
                        prompt = FALSE,
                        lwd = 12,
                        direction = "upwards",
                        subtitle = "",
                        x = -180,
                        y = 45)

# then add custom tick marks
lines(x = rep(-176.525, 2), y = c(45, 147)) # draw vertical line
Y <- cbind(seq(45, 147, length.out = 5), # define x pos for ticks
           seq(45, 147, length.out = 5))
X <- cbind(rep(-176.525, 5), # define y pos for ticks
           rep(-173.925, 5))
for(i in 1:nrow(Y)) lines(X[i,], Y[i,]) # draw ticks
ticks <- seq(contmap$lims[1], contmap$lims[2], length.out = 5) # get tick values
text(x = X[,2], y = Y[,2], round(ticks, 1), pos = 4, cex = 0.8) # draw tick values
rm(X, Y, i, ticks)

# now to plot tip points coloured by flower symmetry

# assign plot to object
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
# from pp object will pull out x and y coordinates to plot points

#function below by GM to offset x and y points, huzzah
offset_xx_yy <- function(xx, yy, offset) {
  angle <- atan2(yy, xx)
  data.frame(
    xx = xx + offset * cos(angle),
    yy = yy + offset * sin(angle)
  )
}

xx_yy <- offset_xx_yy(
  xx = pp$xx[1:ape::Ntip(allotb)],
  yy = pp$yy[1:ape::Ntip(allotb)],
  offset = 4.5
)

# add flower symmetry points
points(xx_yy$xx,
       xx_yy$yy,
       pch = 16, cex = 2,
       col = cols[symV[allotb$tip.label]])

legend(x = 130, y = 150, legend = c("actinomorphic", "zygomorphic"), col = cols, 
       bty = "n", cex = 0.8, title = "Flower symmetry", pch = 19)

#** clade labelling ----

# label larger orders
# first calculate order tip positions
orders_to_label <- spmean_long_sub %>%
  dplyr::select(position, order) %>%
  dplyr::group_by(order) %>%
  tidyr::nest() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mintip = purrr::map(data, ~{min(.x$position)})) %>%
  dplyr::mutate(maxtip = purrr::map(data, ~{max(.x$position)})) %>%
  tidyr::unnest(cols = c(data, mintip, maxtip)) %>%
  dplyr::select(order, mintip, maxtip) %>%
  dplyr::distinct() %>%
  dplyr::mutate(tiprange = maxtip - mintip)

# then filter out smaller orders AND a few that overlap for labelling, 
# leaves 21 of 43 orders 
orders_to_label <- orders_to_label %>%
  dplyr::filter(tiprange > 10 & !(order %in% c("Dipsacales", "Saxifragales", "Sapindales")))
# and replace overlapping orders with abbreviations
orders_to_label$order <- gsub("Santalales", "Santal.", orders_to_label$order)
orders_to_label$order <- gsub("Apiales", "Api.", orders_to_label$order)
orders_to_label$order <- gsub("Rosales", "Ros.", orders_to_label$order)
orders_to_label$order <- gsub("Magnoliales", "Magnoli.", orders_to_label$order)
orders_to_label$order <- gsub("Brassicales", "Brassic.", orders_to_label$order)
orders_to_label$order <- gsub("Boraginales", "Boragin.", orders_to_label$order)

# source custom arc labelling function
source("scripts/functions/arclabel.R")

# loop through and draw labels on phylogeny for larger orders
for(i in 1:length(orders_to_label$order)) {
  arclabel(text = orders_to_label$order[i],
           tips = c(orders_to_label$mintip[i], orders_to_label$maxtip[i]),
           cex = 1,
           ln.offset = 1.05,
           lab.offset = 1.069)
}
rm(i, orders_to_label)

# clade labelling using custom arclabel function for circular (fan) tree

fan_cladelabels <- function(offset = 1, cex = 1.5){ # use offset argument to move labels closer (<1) or further away (>1) from tree
  source("scripts/functions/arclabel.R") # get arclabel function
  arclabel(text = "ANA", tips = c(1432, 1433),
           lwd = 15, cex = cex, col = "#636363",
           ln.offset = offset + .065, lab.offset = offset + .13,
           orientation = "perpendicular")
  arclabel(text = "Magnoliids", tips = c(1388, 1431), 
           lwd = 10, cex = cex-0.3, col = "#dadada",
           ln.offset = offset + .065, lab.offset = offset + .10)
  arclabel(text = "Monocots", tips = c(1135, 1387), 
           lwd = 10, cex = cex, col = "#636363",
           ln.offset = offset + .065, lab.offset = offset + .10)
  arclabel(text = "Commelinids", tips = c(1135, 1229), 
           lwd = 4, cex = cex, col = "#dadada",
           ln.offset = offset + .055, lab.offset = offset + .10)
  arclabel(text = "Eudicots", tips = c(1, 1134),
           lwd = 10, cex = cex, col = "#636363",
           ln.offset = offset + .065, lab.offset = offset + .10)
  arclabel(text = "Rosids", tips = c(1, 376),
           lwd = 4, cex = cex, col = "#dadada",
           ln.offset = offset + .055, lab.offset = offset + .10)
  arclabel(text = "Asterids", tips = c(399, 956), 
           lwd = 4, cex = cex, col = "#dadada",
           ln.offset = offset + .055, lab.offset = offset + .10)
}

fan_cladelabels(offset = 1.05)

dev.off()

rm(spmean_long_sub, pp, contmap, cols, fan_cladelabels, offset_xx_yy, xx_yy,
   arclabel, spmean_long, spmean_long_subV, symV, allotb, gbotb)
