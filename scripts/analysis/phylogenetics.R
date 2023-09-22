# script to generate phylogeny for taxa from Smith and Brown (2018) tree
# and run some basic phylogenetic analyses to see if floral longevity and
# floral symmetry are co-evolving traits

# TO INVESTIGATE
# - randomly sample one species per genus. Same results?
# - use shorter Smith and Brown tree (GBOTB.tre with 79,881 tips)
# - use Smith and Brown without Qian and Jin?? (ALLOTB.tre with 353,185 tips)

#### get trees and taxonomic match data ####

# read in short Smith and Brown tree (GBOTB.tre with 79,881 tips)
gbotb <- ape::read.tree("data_input/GBOTB.tre")
# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# read in taxonomic name matching key
source("scripts/prepdata/phylo_names_match.R")

#### prune trees to matched names ####

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

#### run analyses ####

# summarise longevity per ACCEPTED species
# take mean per species for longevity
pgls <- sym_long %>%
  dplyr::group_by(species = Accepted_name, sym_species) %>%
  dplyr::summarise(spmean_long_days = mean(mean_long_days)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(species)) %>%
  as.data.frame()
# add underscore to match tip labels
pgls$species <- gsub(" ", "_", pgls$species)

# match allotb and gbotb names to this data
phylo_names_match <- phylo_names_match %>%
  dplyr::select(species, genus:match_level_gbotb) %>%
  dplyr::distinct()
# above df has 3 rows more than pgls because of different og name matches giving
# different match levels, will leave these duplicates for now and should be
# fixed when I choose one species per genus
pgls <- pgls %>%
  dplyr::left_join(phylo_names_match, by = "species")
rm(phylo_names_match)

# redefine symmetry as 0 and 1
pgls$sym_species <- gsub("zygomorphic", "1", pgls$sym_species)
pgls$sym_species <- gsub("actinomorphic", "0", pgls$sym_species)
table(pgls$sym_species)
# 986 actinomorphic to 469 zygomorphic taxa (will change with subsampling, record in results)

#* PGLS without subsampling, ALLOTB ----

# remove allotb duplicate matches from data, randomly choose best matches
pgls_allotb <- pgls
pgls_allotb$allotbgenusmatch <- ifelse(pgls_allotb$match_level_allotb %in% c("direct_genus_accepted",
                                                                             "direct_genus_original",
                                                                             "closest_genus",
                                                                             "genus_synonym"),
                                       "1", "0")
pgls_allotb <- pgls_allotb %>%
  dplyr::group_by(allotb) %>% # group by allotb species
  dplyr::slice_min(order_by = allotbgenusmatch, n = 1) %>% # if possible choose taxa without genus level match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::ungroup()

# reorder data so species order matches order of tips in tree
pgls_allotb <- as.data.frame(allotb$tip.label) %>%
  dplyr::left_join(pgls_allotb, by = c("allotb$tip.label" = "allotb"))

# taxon_name to row names
rownames(pgls_allotb) <- pgls_allotb[,1]
pgls_allotb[,1] <- NULL

hist(pgls_allotb$spmean_long_days)
hist(log(pgls_allotb$spmean_long_days))
# H says should log transform because of assumptions of Brownian motion

# boxplot of longevity by symmetry, allotb species
ggplot(data = pgls_allotb, aes(x = log(spmean_long_days), y = sym_species, fill = sym_species)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (log days)") +
  ylab("") +
  scale_y_discrete(labels = c("0" = "actinomorphic", "1" = "zygomorphic"))
ggsave("figures/symmetry_longevity_allotbspecies_boxplot.pdf", width = 9, height = 5)

# t-test of longevity by symmetry
ttests <- list()
ttests$allotbspecies <- t.test(log(pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "1"]), 
                               log(pgls_allotb$spmean_long_days[pgls_allotb$sym_species == "0"]))
ttests$allotbspecies
# t = 4.3965, df = 865.78, p-value = 1.237e-05
# mean of x   mean of y 
# 1.0370836 0.7698605 

pglsresults <- list()
pglsresults$allotbspecies <- nlme::gls(log(spmean_long_days) ~ sym_species,
                                       correlation = ape::corBrownian(phy = allotb),
                                       data = pgls_allotb, method = "ML")
summary(pglsresults$allotbspecies)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ sym_species 
# Data: pgls_allotb 
# AIC      BIC    logLik
# 4598.858 4614.661 -2296.429
# 
# Correlation Structure: corBrownian
# Formula: ~1 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  0.8678011 0.7586449 1.143883  0.2529
# sym_species1 0.3165041 0.0894922 3.536666  0.0004
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.021
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.08527600 -0.27997534 -0.01706794  0.22158885  0.80360904 
# 
# Residual standard error: 3.099563 
# Degrees of freedom: 1433 total; 1431 residual
plot(log(spmean_long_days) ~ sym_species, data = pgls_allotb)

#* PGLS without subsampling, GBOTB ----

# remove gbotb duplicates from data, randomly choose best matches
pgls_gbotb <- pgls
# build column of match ranking for gbotb matches
pgls_gbotb$gbotb_matchrank <- ifelse(pgls$match_level_gbotb %in% c("direct_accepted", 
                                                               "direct_accepted_nosubsp",
                                                               "manual_misspelling"), 
                                     "1", 
                                     ifelse(pgls$match_level_gbotb %in% c("direct_original",
                                                                      "manual_synonym"),
                                            "2",
                                            ifelse(pgls$match_level_gbotb == "direct_genus_accepted", 
                                              "3",
                                              ifelse(pgls$match_level_gbotb %in% c("direct_genus_original",
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

hist(pgls_gbotb$spmean_long_days)
hist(log(pgls_gbotb$spmean_long_days))
# H says should log transform because of assumptions of Brownian motion

# boxplot of longevity by symmetry, allotb species
ggplot(data = pgls_gbotb, aes(x = log(spmean_long_days), y = sym_species, fill = sym_species)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (log days)") +
  ylab("") +
  scale_y_discrete(labels = c("0" = "actinomorphic", "1" = "zygomorphic"))
ggsave("figures/symmetry_longevity_gbotbspecies_boxplot.pdf", width = 9, height = 5)

# t-test of longevity by symmetry
ttests$gbotbspecies <- t.test(log(pgls_gbotb$spmean_long_days[pgls_gbotb$sym_species == "1"]), 
                              log(pgls_gbotb$spmean_long_days[pgls_gbotb$sym_species == "0"]))
ttests$gbotbspecies
# t = 4.1448, df = 665.13, p-value = 3.841e-05
# mean of x   mean of y 
# 1.0814189   0.7990447 

# PGLS of longevity by symmetry with phylogeny considered
pglsresults$gbotbspecies <- nlme::gls(log(spmean_long_days) ~ sym_species,
                                      correlation = ape::corBrownian(phy = gbotb),
                                      data = pgls_gbotb, method = "ML")
summary(pglsresults$gbotbspecies)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ sym_species 
# Data: pgls_gbotb 
# AIC      BIC    logLik
# 4037.224 4052.462 -2015.612
# 
# Correlation Structure: corBrownian
# Formula: ~1 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  0.8420417 0.8449287 0.996583  0.3192
# sym_species1 0.3897138 0.1058828 3.680614  0.0002
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.022
# 
# Standardized residuals:
#   Min           Q1          Med           Q3          Max 
# -0.976728059 -0.246379808  0.006146746  0.208960294  0.736351219 
# 
# Residual standard error: 3.417657 
# Degrees of freedom: 1187 total; 1185 residual
plot(log(spmean_long_days) ~ sym_species, data = pgls_gbotb)

#* PGLS with subsampling ----

# then loop to subsample one species per genus for ALLOTB analyses
# first create column of allotb genus
pgls$allotb_genus <- gsub("_.*", "", pgls$allotb)
# and build column of match ranking for allotb matches
pgls$allotb_matchrank <- ifelse(pgls$match_level_allotb %in% c("direct_accepted", 
                                                               "direct_accepted_nosubsp",
                                                               "manual_misspelling"), 
                                "1", 
                                ifelse(pgls$match_level_allotb %in% c("direct_original",
                                                                      "manual_synonym"), 
                                       "2",
                                       ifelse(pgls$match_level_allotb == "direct_genus_accepted", 
                                              "3",
                                              ifelse(pgls$match_level_allotb %in% c("direct_genus_original",
                                                                                    "closest_genus",
                                                                                    "genus_synonym"), 
                                                     "4", "5"))))
table(pgls$allotb_matchrank)

# now loop!
for (n in 1:10){
  
  # first build data frame, randomly sampling one taxon per genus
  pgls_onepergenus <- pgls %>%
    dplyr::group_by(allotb_genus) %>%
    dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # first choose taxa with best taxonomic matches in tree (THOUGH THIS WILL MEAN SOME DATA POINTS NEVER USED EVEN WITH RANDOMISATION)
    dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
    dplyr::ungroup()
  # this selects 804 taxa, one in each genus
  
  # prune tree to these 804
  tree_nomissing <- ape::drop.tip(allotb, allotb$tip.label[-match(pgls_onepergenus$allotb, allotb$tip.label)])
  length(tree_nomissing$tip.label)
  
  # reorder data so species order matches order of tips in tree
  pgls_onepergenus <- as.data.frame(tree_nomissing$tip.label) %>%
    dplyr::left_join(pgls_onepergenus, by = c("tree_nomissing$tip.label" = "allotb"))
  #pgls_onepergenus$sym_species <- forcats::as_factor(pgls_onepergenus$sym_species)
  
  # taxon_name to row names
  rownames(pgls_onepergenus) <- pgls_onepergenus[,1]
  pgls_onepergenus[,1] <- NULL
  
  table(pgls_onepergenus$sym_species)
  # ~566 actinomorphic to ~238 zygomorphic taxa (will change with subsampling, record in results)
  
  # might need to find way to automatically check and report on variable
  # distribution in different subsamples, can't eyeball easily
  
  pglsresults[[paste0("subsample", n)]] <- nlme::gls(log(spmean_long_days) ~ sym_species, 
                                correlation = ape::corBrownian(phy = tree_nomissing),
                                data = pgls_onepergenus, method = "ML")
  
}

# next - export results from all subsampling
# find ways to check model assumptions??
summary(pglsresults$subsample1)
plot(pglsresults$subsample1$residuals)
summary(pglsresults$subsample2)
summary(pglsresults$subsample3)
summary(pglsresults$subsample4)
summary(pglsresults$subsample5)
summary(pglsresults$subsample6)
summary(pglsresults$subsample7)
summary(pglsresults$subsample8)
summary(pglsresults$subsample9)
summary(pglsresults$subsample10)

#### analyses with whole data ####


pglsresults$``


#### model longevity evolution ####

# first summarise longevity mean per species
longevspmean <- sym_long %>%
  dplyr::group_by(species = Accepted_name) %>%
  dplyr::summarise(spmean_long_days = mean(mean_long_days)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(species))
# join symmetry data back on
sym <- sym_long %>%
  dplyr::select(species = Accepted_name, sym_species) %>%
  dplyr::filter(!is.na(species)) %>%
  dplyr::distinct()
longevspmean <- longevspmean %>%
  dplyr::left_join(sym, by = "species")
rm(sym)
longevspmean$sym_species[is.na(longevspmean$sym_species)] <- "blank" # add dummy category for NA symmetry
# add underscore to match tip labels
longevspmean$species <- gsub(" ", "_", longevspmean$species)
# change longevity data into a named vector for phytools
longevspmeanV <- log(longevspmean$spmean_long_days)
names(longevspmeanV) <- longevspmean$species
# and prep symmetry data also
symV <- as.factor(longevspmean$sym_species)
names(symV) <- longevspmean$species

# visualise longevity evolution across phylogeny
contmap <- phytools::contMap(tree_TPL$scenario.3, longevspmeanV, plot = FALSE)
# re-colour contmap with viridis scale
contmap <- phytools::setMap(contmap, c("#440154FF", "#46337EFF", "#365C8DFF", "#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF"))

# make sure discrete character is in the order of tree
symV <- symV[contmap$tree$tip.label]
# set factor colours
cols <- setNames(c("#7570B3", "#D95F02"), levels(symV))

# dummy plot to get locations to draw tip labels coloured by symmetry
plot(contmap, fsize = c(0.5, 0.7))
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# have to figure out type = "fan" later, doesn't draw tips in right direction
# but can probs do tip points easily

# and plot it! v tall plot with tip labels
pdf(file = "figures/contmap_spmeanlongevity.pdf", width = 20, height = 80)
plot(contmap, legend = 0.7*max(nodeHeights(tree_TPL$scenario.3)), sig = 1, 
     lwd = 4, outline = FALSE, ftype = "off", #type = "fan",
     xlim = lastPP$x.lim, #ylim = lastPP$y.lim,
     leg.txt = "Floral longevity (log mean # days)")
for(i in 1:length(symV)) {
  text(lastPP$xx[i], lastPP$yy[i], contmap$tree$tip.label[i],
       pos = 4, cex = 0.6, col = cols[symV[i]], font = 3)
}
# insert legend
legend(x = "bottomright", legend = names(cols), bg = "white",
       fill = cols, cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Floral symmetry")
dev.off()

# ultimately think contMap plot would look best in fan style, with points at tips 
# to indicate if taxon is actinomorphic or zygomorphic and round clade labels
# will build this LATER if decide it's worth it
rm(symV, longevspmean, lastPP, contmap, longevspmeanV, i, cols)

rm(tree_TPL)
