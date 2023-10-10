# this script adapts Smith and Brown (2018) evolutionary trees to sample taxa,
# runs phylogenetic generalised least squares regressions (PGLS) to see if 
# flowers have evolved to last longer or shorter in taxa with zygomorphic 
# flowers, and test the phylogenetic signal and evolutionary rate of floral 
# symmetry and floral longevity evolution

#### get trees and match taxonomy ####

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
pgls_data <- sym_long %>%
  dplyr::group_by(species = Accepted_name, sym_species) %>%
  dplyr::summarise(spmean_long_days = mean(mean_long_days)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(species)) %>%
  as.data.frame()
# add underscore to match tip labels
pgls_data$species <- gsub(" ", "_", pgls_data$species)

# match allotb and gbotb names to this data
phylo_names_match <- phylo_names_match %>%
  dplyr::select(species, genus:match_level_gbotb) %>%
  dplyr::distinct()
# above df has 3 rows more than pgls_data because of different og name matches giving
# different match levels, will leave these duplicates for now and should be
# fixed when I choose one species per genus
pgls_data <- pgls_data %>%
  dplyr::left_join(phylo_names_match, by = "species")
rm(phylo_names_match)

# redefine symmetry as 0 and 1
pgls_data$sym_species <- gsub("zygomorphic", "1", pgls_data$sym_species)
pgls_data$sym_species <- gsub("actinomorphic", "0", pgls_data$sym_species)
table(pgls_data$sym_species)
# 986 actinomorphic to 469 zygomorphic taxa (will change with subsampling)

# define heirarchy of taxonomic matches so can sample best matches first
pgls_data$allotb_matchrank <- ifelse(pgls_data$match_level_allotb %in% c("direct_accepted", 
                                                                         "direct_accepted_nosubsp",
                                                                         "manual_misspelling"), 
                                     "1", 
                                     ifelse(pgls_data$match_level_allotb %in% c("direct_original",
                                                                                "manual_synonym"), 
                                            "2",
                                            ifelse(pgls_data$match_level_allotb == "direct_genus_accepted", 
                                                   "3",
                                                   ifelse(pgls_data$match_level_allotb %in% c("direct_genus_original",
                                                                                              "closest_genus",
                                                                                              "genus_synonym"), 
                                                          "4", "5"))))
table(pgls_data$allotb_matchrank)

#* PGLS without subsampling, ALLOTB ----

# remove allotb duplicate matches from data, randomly choose best matches
pgls_allotb <- pgls_data

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

# boxplot of longevity by symmetry, allotb species
ggplot(data = pgls_allotb, aes(x = log(spmean_long_days), y = sym_species, fill = sym_species)) +
  geom_boxplot() +
  scale_fill_discrete(type = setNames(my_colours$symmetry, c(1, 0))) +
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
# t = 4.4473, df = 864.12, p-value = 9.831e-06
# mean of x   mean of y 
# 1.0397000 0.7696431 

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
# 4594.801 4610.604 -2294.401
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  0.8684088 0.7575712 1.146307  0.2519
# sym_species1 0.3132770 0.0894827 3.500979  0.0005
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.021
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -1.08700984 -0.28056831 -0.01624582  0.22079506  0.80455117 
# 
# Residual standard error: 3.095178 
# Degrees of freedom: 1433 total; 1431 residual

rm(spp, pgls_allotb)

#* PGLS without subsampling, GBOTB ----

# remove gbotb duplicates from data, randomly choose best matches
pgls_gbotb <- pgls_data
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

# boxplot of longevity by symmetry, allotb species
ggplot(data = pgls_gbotb, aes(x = log(spmean_long_days), y = sym_species, fill = sym_species)) +
  geom_boxplot() +
  scale_fill_discrete(type = setNames(my_colours$symmetry, c(1, 0))) +
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
# t = 4.2181, df = 665.12, p-value = 2.806e-05
# mean of x   mean of y 
# 1.0842894 0.7983383 

# PGLS of longevity by symmetry with phylogeny considered
pgls_models$gbotbspecies <- nlme::gls(log(spmean_long_days) ~ sym_species,
                                      correlation = ape::corBrownian(phy = gbotb,
                                                                     form = ~spp),
                                      data = pgls_gbotb, method = "ML")
summary(pgls_models$gbotbspecies)
# Generalized least squares fit by maximum likelihood
# Model: log(spmean_long_days) ~ sym_species 
# Data: pgls_gbotb 
# AIC     BIC    logLik
# 4060.513 4075.75 -2027.256
# 
# Correlation Structure: corBrownian
# Formula: ~spp 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept)  0.8170607 0.8532520 0.957584  0.3385
# sym_species1 0.3710309 0.1050007 3.533602  0.0004
# 
# Correlation: 
#   (Intr)
# sym_species1 -0.022
# 
# Standardized residuals:
#   Min          Q1         Med          Q3         Max 
# -0.95995534 -0.23673665  0.01692324  0.21957171  0.67612298 
# 
# Residual standard error: 3.451349 
# Degrees of freedom: 1187 total; 1185 residual

rm(spp, pgls_gbotb)

#* PGLS with subsampling ----

# loop to subsample one species per genus for ALLOTB analyses
# first create column of allotb genus
pgls_data$allotb_genus <- gsub("_.*", "", pgls_data$allotb)

# loop through to subsample one per genus and run ttest and PGLS for each
for (n in 1:50){
  
  # first build data frame, randomly sampling one taxon per genus
  pgls_onepergenus <- pgls_data %>%
    dplyr::group_by(allotb_genus) %>%
    dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # first choose taxa with best taxonomic matches in tree
    dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
    dplyr::ungroup()
  # this selects 804 taxa, one in each genus
  
  # first run straight t-test to get means and sample sizes in each group
  ttests[[paste0("subsample", n)]] <- t.test(log(pgls_onepergenus$spmean_long_days[pgls_onepergenus$sym_species == "1"]),
                                             log(pgls_onepergenus$spmean_long_days[pgls_onepergenus$sym_species == "0"]))
  
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

# export results from all PGLS models in simple table with subsamples and their mean

#create data frame for regression results table
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

#** check model assumptions ####

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


#### model longevity evolution ####

# prep data for modelling
longevspmean <- pgls_data %>%
  dplyr::select(allotb, sym_species, spmean_long_days, allotb_matchrank) %>%
  dplyr::group_by(allotb) %>% # group by allotb species
  dplyr::slice_min(order_by = allotb_matchrank, n = 1) %>% # choose best possible taxonomic match
  dplyr::slice_sample(n = 1) %>% # then randomly choose one of these
  dplyr::ungroup()

# change longevity data into a named vector for phytools
longevspmeanV <- log(longevspmean$spmean_long_days)
names(longevspmeanV) <- longevspmean$allotb
# and prep symmetry data also
symV <- as.factor(longevspmean$sym_species)
names(symV) <- longevspmean$allotb

# visualise longevity evolution across phylogeny
contmap <- phytools::contMap(allotb, longevspmeanV, plot = FALSE)
# re-colour contmap with custom scale
contmap <- phytools::setMap(contmap, my_colours$longevity)

# make sure discrete character is in the order of tree
symV <- symV[contmap$tree$tip.label]
# set factor colours
cols <- setNames(c(my_colours$symmetry[2], my_colours$symmetry[1]), c("0", "1"))

# dummy plot to get locations to draw tip labels coloured by symmetry
plot(contmap, fsize = c(0.5, 0.7))
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# and plot it! v tall plot with tip labels
pdf(file = "figures/contmap_spmeanlongevity.pdf", width = 20, height = 120)
plot(contmap, legend = 0.7*max(nodeHeights(allotb)), sig = 1, 
     lwd = 4, outline = FALSE, ftype = "off", #type = "fan",
     xlim = lastPP$x.lim, #ylim = lastPP$y.lim,
     leg.txt = "Floral longevity (log mean # days)")
for(i in 1:length(symV)) {
  text(lastPP$xx[i], lastPP$yy[i], contmap$tree$tip.label[i],
       pos = 4, cex = 0.6, col = cols[symV[i]], font = 3)
}
# insert legend
legend(x = "bottomright", legend = c("actinomorphic", "zygomorphic"), bg = "white",
       fill = cols, cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Floral symmetry")
dev.off()
rm(lastPP)

#### CIRCULAR PLOT ####

# build fan style plot for main figure
#pdf(file = "figures/allotb_longevity_contMap_fan.pdf", width = 20, height = 20)

plot(contmap, type = "fan", legend = FALSE, lwd = 4, outline = FALSE, 
     ftype = "off", xlim = c(-180, 140))

# label Cretaceous period, 139 to 66 mya shown here
plotrix::draw.circle(0, 0, radius = max(nodeHeights(allotb)) - 66, 
                     col = "#dadada", lty = 0)

# plot contMap again
par(new = TRUE) # hack to force below to plot on top of above 
plot(contmap, type = "fan", legend = FALSE, lwd = 4, outline = FALSE, 
     ftype = "off", xlim = c(-180, 140), add = TRUE)

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
                        y = -50)

# then add custom tick marks
lines(x = rep(-175.525, 2), y = c(-50, 52)) # draw vertical line
Y <- cbind(seq(-50, 52, length.out = 5), # define x pos for ticks
           seq(-50, 52, length.out = 5))
X <- cbind(rep(-175.525, 5), # define y pos for ticks
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
  offset = 5
)

# add flower symmetry points
points(xx_yy$xx,
       xx_yy$yy,
       pch = 15, cex = 0.5,
       col = cols[symV[allotb$tip.label]])
# will need to use trigonometry if want to offset these points at all

legend(x = "topright", legend = c("actinomorphic", "zygomorphic"), col = cols, 
       bty = "n", cex = 0.8, title = "Flower symmetry", pch = 15)

#dev.off()

# TO DO
# - label clades
# - offset or enlarge tip points?


rm(symV, longevspmean, lastPP, contmap, longevspmeanV, i, cols)


