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

#### PGLS analyses ####

#* ALLOTB PGLS ----

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

#* GBOTB PGLS ----

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
# add orders and clades for labelling in phylogeny
familyorderclade <- readr::read_csv("data_input/family_order_clade.csv")
spmean_long_sub <- spmean_long_sub %>%
  dplyr::left_join(familyorderclade, by = "family")
rm(familyorderclade)

# rename some paraphyletic family labels so they display properly
spmean_long_sub$family[spmean_long_sub$position %in% 589:592] <- "Lamiaceae (a)"
spmean_long_sub$family[spmean_long_sub$position %in% 426:469] <- "Lamiaceae (b)"
spmean_long_sub$family[spmean_long_sub$position == 205] <- "Fabaceae (a)"
spmean_long_sub$family[spmean_long_sub$position %in% 95:201] <- "Fabaceae (b)"

# change longevity data into a named vector for phytools
spmean_long_subV <- log(spmean_long_sub$spmean_long_days)
names(spmean_long_subV) <- spmean_long_sub$allotb
# and prep symmetry data also
symV <- as.factor(spmean_long_sub$sym_species)
names(symV) <- spmean_long_sub$allotb

# simulate longevity evolution across phylogeny to visualise
contmap <- phytools::contMap(allotb, spmean_long_subV, plot = FALSE)
# re-colour contmap with custom scale
contmap <- phytools::setMap(contmap, my_colours$longevity)

# make sure discrete character is in the order of tree
symV <- symV[contmap$tree$tip.label]
# set factor colours
cols <- setNames(c(my_colours$symmetry[2], my_colours$symmetry[1]), c("0", "1"))


#* tall plot ----

# dummy plot to get locations to draw tip labels coloured by symmetry
plot(contmap, fsize = c(0.5, 0.7))
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# export figure
pdf(file = "figures/contmap_spmeanlongevity.pdf", width = 20, height = 120)
plot(contmap, legend = 0.7*max(nodeHeights(allotb)), sig = 1, 
     lwd = 4, outline = FALSE, ftype = "off", #type = "fan",
     xlim = lastPP$x.lim + 20, #ylim = lastPP$y.lim,
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
    text(x = xpos+3, y = tip_pos$y0[[i]] + ((tip_pos$y1[[i]] - tip_pos$y0[[i]])/2), 
         paste(tip_pos$family[[i]], sep = ""), adj = 0, srt = 0, cex = 0.4)
  }
  rm(tip_pos)
}

family_labels(xpos = 159)

# function to work out highest and lowest tip numbers for each order,
# then loop through these to draw segment and text labels for each order
order_labels <- function(xpos = 200){
  tip_pos <- spmean_long_sub %>%
    dplyr::group_by(order) %>%
    dplyr::summarise(y0 = min(position), y1 = max(position))
  for (i in 1:nrow(tip_pos)){
    segments(x0 = xpos, y0 = tip_pos$y0[[i]], x1 = xpos, y1 = tip_pos$y1[[i]], lwd = 2)
    text(x = xpos + 3, y = tip_pos$y0[[i]] + ((tip_pos$y1[[i]] - tip_pos$y0[[i]])/2), 
         paste(tip_pos$order[[i]], sep = ""), adj = 0, srt = 0, cex = 0.4)
  }
  rm(tip_pos)
}

order_labels(xpos = 170)

# TO DO - label clades??!!!

# insert legend
legend(x = "bottomright", legend = c("actinomorphic", "zygomorphic"), bg = "white",
       fill = cols, cex = 3, pt.lwd = 0.001, bty = "n",
       title = "Floral symmetry")
dev.off()

rm(lastPP, family_labels, order_labels, i)

#* circular plot ----

# build fan style plot for main figure
pdf(file = "figures/allotb_longevity_contMap_fan.pdf", width = 15, height = 15)

plot(contmap, type = "fan", legend = FALSE, lwd = 4, outline = FALSE, 
     ftype = "off", xlim = c(-185, 150), rotate.tree = 180)

# label Cretaceous period, 139 to 66 mya shown here
plotrix::draw.circle(0, 0, radius = max(nodeHeights(allotb)) - 66, 
                     col = "#f8f6f7", lty = 0)

# plot contMap again
par(new = TRUE) # hack to force below to plot on top of above 
plot(contmap, type = "fan", legend = FALSE, lwd = 4, outline = FALSE, 
     ftype = "off", xlim = c(-185, 150), rotate.tree = 180, add = TRUE)

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
  offset = 2
)

# add flower symmetry points
points(xx_yy$xx,
       xx_yy$yy,
       pch = 15, cex = 1,
       col = cols[symV[allotb$tip.label]])

legend(x = 130, y = 150, legend = c("actinomorphic", "zygomorphic"), col = cols, 
       bty = "n", cex = 0.8, title = "Flower symmetry", pch = 15)

# greyscale clade labelling using custom arclabel function for circular (fan) tree

cladelabels_fan <- function(offset = 1){ # use offset argument to move labels closer (<1) or further away (>1) from tree
  source("scripts/functions/arclabel.R") # get arclabel function
  arclabel(text = "ANA", tips = c(1432, 1433),
           lwd = 20, cex = 1.6, col = "#bdbdbd",
           ln.offset = offset + .07, lab.offset = offset + .11,
           orientation = "perpendicular",
           lend = "butt")
  arclabel(text = "Magnoliids", tips = c(1388, 1431), 
           lwd = 20, cex = 1.6, col = "#636363",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Monocots", tips = c(1135, 1387), 
           lwd = 20, cex = 1.6, col = "#bdbdbd",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Commelinids", tips = c(1135, 1229), 
           lwd = 15, cex = 1.6, col = "#636363",
           ln.offset = offset + .06, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Eudicots", tips = c(1, 1134),
           lwd = 20, cex = 1.6, col = "#636363",
           ln.offset = offset + .07, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Rosids", tips = c(1, 376),
           lwd = 15, cex = 1.6, col = "#bdbdbd",
           ln.offset = offset + .06, lab.offset = offset + .11,
           lend = "butt")
  arclabel(text = "Asterids", tips = c(399, 956), 
           lwd = 15, cex = 1.6, col = "#bdbdbd",
           ln.offset = offset + .06, lab.offset = offset + .11,
           lend = "butt")
}

cladelabels_fan(offset = 1.05)

dev.off()

# TO DO
# - label clades
# - offset or enlarge tip points?


rm(symV, spmean_long_sub, lastPP, contmap, spmean_long_subV, i, cols)


