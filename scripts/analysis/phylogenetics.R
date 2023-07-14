# script to generate phylogeny for taxa from Smith and Brown (2018) tree
# and run some basic phylogenetic analyses to see if floral longevity and
# floral symmetry are co-evolving traits

# TO INVESTIGATE
# - randomly sample one species per genus. Same results?
# - use shortwe Smith and Brown tree (GBOTB.tre with 79,881 tips)
# - use Smith and Brown without Qian and Jin?? (ALLOTB.tre with 353,185 tips)

# first use V.PhyloMaker2 package to get phylogeny from Smith and Brown GBOTB.extended tree
# prep data
for_phylo <- sym_long %>%
  dplyr::select(species = Accepted_name, genus = Accepted_genus, family = Accepted_family) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(species))

# compare 3 dif name references, leave scenario as default for now (ask Hervé later)
# ultimately best might be scenario 2 which randomizes, with high randomisation?
tree_TPL <- V.PhyloMaker2::phylo.maker(for_phylo)
# [1] "Taxonomic classification not consistent between sp.list and tree."
# genus family_in_sp.list family_in_tree
# 250  Emmotum       Icacinaceae Metteniusaceae
# 726 Viburnum       Viburnaceae      Adoxaceae
# [1] NA

#tree_LCVP <- V.PhyloMaker2::phylo.maker(for_phylo, tree = GBOTB.extended.LCVP)
# [1] "Taxonomic classification not consistent between sp.list and tree."
# genus family_in_sp.list family_in_tree
# 249  Emmotum       Icacinaceae Metteniusaceae
# 721 Viburnum       Viburnaceae      Adoxaceae
# [1] "Note: 1 taxa fail to be binded to the tree,"
# [1] NA
# Error in 1:n : argument of length 0
# doesn't work! presumably because of above mysterious error

#tree_WP <- V.PhyloMaker2::phylo.maker(for_phylo, tree = GBOTB.extended.WP)
# [1] "Taxonomic classification not consistent between sp.list and tree."
# genus family_in_sp.list family_in_tree
# 248  Emmotum       Icacinaceae Metteniusaceae
# 719 Viburnum       Viburnaceae      Adoxaceae
# [1] "Note: 1 taxa fail to be binded to the tree,"
# [1] NA
# Error in 1:(a1 - 1) : NA/NaN argument
# ALSO doesn't work! different mysterious error

plot(tree_TPL$scenario.3, type = "fan", show.tip.label = FALSE)
# it's a tree! who knows if it's right! crazy.

table(tree_TPL$species.list$status)
# not sure exactly what bind and prune mean? Maybe bind means it just bound this
# straight to an available tip, whereas prune means it had to use genus info???

#### PHYLOGENETIC LOGISTIC REGRESSION ####
# subset data to variables of interest for phylogenetic logistic regression
# take mean per species for longevity
pgls <- sym_long %>%
  dplyr::group_by(species = Accepted_name, sym_species) %>%
  dplyr::summarise(spmean_long_days = mean(mean_long_days)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(species)) %>%
  as.data.frame()
# add underscore to match tip labels
pgls$species <- gsub(" ", "_", pgls$species)

# drop missing data tips from tree
# lose 693 tips at this stage
to_drop <- pgls %>%
  dplyr::filter(is.na(spmean_long_days)|!(sym_species %in% c("actinomorphic", "zygomorphic")))
tree_nomissing <- ape::drop.tip(tree_TPL$scenario.3, to_drop$species)
plot(tree_nomissing, type = "fan", show.tip.label = FALSE)
rm(to_drop)
# and remove missing data taxa from morphological data
pgls <- pgls %>%
  dplyr::filter(is.na(spmean_long_days)|sym_species %in% c("actinomorphic", "zygomorphic"))
# 757 species obs remain

# taxon_name to row names
rownames(pgls) <- pgls[,1]
pgls[,1] <- NULL

# # redefine as 0 and 1
pgls$sym_species <- gsub("zygomorphic", "1", pgls$sym_species)
pgls$sym_species <- gsub("actinomorphic", "0", pgls$sym_species)
table(pgls$sym_species)
# 499 actinomorphic taxa to 258 zygomorphic taxa

# double check distribution of continuous variables
plot(pgls$spmean_long_days) # some outlying high values
hist(pgls$spmean_long_days) # few outlying high values
# generally left-biased distribution, think okay for now?

# quick boxplot of this data subset to compare means between symmetry
boxplot(spmean_long_days ~ sym_species, data = pgls)
# marginally longer longevity for zygomorphy but minimal really, doubt I'll see anythign in analysis

#### run model ####
# below adapted from Joly and Schoen (2021)
# works without polymorphic or missing data
# Model fit with Ives and Garlan optimisation

PGLS_symlong <- phylolm::phyloglm(sym_species ~ spmean_long_days, 
                                  data = pgls, 
                                  phy = tree_nomissing,
                                  method = "logistic_IG10", 
                                  boot = 100)

summary(PGLS_symlong)
# actually kind of close! p = 0.004, with preliminary messy data set
# ultimately probs want to run this as a phylogenetic t-test??

rm(PGLS_symlong, for_phylo)

#### PHYLOGENETIC ANOVA ####

sym <- pgls$sym_species
names(sym) <- rownames(pgls)

long <- pgls$spmean_long_days
names(long) <- rownames(pgls)

anova <- phytools::phylANOVA(tree_nomissing, sym, long)

print(anova)

anova$Pf

# hmm Pf = 0.393, much much higher than phylogenetic logistic regression, wonder why?

rm(anova, tree_nomissing, sym, long, pgls)

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
cols <- setNames(RColorBrewer::brewer.pal(n = length(levels(symV)), "Dark2"),
                 levels(symV))

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
