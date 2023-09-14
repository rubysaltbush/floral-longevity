# script to generate phylogeny for taxa from Smith and Brown (2018) tree
# and run some basic phylogenetic analyses to see if floral longevity and
# floral symmetry are co-evolving traits

# TO INVESTIGATE
# - randomly sample one species per genus. Same results?
# - use shorter Smith and Brown tree (GBOTB.tre with 79,881 tips)
# - use Smith and Brown without Qian and Jin?? (ALLOTB.tre with 353,185 tips)

# read in short Smith and Brown tree (GBOTB.tre with 79,881 tips)
gbotb <- ape::read.tree("data_input/GBOTB.tre")
# read in long Smith and Brown tree (ALLOTB.tre with 353,185 tips)
allotb <- ape::read.tree("data_input/ALLOTB.tre")

# read in taxonomic name matching key
source("scripts/prepdata/phylo_names_match.R")

plot(tree_TPL$scenario.3, type = "fan", show.tip.label = FALSE)
# it's a tree! who knows if it's right! crazy.
# Hervé unimpressed by number of polytomies

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

# drop missing data tips from tree NO MISSING DATA AS ALREADY FILTERED
# lose 0 tips at this stage
to_drop <- pgls %>%
  dplyr::filter(is.na(spmean_long_days)|!(sym_species %in% c("actinomorphic", "zygomorphic")))
tree_nomissing <- ape::drop.tip(tree_TPL$scenario.3, to_drop$species)
plot(tree_nomissing, type = "fan", show.tip.label = FALSE)
rm(to_drop)
# and remove missing data taxa from morphological data
pgls <- pgls %>%
  dplyr::filter(is.na(spmean_long_days)|sym_species %in% c("actinomorphic", "zygomorphic"))
# 1452 species obs remain
# reorder data so species order matches order of tips in tree
pgls <- as.data.frame(tree_nomissing$tip.label) %>%
  dplyr::left_join(pgls, by = c("tree_nomissing$tip.label" = "species"))
pgls$sym_species <- forcats::as_factor(pgls$sym_species)

# taxon_name to row names
rownames(pgls) <- pgls[,1]
pgls[,1] <- NULL

# # redefine as 0 and 1
pgls$sym_species <- gsub("zygomorphic", "1", pgls$sym_species)
pgls$sym_species <- gsub("actinomorphic", "0", pgls$sym_species)
table(pgls$sym_species)
# 985 actinomorphic taxa to 467 zygomorphic taxa

# double check distribution of continuous variables
plot(pgls$spmean_long_days) # some outlying high values
hist(pgls$spmean_long_days) # few outlying high values
# generally left-biased distribution, think okay for now?

# quick boxplot of this data subset to compare means between symmetry
boxplot(spmean_long_days ~ sym_species, data = pgls)
# marginally longer longevity for zygomorphy

#### run PGLS ####

PGLS_symlong <- nlme::gls(spmean_long_days ~ sym_species, 
                          correlation = ape::corBrownian(phy = tree_nomissing),
                          data = pgls, method = "ML")
anova(PGLS_symlong)

coef(PGLS_symlong)
# Warning message:
#   In Initialize.corPhyl(X[[i]], ...) :
#   No covariate specified, species will be taken as ordered in the data frame. 
#   To avoid this message, specify a covariate containing the species names with the 'form' argument.
# will have to work out how to interpret this, with Luke Harmon's book maybe??

#### run model ####
# below adapted from Joly and Schoen (2021)
# works without polymorphic or missing data
# Model fit with Ives and Garlan optimisation

PGLogS_symlong <- phylolm::phyloglm(sym_species ~ spmean_long_days, 
                                  data = pgls, 
                                  phy = tree_nomissing,
                                  method = "logistic_IG10", 
                                  boot = 500) # wow takes a long time with more bootstraps
# Warning message:
# In phylolm::phyloglm(sym_species ~ spmean_long_days, data = pgls,  :
                       #phyloglm failed to converge.
# even with 500 bootstraps! yikes.

summary(PGLogS_symlong)

# p = 0.0003 !

rm(PGLS_symlong, for_phylo)

#### PHYLOGENETIC ANOVA ####

sym <- pgls$sym_species
names(sym) <- rownames(pgls)

long <- pgls$spmean_long_days
names(long) <- rownames(pgls)

anova <- phytools::phylANOVA(tree_nomissing, sym, long)

print(anova)

anova$Pf

# hmm Pf = 0.375, much much higher than phylogenetic logistic regression, wonder why?

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
