# script to save consistent colours for plotting

my_colours <- list()

# continuous colour scale for floral longevity, colourblind safe
my_colours$longevity <- c("#648fff", "#ffffda", "#ffb000")

# two colours that will stand out against this longevity scale for symmetry
my_colours$symmetry <- c("#977A9F", "#fa9fb5")
names(my_colours$symmetry) <- c("zygomorphic", "actinomorphic")

