# script to save consistent colours for plotting

my_colours <- list()

# continuous colour scale for floral longevity, colourblind safe
my_colours$longevity <- c("#feedde", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603")

# two colours that will stand out against this longevity scale for symmetry
my_colours$symmetry <- c("#977A9F", "#FEEF8B")
names(my_colours$symmetry) <- c("zygomorphic", "actinomorphic")
