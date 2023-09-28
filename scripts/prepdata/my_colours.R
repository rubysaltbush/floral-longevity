# script to save consistent colours for plotting

my_colours <- list()

# continuous colour scale for floral longevity, colourblind safe
my_colours$longevity <- c("#ffffcc", "#c2e699", "#78c679", "#31a354", "#006837")

# two colours that will stand out against this longevity scale for symmetry
my_colours$symmetry <- c("#977A9F", "#FEEF8B")
names(my_colours$symmetry) <- c("zygomorphic", "actinomorphic")

