# uncomment below script to install all necessary packages
#source("scripts/install_dependencies.R")

library(tidyverse)
library(ggpubr)
library(rtry)
library(TNRS)
library(ape)
library(phytools)
library(nlme)

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_csv.R")

# TO DO

# -run tests for floral longevity and symmetry phylogenetic signal and evolutionary rates?

#### DATA ####

# read in data from different sources
source("scripts/prepdata/final_data.R")

# custom colour scales
source("scripts/prepdata/my_colours.R")

#### ANALYSES ####

# first just compare mean longevity by symmetry and describe raw data
source("scripts/analysis/compare_longevity.R")

# would be interesting to see how phylogeny changes this - does the relationship
# get weaker or stronger when evolutionary history is considered?
source("scripts/analysis/phylogenetics.R")

