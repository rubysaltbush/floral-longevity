# uncomment below script to install all necessary packages
#source("scripts/install_dependencies.R")

library(tidyverse)
library(ggpubr)
library(rtry)
library(TNRS)
library(ape)
library(phylolm) # might not need this? only for phylogenetic logistic regression I think
library(phytools)
library(nlme)

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_csv.R")

# My rough plan: 
#~~ -pull together symmetry data from a few published sources as well as the ~~
#~~  PROTEUS database and some previous work of mine ~~
#~~ -align the taxonomy of symmetry data and Marcos' longevity data to WCVP ~~
#~~ -score symmetry data for any missing species ~~
#~~ -match longevity/symmetry data to the Smith and Brown (2018) angiosperm phylogeny ~~
#~~ -check for any sampling bias of floral longevity and SUBSAMPLE ~~
#~~ -run phylogenetic least squares regression to see whether floral longevity is longer 
#  or shorter in actinomorphic vs. zygomorphic flowers ~~
# -run tests for floral longevity and symmetry phylogenetic signal and evolutionary rates?
# -run (phylogenetic?) GLM with site as a fixed effect for the community studies 
#  that you sent me already, to see how local longevity/symmetry relationships 
#  compare to the global relationship

# could get 6 more species from Ke, Y., Zhang, F.-P., Zhang, Y.-B., Li, W., Wang, Q., Yang, D., Zhang, J.-L., & Cao, K.-F. (2023). Convergent relationships between flower economics and hydraulic traits across aquatic and terrestrial herbaceous plants. Plant Diversity. https://doi.org/10.1016/j.pld.2023.01.006
# probs not worth hassle??? any new families? yep, Nymphaeaceae and Menyanthaceae, but no published data with paper gah

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

