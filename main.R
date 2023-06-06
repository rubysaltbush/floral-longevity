library(tidyverse)
library(ggpubr)

# My rough plan: 
# -pull together symmetry data from a few published sources as well as the 
#  PROTEUS database and some previous work of mine 
# -align the taxonomy of symmetry data and Marcos' longevity data to the World 
#  Checklist of Vascular Plants
# -score symmetry data for any missing species
# -match longevity/symmetry data to the Smith and Brown (2018) angiosperm phylogeny
# -check for any sampling bias of floral longevity
# -run tests for floral longevity phylogenetic signal and evolutionary rate
# -run phylogenetic logistic regression to see whether floral longevity is longer 
#  or shorter in actinomorphic vs. zygomorphic flowers
# -run (phylogenetic?) GLM with site as a fixed effect for the community studies 
#  that you sent me already, to see how local longevity/symmetry relationships 
#  compare to the global relationship

# read in data from different sources
source("scripts/data_prep/final_data.R")
