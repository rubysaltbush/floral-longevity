# uncomment below script to install all necessary packages
#source("scripts/install_dependencies.R")

library(tidyverse)
library(ggpubr)
library(rtry)
library(TNRS)
library(ape)
#library(phylolm) # might not need this? only for phylogenetic logistic regression I think
library(phytools)
library(nlme)

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_csv.R")

# TO DO

# -run tests for floral longevity and symmetry phylogenetic signal and evolutionary rates?
# -run (phylogenetic?) GLM with site as a fixed effect for the community studies 
#  that you sent me already, to see how local longevity/symmetry relationships 
#  compare to the global relationship

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

