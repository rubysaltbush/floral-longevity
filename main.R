# uncomment below script to install all necessary packages
#source("scripts/install_dependencies.R")

library(tidyverse)
library(ggpubr)
library(rtry)
library(TNRS)
library(ape)
library(phytools)
library(nlme)
library(caper)

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_csv.R")

#### DATA ####

# custom colour scales
source("scripts/prepdata/my_colours.R")

# read in data from different sources
source("scripts/prepdata/final_data.R")

#### ANALYSES ####

# first just compare mean longevity by symmetry and describe raw data
source("scripts/analysis/compare_longevity.R")

# main analysis script, considering phylogenetic relatedness, also with 
# figures to show longevity and symmetry on phylogeny
source("scripts/analysis/phylogenetics.R")

