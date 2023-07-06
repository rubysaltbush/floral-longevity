library(tidyverse) #install.packages("tidyverse")
library(ggpubr) #install.packages("ggpubr")
library(rtry) #install.packages("rtry")
library(TNRS) #install.packages("TNRS")
library(ape) #install.packages("ape")
library(V.PhyloMaker2) #library(devtools) devtools::install_github("jinyizju/V.PhyloMaker2")
library(phylolm) #install.packages("phylolm")
library(phytools) #install.packages("phytools")

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_RDS.R")

# My rough plan: 
#~~ -pull together symmetry data from a few published sources as well as the~~
#~~  PROTEUS database and some previous work of mine~~
#~~ -align the taxonomy of symmetry data and Marcos' longevity data to the World~~
#~~  Checklist of Vascular Plants~~
# -score symmetry data for any missing species
#~~ -match longevity/symmetry data to the Smith and Brown (2018) angiosperm phylogeny
# -check for any sampling bias of floral longevity
# -run tests for floral longevity phylogenetic signal and evolutionary rate
# -run phylogenetic logistic regression to see whether floral longevity is longer 
#  or shorter in actinomorphic vs. zygomorphic flowers
# -run (phylogenetic?) GLM with site as a fixed effect for the community studies 
#  that you sent me already, to see how local longevity/symmetry relationships 
#  compare to the global relationship

# read in data from different sources
source("scripts/prepdata/final_data.R")

# NOW TO DO SOME PRELIMINARY ANALYSIS!!!
# first just compare mean longevity by symmetry
source("analysis/compare_longevity.R")

# would be interesting to see how phylogeny changes this - does the relationship
# get weaker or stronger when evlutionary history is considered?
source("analysis/phylogenetics.R")