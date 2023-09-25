library(tidyverse) #install.packages("tidyverse")
library(ggpubr) #install.packages("ggpubr")
library(rtry) #install.packages("rtry")
library(TNRS) #install.packages("TNRS")
library(ape) #install.packages("ape")
library(phylolm) #install.packages("phylolm")
library(phytools) #install.packages("phytools")
library(nlme) #install.packages("nlme")

# function to cache pre-prepared R data. If RDS already in cache will read data
source("scripts/functions/cache_csv.R")

# My rough plan: 
#~~ -pull together symmetry data from a few published sources as well as the ~~
#~~  PROTEUS database and some previous work of mine ~~
#~~ -align the taxonomy of symmetry data and Marcos' longevity data to WCVP ~~
#~~ -score symmetry data for any missing species ~~
# -match longevity/symmetry data to the Smith and Brown (2018) angiosperm phylogeny
# -check for any sampling bias of floral longevity and SUBSAMPLE
# -run tests for floral longevity and symmetry phylogenetic signal and evolutionary rates
# -run phylogenetic logistic regression to see whether floral longevity is longer 
#  or shorter in actinomorphic vs. zygomorphic flowers
# -run (phylogenetic?) GLM with site as a fixed effect for the community studies 
#  that you sent me already, to see how local longevity/symmetry relationships 
#  compare to the global relationship

# read in data from different sources
source("scripts/prepdata/final_data.R")

# NOW TO DO SOME PRELIMINARY ANALYSIS!!!
# first just compare mean longevity by symmetry
source("scripts/analysis/compare_longevity.R")

# would be interesting to see how phylogeny changes this - does the relationship
# get weaker or stronger when evolutionary history is considered?
source("scripts/analysis/phylogenetics.R")

