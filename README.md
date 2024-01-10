# floral-longevity

Data and analysis code for the paper *Zygomorphic flowers last longer: the evolution of floral symmetry and floral longevity*

Authors: R.E. Stephens (1,2), R.V. Gallagher (1,3), M. Méndez (4), H. Sauquet (2,5)
+ Corresponding author: Ruby E. Stephens, stephenseruby@gmail.com 

Author addresses:

1. School of Natural Sciences, Macquarie University, Ryde, NSW, Australia 
2. Royal Botanic Gardens and Domain Trust, Sydney, NSW, Australia
3. Hawkesbury Institute for the Environment, University of Western Sydney, Richmond, NSW, Australia
4. Area of Biodiversity and Conservation, Universidad Rey Juan Carlos, Madrid, Spain
5. Ecology & Evolution Research Centre, University of New South Wales, Sydney, NSW, Australia


## Re-running analysis

Using RStudio open [main.R](https://github.com/rubysaltbush/floral-longevity/blob/main/main.R) 
and run scripts in order given in this main script.

Necessary packages can be installed by running [install_dependencies.R](https://github.com/rubysaltbush/floral-longevity/blob/main/scripts/install_dependencies.R).
Check package versions mentioned in this script are consistent with your installed versions. 
Originally run in R version 4.3.0

The code caches several steps of the data cleaning and results.

Feel free to re-use scripts and functions for your own analyses, e.g.

[arclabel.R](https://github.com/rubysaltbush/floral-longevity/blob/main/scripts/functions/arclabel.R) adds custom labels to circular phylogenies by providing the tips to draw the label between

