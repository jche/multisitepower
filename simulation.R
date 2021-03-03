
#####
# Power analysis for individual site effects in multisite trial
# 
# uses blkvar package: devtools::install_github("https://github.com/lmiratrix/blkvar")
#####

# load packages
require(plyr)
require(tidyverse)
require(blkvar)

library(uuid)
library(glue)
library(tictoc)

# setup
uid = str_sub(UUIDgenerate(), start= -12)
FILE_NAME = paste("simulation_results_paper/", re.distribution, "_", J, "_", site.size,"_",r.squared*100, "_", uid, sep = "")
scat("Saving to File '%s'\n", FILE_NAME)


# set simulation parameters: number of sites, R^2, etc.
J <- 50
site.size <- 50
ICC <- 0.30
r.squared <- 0



# run simulation

simulation_function <- function(r) {
  # generate dataset
  generate_multilevel_data()
  
  # fit hierarchical models
  fit_model()
  
  # evaluate model performance
  evaluate_model()
  
  # output results
  output_results()
  
}

# some sort of apply() function here


# write results to csv



