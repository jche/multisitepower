
#####
# Power analysis for individual site effects in multisite trial,
# comparing MLMs to just using each single site
# 
# uses blkvar package: devtools::install_github("https://github.com/lmiratrix/blkvar")
#####

# technical packages
library( arm )   # bayesian-type functions for lmer, note this loads MASS package
# devtools::install_github("lmiratrix/blkvar")
library( blkvar )
library(rstan)   # bayesian models!

# admin packages
require(tidyverse)
require(glue)
require(tictoc)

# source required functions
source("simulation_functions.R")

# turn off summarize messages
options(dplyr.summarise.inform = FALSE)

# set simulation seed
set.seed(02141)

# initialize parallel backend
#  - good tutorial: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
PARALLEL <- TRUE
if (PARALLEL) {
  require(foreach)
  require(doParallel)
  workers <- makeCluster(7, type="PSOCK")
  print(workers)
  
  registerDoParallel(cl=workers)
  getDoParRegistered()
  getDoParWorkers()
}


#####
# run power simulation
#####

# simulation settings
df_sim <- expand_grid(
  nbar= c(25, 50, 75),
  J   = c(25, 50, 75),
  ICC = c(0, 0.3, 0.6),
  tau = c(0.01, 0.2, 0.5),
  tx_var = c(0.3)
)

if (PARALLEL) {
  tic()
  sim <- foreach(nbar= df_sim$nbar,
                 J   = df_sim$J,
                 ICC = df_sim$ICC,
                 tau = df_sim$tau,
                 tx_var = df_sim$tx_var,
                 .combine = rbind,
                 .packages=c("blkvar", "tidyverse", "glue", "arm", "rstan")) %dopar% {
                   
                   # load functions
                   source("simulation_functions.R")
                   
                   # power_sim
                   power_sim(nbar=nbar, J=J, ICC=ICC, tau=tau, tx_var=tx_var,
                             NUMSIM = 10,
                             WRITE_CSV = F) %>%
                      mutate(nbar = nbar, J=J, ICC=ICC, tau=tau, tx_var=tx_var,
                             .before = 1)
                 }
  stopCluster(cl = workers)
  toc()
  
  # 435.58 seconds for small sim, 10 runs x 16 factors
  # 1774.91 seconds for full sim, 10 runs x 81 factors
} else {
  # run simulation: store power_sim() results in df_sim as list column
  tic()
  df_sim <- df_sim %>%
    rowwise() %>%
    mutate(data = list(power_sim(nbar, J, tau, ICC, tx_var, NUMSIM = 10, 
                                 WRITE_CSV = F)))
  sim <- df_sim %>%
    unnest(cols = data)
  toc()
  
  # 956.9 seconds, 10 runs x 16 factors
}


#####
# save results
#####

FNAME <- "sim_results_bayes3"
fname <- glue("results/", FNAME, ".csv")

if (file.exists(fname)) {
  # ASSUMING that all sim settings are run equally
  max_runID <- read_csv(fname) %>%
    pull(runID) %>%
    max()

  sim %>%
    mutate(runID = runID + max_runID) %>%
    write_csv(fname, append=T)
} else {
  write_csv(hits, fname)
}

