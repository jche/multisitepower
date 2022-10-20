
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
require(uuid)

# source required functions
source("simulation_functions.R")

# turn off summarize messages
options(dplyr.summarise.inform = FALSE)

# set simulation seed
# set.seed(02141)

# initialize parallel backend
#  - good tutorial: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
PARALLEL <- F
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
NUMSIM <- 500   # J * number of simulations = NUMSIM; must be greater than max(J)
# df_sim <- expand_grid(
#   nbar   = c(50,100,250),
#   J      = c(5, 10, 25, 50, 100),
#   ICC    = c(0.2),
#   tau    = c(0, 0.2, 0.5),
#   tx_sd  = c(0.1, 0.2, 0.3)
# )
df_sim <- expand_grid(
  nbar   = seq(5,200,5),
  # nbar   = c(25, 50),
  J      = c(25),
  ICC    = c(0.2),
  tau    = c(0.2),
  tx_sd  = c(0.3)
)

if (PARALLEL) {
  tic()
  sim <- foreach(nbar= df_sim$nbar,
                 J   = df_sim$J,
                 ICC = df_sim$ICC,
                 tau = df_sim$tau,
                 tx_sd = df_sim$tx_sd,
                 .combine = rbind,
                 .packages=c("blkvar", "tidyverse", "glue", "arm", "rstan")) %dopar% {
                   
                   # load functions
                   source("simulation_functions.R")
                   
                   # power_sim
                   power_sim(nbar=nbar, J=J, ICC=ICC, tau=tau, tx_sd=tx_sd,
                             variable.n = T,   # vary site sizes!
                             NUMSIM = NUMSIM/J,
                             WRITE_CSV = F) %>%
                      mutate(nbar = nbar, J=J, ICC=ICC, tau=tau, tx_sd=tx_sd,
                             .before = 1)
                 }
  stopCluster(cl = workers)
  toc()
  
  # aggregate site-level results
  #  - NOTE: this is patchy, it's only run for the single-site estimates
  #    and NOT for the overall estimates
  sim_sites <- sim %>%
    unnest(cols = data)
  
  # 435.58 seconds for small sim, 10 runs x 16 factors
  # 1774.91 seconds for full sim, 10 runs x 81 factors
} else {
  # run simulation: store power_sim() results in df_sim as list column
  tic()
  sim <- df_sim %>%
    rowwise() %>%
    mutate(data = list(power_sim(nbar, J, tau, ICC, tx_sd, 
                                 variable.n = T,   # vary site sizes!
                                 NUMSIM = NUMSIM/J, 
                                 WRITE_CSV = F))) # %>%
    # unnest(cols = data)
  
  # aggregate site-level results
  sim_sites <- sim %>%
    unnest(cols = data) %>%
    filter(names(data) == "sites") %>%
    unnest(cols = data)
  # # aggregate overall results
  # sim_overall <- sim %>%
  #   unnest(cols = data) %>%
  #   filter(names(data) == "overall") %>%
  #   unnest(cols = data)
  # toc()
}


#####
# save results
#####

FNAME <- "sree_sims"
UUID  <- substr(UUIDgenerate(), 25, 36)

# save site-level results
fname <- glue("results_sree/{FNAME}-{UUID}.csv")
# fname <- glue("revised_results/{FNAME}.csv")
if (file.exists(fname)) {
  # ASSUMING that all sim settings are run equally
  max_runID <- read_csv(fname) %>%
    pull(runID) %>%
    max()

  sim_sites %>%
    mutate(runID = runID + max_runID) %>%
    write_csv(fname, append=T)
} else {
  write_csv(sim_sites, fname)
}

# # save overall results
# fname_overall <- glue("revised_results/{FNAME}_overall-{UUID}.csv")
# # fname_overall <- glue("revised_results/{FNAME}_overall.csv")
# if (file.exists(fname_overall)) {
#   # ASSUMING that all sim settings are run equally
#   max_runID <- read_csv(fname_overall) %>%
#     pull(runID) %>%
#     max()
# 
#   sim_overall %>%
#     mutate(runID = runID + max_runID) %>%
#     write_csv(fname_overall, append=T)
# } else {
#   write_csv(sim_overall, fname_overall)
# }
