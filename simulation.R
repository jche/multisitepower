
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
NUMSIM <- 10
df_sim <- expand_grid(
  # nbar   = seq(25, 300, by=25),
  nbar   = c(10, 25, 100, 300),
  J      = c(10, 25, 100, 300),
  ICC    = c(0.2),
  tau    = c(0, 0.2, 0.5),
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
                             variable.n = T,   # vary site sizes!
                             NUMSIM = NUMSIM,
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
    mutate(data = list(power_sim(nbar, J, tau, ICC, tx_var, 
                                 variable.n = T,   # vary site sizes!
                                 NUMSIM = NUMSIM, 
                                 WRITE_CSV = F)))
  # aggregate site-level results
  sim_sites <- df_sim %>%
    unnest(cols = data) %>%
    filter(names(data) == "sites") %>%
    unnest(cols = data)
  # aggregate overall results
  sim_overall <- df_sim %>%
    unnest(cols = data) %>%
    filter(names(data) == "overall") %>%
    unnest(cols = data)
  toc()
}


#####
# save results
#####

# FNAME <- "sim_ex_tau02"
FNAME <- "final_sim"
UUID  <- substr(UUIDgenerate(), 25, 36)

# save site-level results
fname_sites <- glue("results/{FNAME}-{UUID}.csv")
if (file.exists(fname_sites)) {
  # ASSUMING that all sim settings are run equally
  max_runID <- read_csv(fname_sites) %>%
    pull(runID) %>%
    max()

  sim_sites %>%
    mutate(runID = runID + max_runID) %>%
    write_csv(fname_sites, append=T)
} else {
  write_csv(sim_sites, fname_sites)
}

# save overall results
fname_overall <- glue("results/{FNAME}_overall-{UUID}.csv")
if (file.exists(fname_overall)) {
  # ASSUMING that all sim settings are run equally
  max_runID <- read_csv(fname_overall) %>%
    pull(runID) %>%
    max()
  
  sim_overall %>%
    mutate(runID = runID + max_runID) %>%
    write_csv(fname_overall, append=T)
} else {
  write_csv(sim_overall, fname_overall)
}
