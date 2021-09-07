
# checking whether increases in power for some $\tau_j$ values is just due to
# them coming from simulations with higher average $\tau_j$!

# Answer: yes! the higher tau_j values do, on average, come from simulation runs that have higher tau_bars...
#  - but if we remove the site itself, the average of other sites remains the same

library( arm )   # note: loads MASS package
library( blkvar )

require(tidyverse)
require(glue)
require(tictoc)


#####
# Load simulation results
#####

FNAME <- "sim_results_ICC_fixed"

fname <- glue("results/", FNAME, ".csv")
hits <- read_csv(fname) %>%
  mutate(n_bar = as.factor(n_bar),
         ICC = as.factor(ICC),
         tau = as.factor(tau))


#####
# checking for issues
#####

# compute tau-bar for each simulation
#  - then match up tau_j values with tau_bar values
#  - see, for each n_bar, ICC, tau grouping AND tau_j value, what the average tau_bar is

# compute tau_bar for each simulation run
tau_bar_key <- hits %>%
  group_by(n_bar, J, ICC, tau, tx_var, runID) %>%
  summarize(tau_bar = mean(ATE))
ggplot(tau_bar_key) +
  geom_histogram(aes(x=tau_bar)) +
  facet_wrap(~tau) +
  geom_vline(aes(xintercept=as.numeric(as.character(tau))), color="red")

# match up each true tau_j with its simulation's tau_bar
hits %>%
  left_join(tau_bar_key, by=c("n_bar", "J", "ICC", "tau", "tx_var", "runID")) %>%
  group_by(n_bar, J, ICC, tau, tx_var, ATE) %>%
  summarize(mn_tau_bar = mean(tau_bar)) %>%
  filter(tx_var == 0.3) %>%
  # filter(n_bar == 25, ICC == 0, tx_var == 0.3) %>%
  ggplot(aes(x=ATE, y=mn_tau_bar)) +
  geom_point() +
  facet_wrap(~tau)




#####
# checking for issues, v2 (using tau-bar-minus-j)
#####

# compute tau-bar-MINUS-J for each simulation
#  - then match up tau_j values with tau_bar values
#  - see, for each n_bar, ICC, tau grouping AND tau_j value, what the average tau_bar_MINUS_J is

# compute tau_bar_mj for each SITE within each simulation run
helper_fun <- function(ATEs, sids) {
  map_dbl(sids, function(x) mean(ATEs[-x]))
}

tau_bar_mj_key <- hits %>%
  filter(n_bar == 100, ICC == 0.9, tau == 0.01) %>%
  group_by(n_bar, J, ICC, tau, tx_var, runID) %>%
  mutate(tau_bar_mj = helper_fun(ATE, sid))

lm(tau_bar_mj ~ ATE, data=tau_bar_mj_key) %>% summary

tau_bar_mj_key %>%
  ggplot(aes(x=ATE, y=tau_bar_mj)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept=0)

# we see that the most extreme values slope downward,
# but the average values don't actually have a slope

# also, the slope is from extreme simulations!
#  - as we remove higher and higher tau_j values, the tau_bar_mj value decreases lol
#  - EACH SIM PRODUCES A DOWNWARD-SLOPING LINE
#  - the line isn't exactly straight because...?



