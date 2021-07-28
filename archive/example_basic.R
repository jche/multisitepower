
#####
# basic example: 1-site RCT
#
# Power: if there is a nonzero effect of size x, what is my chance of detecting it?
# MDE: if I want 80% power, how big does the effect have to be?
#
# 4 pieces of any power analysis => What to do for MDE:
#
#  1) effect size (standardized) => find this
#  2) sample size:               => vary this
#  3) significance level         => set to 0.05
#  4) power level                => set to 0.80
#####

require(tidyverse)
require(glue)

# function: (# obs, effect size) => (reject null? T/F)
one_sim <- function(n, tau) {
  # generate data: constant treatment effect
  #  - note: sd=1, so TAU is already in Cohen's D units
  dat <- tibble(
    y0 = rnorm(n),
    tx = sample(0:1, n, replace=T),
    yobs = y0 + tau*tx
  )
  
  # compute t-test
  lm <- lm(yobs ~ tx, data=dat)
  summary(lm)$coefficients["tx",4] < 0.05
}

# function: (# obs, effect size) => (power)
power_sim <- function(n, tau, NUMSIM = 250) {
  cat(glue("Working on n = {n}, tau = {tau}

"))
  (1:NUMSIM) %>%
    map_lgl(~one_sim(n, tau)) %>%
    mean()
}

# run power simulation:
#  - n = number of observations
#  - tau = true effect size
#  - (real power sim would have more knobs: J, site.size, ICC/variances, etc.)
df_sim <- expand_grid(
  n = c(25, 50, 75, 100),
  tau = c(0.01, 0.2, 0.5, 0.8)
)
df_sim_res <- df_sim %>%
  rowwise() %>%
  mutate(power = power_sim(n, tau))

# visualize results
ggplot(df_sim_res, aes(x=n, y=power, color=tau)) +
  geom_point() +
  geom_line(aes(group=tau)) +
  geom_hline(aes(yintercept=0.8), color="red", lty="dashed") +
  coord_cartesian(ylim=c(0,1))


# Note on "effect size" when we have an effect distribution:
#  - each site has a (random) effect size
#     - the simple idea is to aggregate "effect sizes" across simulations 
#       (within effect-size bins) to compute power for that effect-size bin
#
#  - wrinkle: detection of an effect for each site *technically* depends on all of 
#    the other sites in a given simulation (via the s.e. estimate)
#     - so just averaging all sites within the same effect-size bin (across simulations) 
#       doesn't tell the whole story, but it's close enough, I think




