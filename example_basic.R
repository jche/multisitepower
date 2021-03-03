
# basic example: 1-site RCT

# Power: P(detect effect | true nonzero effect)
# MDE:   if I want 80% power, how big does the effect have to be?

# 4 pieces of any power analysis => What to do for MDE:
#  1) effect size (standardized) => find this
#  2) sample size:               => vary this
#  3) significance level         => 0.05
#  4) power level                => 0.80

require(tidyverse)

# f: (# obs, effect size) => (reject null? T/F)
one_sim <- function(n, tau) {
  # generate data: constant treatment effect
  #  - note: sd=1, so TAU is already in Cohen's D units
  dat <- tibble(
    y0 = rnorm(n),
    tx = sample(0:1, n, replace=T),
    yobs = y0 + tau*tx
  )
  
  # compute t test
  lm <- lm(yobs ~ tx, data=dat)
  summary(lm)$coefficients["tx",4] < 0.05
}

# f: (# obs, effect size) => (power)
power_sim <- function(n, tau, NUMSIM = 250) {
  cat(glue("Working on n = {n}, tau = {tau}

"))
  (1:NUMSIM) %>%
    map_lgl(~one_sim(n, tau)) %>%
    mean()
}

# run power simulation
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


# the real power analysis has many more knobs for the DGP (this one only has n and tau):
#  - J, site.size (or R^2?), ICC (or variances?)
#
# Q: what is an "effect size" when we have an effect distribution?
#  - each site has an effect size, but detection of an effect for each site
#    also depends on all of the other sites!
#  - so just having averaging over all sites with the same "effect size" (across simulations)
#    doesn't tell the whole story!
# 
