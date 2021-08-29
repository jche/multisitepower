
#####
# Check CI coverage results
#####

library( arm )   # note: loads MASS package
library( blkvar )

require(tidyverse)
require(glue)
require(tictoc)


#####
# Load (small) simulation results
#####

# FNAME <- "sim_results_ICC_fixed"
FNAME <- "sim_results_fixed"   # added intercept variance

fname <- glue("results/", FNAME, ".csv")
hits <- read_csv(fname) %>%
  mutate(n_bar = as.factor(n_bar),
         ICC = as.factor(ICC),
         tau = as.factor(tau))

#####
# check EB coverage: averaged across tau_j ~ N(tau, omega)
#####

coverage_df <- hits %>%
  mutate(CI_low = ATE_hat - 1.96*SE,
         CI_high = ATE_hat + 1.96*SE,
         covered = CI_low <= ATE & ATE <= CI_high)

# EB coverage plot
# Want to see: 95% coverage across all settings
#  - actually see: coverage depending on informativeness of the data...
coverage_df %>%
  group_by(n_bar, J, ICC, tau, tx_var) %>%
  summarize(coverage = mean(covered)) %>%
  ggplot() +
  geom_point(aes(x=tau, y=coverage)) +
  geom_hline(yintercept=0.95, lty="dashed") +
  facet_grid(ICC ~ n_bar, labeller = label_both) +
  coord_cartesian(ylim = c(0.9, 1))

# Conditional coverage plot
# Want to see: curvature around true tau, with overcoverage near tau and undercoverage far out
#  - actually see: right shape, wrong height (no overcoverage)
coverage_df %>%
  mutate(tau = as.numeric(as.character(tau))) %>%
  filter(ICC == 0, ATE > tau-1, ATE < tau+1) %>%
  group_by(n_bar, J, ICC, tau, tx_var, ATE) %>%
  summarize(coverage = mean(covered)) %>%
  mutate(tau = as.factor(tau)) %>%
  ggplot(aes(color=tau)) +
  geom_line(aes(x=ATE, y=coverage, group=tau)) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  geom_vline(aes(xintercept = as.numeric(as.character(tau)))) +
  facet_grid(tau ~ n_bar, labeller = label_both)

# Weird distribution of hits plot, for illustrative purposes
# As information increases, distribution of misses approaches distribution of hits
#  - for low information settings, misses tend to be the extreme sites!
coverage_df %>%
  filter(ICC == 0) %>%
  # filter(ICC == 0.9) %>%
  ggplot() +
  geom_density(aes(x=ATE, color = covered, group=covered)) +
  geom_vline(aes(xintercept = as.numeric(as.character(tau)))) +
  facet_grid(n_bar~tau, labeller = label_both)


# ISSUE: our "fix" isn't fixing things... we're achieving 95% coverage in low-info cases,
#  but 100% coverage in high-info cases...?
#
#  - se.ranef -> 0 as site size increases
#  - se.fixef -> 0 as J increases
#
# - before: correct coverage for low-J, high nbar
#   and underoverage for low-J, low nbar: this is because...?
# - now, we get correct coverage for low-J, low nbar
#   and overcoverage for low-J, high nbar: this is because...?
#
# - expectation: overcoverage for sites close to tau, undercoverage for sites away from tau
#    - EB coverage SHOULD be consistently 95% for every setting, though!
# - reality: see the first phenomenon, DON'T SEE the second phenomenon...

# check FIRC paper, to see how they make site CIs: they use tau^2 * sqrt(1 - lambda_j)
#  for cross-site variance tau^2 and EB weight for site j lambda_j
#  - exactly estimate site: EB weight is 1 (no shrinkage), so SE^2 is super small
#  - no cross-site variance: overall ATE is a good estimate, so SE^2 is super small
#
# compared to us: when we exactly estimate the site, we still have high variance because of the standard error of our fixed-effect estimate (which is not quite an estimate of tau^2? or is it?)
#  - we're doing something additive, not multiplicative


