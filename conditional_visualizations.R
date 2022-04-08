
#####
# big idea: instead of conditioning on tau_j, let's try conditioning on tau-hat_j?
#####

require(tidyverse)
require(glue)
require(latex2exp)


#####
# load results
#####

tidy_results <- read_csv("tidy_results.csv")
tidy_results <- tidy_results %>%
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05))
  # mutate(ATEhat_bin = cut_interval(ATEhat, length=0.05) %>%
  #          as.character() %>%
  #          str_extract(., "\\d.*(?=,)") %>%
  #          as.numeric())

#####
# set fixed simulation factors
#####

# 5 potential covariate factors to visualize:
#  * ATE (true)
#  * method
#  * nbar
#  * tx_sd
#  * tau
#  * J

INCLUDE_OVERALL <- F   # include overall estimators?
TAU_FIXED <- 0.2
J_FIXED <- 100

# 5 potential outcomes to visualize:
#  * ATEhat
#  * SE
#  * MDES (one-sided, two-sided)
#  * power (one-sided, two-sided)
#  * coverage (one-sided, two-sided)

### prep summary dfs

MIN_SIZE <- 50

# prep power df
#  - NOTE: power is probability of TRUE rejection, so need to account for that here
power_df <- tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  group_by(nbar, tx_sd, tau, J, method, ATEhat_bin) %>%
  filter(ATE > 0) %>%
  summarize(n = n(),
            power = mean(reject_one)) %>%
  filter(n >= MIN_SIZE)

# prep coverage dfs
coverage2_df <- tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  group_by(nbar, tx_sd, tau, J, method, ATEhat_bin) %>%
  summarize(n = n(),
            coverage_two = mean(covered_two)) %>%
  filter(n >= MIN_SIZE)
coverage1_df <- tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  group_by(nbar, tx_sd, tau, J, method, ATEhat_bin) %>%
  summarize(n = n(),
            coverage_one = mean(covered_one)) %>%
  filter(n >= MIN_SIZE)

#####
# coverage plots
#####

# THINKING ABOUT A BIG IDEA: vertical average != horizontal average
#  - For each tau_j value, MLMs are biased toward tau
#  - For each tauhat_j value, single-site estimates are biased toward tau
#    - idea: "if I happen to estimate 0.4, but my population is mostly
#      around 0.2, then it's probably due to error, so in reality I'm
#      closer to 0.2."
#    - I.e., we need to do shrinkage here!

res_sm <- tidy_results %>%
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(J == J_FIXED, tau == TAU_FIXED, nbar == 250, tx_sd == 0.1) %>%
  group_by(nbar, J, ICC, tau, tx_sd, method, ATEhat_bin) %>%
  mutate(n = n()) %>%
  filter(n >= 50) %>%
  ungroup()

# # unaveraged scatterplot of ATEhat vs. ATE
# res_sm %>% 
#   sample_n(10000) %>%
#   ggplot(aes(x=ATE, y=ATEhat_bin, color=method)) +
#   geom_point() +
#   geom_abline() +
#   facet_grid( ~ method, labeller=label_both)

# vertical averages (per ATE, average ATEhat)
vert_avg <- res_sm %>% 
  group_by(method, ATE) %>%
  summarize(ATEhat_bin = mean(ATEhat_bin))

# horizontal averages(per ATEhat, average ATE)
hor_avg <- res_sm %>% 
  group_by(method, ATEhat_bin) %>%
  summarize(ATE = mean(ATE))

# plot everything at once
res_sm %>% 
  filter(method %in% c("bayesnorm", "single")) %>%
  sample_n(10000) %>%
  ggplot(aes(x=ATE, y=ATEhat_bin, color=method)) +
  geom_point() +
  geom_point(data=vert_avg %>% 
               filter(method %in% c("bayesnorm", "single")), 
             color="black") +
  geom_point(data=hor_avg %>% 
               filter(method %in% c("bayesnorm", "single")), 
             color="red") +
  geom_abline() +
  facet_grid( ~ method, labeller=label_both)


# THESE ARE THE OPPOSITE(ISH) OF THE OTHER COVERAGE PLOTS!
#  - single-site coverage is curved (when the estimate is extreme, it's probably wrong)
#  - MLM coverage is (almost) flat!
# MLM coverage dips down at tau, because shrinkage pulls everything toward tau,
#  meaning that sites with tau_j \neq tau are estimated as approximately tau
#  - this is extra bad when tx_sd is low, and we estimate zero cross-site variance
#  - degenerate models happen more with FIRC than with Bayes

# plot two-sided interval coverage vs. ATE size
coverage2_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  ggplot() +
  geom_line(aes(x=ATEhat_bin, y=coverage_two, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  geom_vline(xintercept=TAU_FIXED) +
  labs(title = "Conditional coverage rates (two-sided CI)",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Coverage",
       x = "Estimated site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/cond_coverage_2sided.png"), width=200, height=125, units="mm")

# plot one-sided interval coverage vs. ATE size
coverage1_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  ggplot() +
  geom_line(aes(x=ATEhat_bin, y=coverage_one, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Conditional coverage rates (one-sided CI)",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Coverage",
       x = "Estimated site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/cond_coverage_1sided.png"), width=200, height=125, units="mm")

# EB coverage is the same


#####
# checking power
#####

# Plot power vs. estimated ATE
#  - NOTE: power is probability of TRUE rejection, so need to account for that here
power_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  filter(method %in% c("firc2", "rirc2", "bayesnorm", "single")) %>%
  mutate(method = case_when(
    method == "firc2" ~ "FIRC",
    method == "rirc2" ~ "RIRC",
    method == "bayesnorm" ~ "Bayes",
    method == "single" ~ "Single")) %>%
  ggplot(aes(x = ATEhat_bin, y = power, col = method, group=method ) ) +
  facet_grid(tx_sd ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = 0.1, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 0.6)) +
  labs(title = glue("Power (\u03B1 = 0.1) vs. true site effect"),
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Power",
       x = "Estimated site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/cond_power_1sided.png"), width=200, height=125, units="mm")




