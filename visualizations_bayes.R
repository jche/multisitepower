
#####
# Visualize simulation results
#####

library( arm )   # note: loads MASS package
library( blkvar )

require(tidyverse)
require(glue)
require(tictoc)
require(latex2exp)


#####
# Load simulation results
#####

source("preprocess_cluster_results.R")
ALPHA <- 0.1

hits <- res %>%
  mutate(nbar = as.factor(nbar),
         J = as.factor(J),
         ICC = as.factor(ICC),
         tau = as.factor(tau))

if ("df_single" %in% colnames(hits)) {
  df_key <- hits %>%
    select(nbar, J, ICC, tau, tx_var, runID, sid, df_single)
}

# Get ATEs per method (long data)
ATEhats <- hits %>%
  group_by(nbar, J, ICC, tau, tx_var, runID) %>%
  pivot_longer(contains("ATEhat"), names_to = "method", values_to = "ATEhat") %>%
  select(ATE, sid, method, ATEhat) %>%
  mutate(method = str_sub(method, 8))

# Get SEs per method (long data)
SEs <- hits %>%
  mutate(SE_rirc = SE_rirc_rand,
         SE_firc = SE_firc_rand) %>%
  # mutate(SE_rirc = sqrt(SE_rirc_fixed^2 + SE_rirc_rand^2),
  #        SE_firc = sqrt(SE_firc_fixed^2 + SE_firc_rand^2)) %>%   # using both SEs
  group_by(nbar, J, ICC, tau, tx_var, runID) %>%
  pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
  select(sid, method, SE) %>%
  mutate(method = str_sub(method, 4)) %>%
  filter(method %in% c("single", "firc", "rirc", "bayesnorm"))

# join data
if (exists("df_key")) {
  tidy_results <- ATEhats %>%
    cbind(SE = SEs$SE) %>%
    # left_join(SEs, by=c("nbar", "J", "ICC", "tau", "tx_var", "runID", "sid", "method")) %>%
    left_join(df_key, by=c("nbar", "J", "ICC", "tau", "tx_var", "runID", "sid")) %>%
    mutate(pvalue_one = ifelse(method == "single", pt(-ATEhat/SE, df=df_single), pnorm(-ATEhat/SE)),
           interval = ifelse(method == "single", qt(1-ALPHA/2, df=df_single)*SE, qnorm(1-ALPHA/2)*SE),
           covered = (ATE <= ATEhat+interval) & (ATE >= ATEhat-interval),
           method = ifelse(method == "single", "Single",
                           ifelse(method == "firc", "FIRC",
                                  ifelse(method == "rirc", "RIRC", "Bayes")))) %>%
    select(-df_single) %>%
    ungroup()
} else {
  tidy_results <- ATEhats %>%
    cbind(SE = SEs$SE) %>%
    mutate(pvalue_one = pnorm(-ATEhat/SE),
           interval = qnorm(1-ALPHA/2)*SE,
           covered = (ATE <= ATEhat+interval) & (ATE >= ATEhat-interval),
           method = ifelse(method == "single", "Single",
                           ifelse(method == "firc", "FIRC",
                                  ifelse(method == "rirc", "RIRC", "Bayes")))) %>%
    ungroup()
}


rm(res)
# rm(hits)
rm(ATEhats)
rm(SEs)


#####
# checking coverage
#####

# plot coverage for a single trial
tidy_results %>% 
  filter(nbar == 25, J == 25, tau == 0, runID == 7) %>%
  # filter(nbar == 25, J == 20, ICC == 0, tau == 0.01, runID == 2) %>%
  ggplot(aes(x = ATE, color = method)) +
  geom_point(aes(y = ATEhat)) +
  geom_errorbar(aes(ymin = ATEhat - qnorm(1-ALPHA/2)*SE, 
                    ymax = ATEhat + qnorm(1-ALPHA/2)*SE)) +
  facet_wrap(~method) +
  geom_abline(slope = 1)

# plot coverage vs. ATE size
tidy_results %>%
  # filter(tau == 0) %>%
  # filter(J != 20, tau == 0) %>%
  filter(method != "Single", tau == 0) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, ATE) %>%
  summarize(coverage = mean(covered)) %>%
  filter(ATE >= as.numeric(as.character(tau)) - 1,
         ATE <= as.numeric(as.character(tau)) + 1) %>%
  ggplot() +
  geom_line(aes(x=ATE, y=coverage, color=method, group=method)) +
  facet_grid(J ~ nbar, labeller=label_both) +
  geom_hline(yintercept=1-ALPHA, lty="dashed") +
  labs(title = "Conditional coverage rates",
       y = "Coverage",
       x = "True site ATE",
       color = "Method")
# ggsave(glue("writeup/images/coverage_plot.png"), width=200, height=125, units="mm")

# plot EB coverage
tidy_results %>%
  # filter(tau == 0) %>%
  filter(method != "Single") %>%
  group_by(nbar, J, ICC, tau, tx_var, method) %>%
  summarize(coverage = mean(covered)) %>%
  ggplot(aes(x=nbar, y=coverage, color=method)) +
  geom_point() +
  geom_line(aes(group=method)) +
  facet_grid( ~ J, labeller=label_both, scales="free") +
  geom_hline(yintercept=1-ALPHA, lty="dashed")
# ggsave(glue("writeup/images/coverage_plot_eb.png"), width=200, height=75, units="mm")


#####
# checking power
#####

J_plot <- 300

# Plot power vs. ATE size, for a specified true tau and ICC
power_df <- tidy_results %>%
  filter(J == J_plot) %>%
  mutate(reject = pvalue_one < ALPHA) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, ATE) %>%
  summarize(power = mean(reject))
power_df %>%
  ggplot(aes(x = ATE, y = power, col = method, group=method ) ) +
  facet_grid(tau ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = ALPHA, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 0.6)) +
  labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true site effect"),
       subtitle = glue("J = {J_plot}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/power_plot_J300.png"), width=200, height=125, units="mm")

# Plot power gaps
power_df_single <- power_df %>%
  filter(method == "Single") %>%
  rename(power_single = power) %>%
  ungroup() %>%
  select(-method)

power_df %>%
  left_join(power_df_single, by=c("nbar", "J", "ICC", "tau", "tx_var", "ATE")) %>%
  mutate(power_gap = power - power_single) %>%
  ggplot(aes(x = ATE, y = power_gap, col = method, group=method ) ) +
  facet_grid(tau ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true site effect"),
       subtitle = glue("J = {J_plot}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/power_plot_J300_diff.png"), width=200, height=125, units="mm")


### FOR POWER ANALYSIS EXAMPLE

ALPHA <- 0.1
tidy_results %>%
  mutate(nbar = as.numeric(as.character(nbar))) %>%
  filter(round(ATE, 2) %in% c(0.2, 0.3, 0.4)) %>%
  mutate(reject = pvalue_one < ALPHA,
         ATE = as.factor(ATE)) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, ATE) %>%
  summarize(power = mean(reject)) %>%
ggplot(aes(x = nbar, y = power, col = method, group=method ) ) +
  facet_grid(~ ATE, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0.8, lty = "dashed" ) +
  labs(title = glue("Power (\u03B1 = {ALPHA}) vs. site size"),
       subtitle = "Overall average effect: 0.2",
       y = "Power",
       x = "Site size",
       color = "Method")
ggsave(glue("writeup/images/power_plot_ex.png"), width=200, height=75, units="mm")


#####
# visualizing power for sites with ATE=0.2
#####

ATE_of_interest <- 0.3

# for sites with ATE=0.2, plot power (across sim settings)
circ <- tidy_results %>% 
  filter( abs(ATE-ATE_of_interest) < 0.01 ) %>%
  group_by( tau, nbar, J, ICC, tx_var, method ) %>%
  summarise( n = n(),
             power = mean( pvalue_one <= ALPHA ) )

circ %>%
  filter(J == 300) %>%
  ggplot( aes( x=nbar, y=power, col=method, group=method ) ) +
  geom_point() +
  geom_line(aes(group=method)) +
  geom_hline(yintercept=1-ALPHA, lty="dashed") +
  facet_grid(~ tau, labeller=label_both) +
  labs( title = glue("Power (\u03B1 = {ALPHA}) to detect a site with ATE={ATE_of_interest}"), 
        x = "nbar" )
ggsave("writeup/images/power_plot_ATE03.png", width=200, height=75, units="mm")

# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
tidy_results %>%
  filter(abs(ATE-ATE_of_interest) < 0.01) %>%
  filter(J == 300) %>%
  ggplot( aes( x=ATEhat ) ) +
  geom_density(aes(color=method, group=method), 
               position="identity", fill="transparent") +
  geom_vline( xintercept = 0.2, col="red" ) +
  coord_cartesian(xlim = c(-1, 1)) +
  facet_grid( tau ~ nbar, labeller = label_both ) +
  labs(title = glue("Estimated ATEs for sites with ATE={ATE_of_interest}"))
ggsave("writeup/images/power_plot_ATE03_dens.png", width=200, height=100, units="mm")


#####
# RMSE plots
#####

# boxplots of RMSE values
tidy_results %>%
  mutate(method = factor(method, levels=c("Single", "Bayes", "FIRC", "RIRC"))) %>%
  # filter(tau == 0) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, runID) %>%
  summarize(rmse = sqrt(mean((ATE - ATEhat)^2))) %>%
  ggplot(aes(x=method, y=rmse, color=method)) +
  geom_boxplot() +
  facet_grid(J ~ nbar, labeller = label_both) +
  labs(y = "RMSE",
       x = "Method") +
  guides(color = "none")
ggsave("writeup/images/rmse_plot.png", width=200, height=150, units="mm")



#####
# check out standard errors
#####

tidy_results %>%
  # filter(tau == 0) %>%
  group_by(nbar, J, ICC, tau, tx_var, method) %>%
  summarize(SE = mean(SE)) %>%
  ggplot(aes(x = nbar, y = SE, color = method)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_grid( ~ J, labeller=label_both)


#####
# troubleshooting
#####

# compute empirical SEs of effect estimates
#  - weird...? for RIRC/Single we're correct, but for Bayes we overestimate and for FIRC we underestimate...
tidy_results %>%
  ungroup() %>%
  group_by(nbar, J, ICC, tau, tx_var, ATE, method) %>%
  summarize(n = n(),
            emp_SE = sd(ATEhat),
            avg_SE = mean(SE)) %>%
  filter(n > 25) %>%
  ggplot() +
  geom_point(aes(x=emp_SE, y=avg_SE, color=nbar)) +
  geom_abline() +
  facet_grid(J~method)


# checking singularities: RIRC singular more often than FIRC
#  - When J=10 or nbar=10, we get singularities
# checking ESS low: always running into issues when nbar=300 and J>=25!
#  - this is because of reparameterization... recall that our parameterization doesn't work well when the data are super informative!
hits %>%
  group_by(nbar, J, tau) %>%
  summarize(prop_singular_firc = mean(is_singular_firc),
            prop_singular_rirc = mean(is_singular_rirc),
            prop_ESS_low = mean(ESS_low)) %>%
  print(n=100)


#####
# overall power
#####

res_overall %>%
  filter(method != "Single") %>%
  mutate(pvalue_one = pnorm(-ATEhat/SE)) %>%
  group_by(nbar, J, ICC, tau, tx_var, method) %>%
  summarize(power = mean(pvalue_one < ALPHA)) %>%
ggplot(aes(x=nbar, y=power, color=method)) +
  geom_point(aes(group=method)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.8), lty = "dashed") +
  labs(y = "Power",
       x = "Average site size",
       color = "Method")
ggsave(glue("writeup/images/power_plot_overall_tau02.png"), width=200, height=75, units="mm")

res_overall %>%
  filter(method != "Single") %>%
  mutate(interval = qnorm(1-ALPHA/2)*SE,
         covered = (tau <= ATEhat+interval) & (tau >= ATEhat-interval)) %>%
  group_by(nbar, J, ICC, tau, tx_var, method) %>%
  summarize(coverage = mean(covered)) %>%
  ggplot(aes(x=nbar, y=coverage, color=method)) +
  geom_point(aes(group=method)) +
  geom_line() +
  geom_hline(aes(yintercept = 1-ALPHA), lty = "dashed") +
  labs(y = "Coverage",
       x = "Average site size",
       color = "Method")
ggsave(glue("writeup/images/coverage_plot_overall_tau02.png"), width=200, height=75, units="mm")

