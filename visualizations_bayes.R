
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

# FNAME <- "sim_results_bayes"
FNAME <- "sim_results_bayes3"

fname <- glue("results/", FNAME, ".csv")
hits <- read_csv(fname) %>%
  mutate(nbar = as.factor(nbar),
         J = as.factor(J),
         ICC = as.factor(ICC),
         tau = as.factor(tau))

# clean data
ATEhats <- hits %>%
  group_by(nbar, J, ICC, tau, tx_var, runID) %>%
  pivot_longer(contains("ATEhat"), names_to = "method", values_to = "ATEhat") %>%
  select(ATE, sid, method, ATEhat) %>%
  mutate(method = str_sub(method, 8))

SEs <- hits %>%
  # mutate(SE_rirc = SE_rirc_rand,
         # SE_firc = SE_firc_rand) %>%
  # mutate(SE_rirc = sqrt(SE_rirc_fixed^2 + SE_rirc_rand^2),
  #        SE_firc = sqrt(SE_firc_fixed^2 + SE_firc_rand^2)) %>%   # using both SEs
  group_by(nbar, J, ICC, tau, tx_var, runID) %>%
  pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
  select(sid, method, SE) %>%
  mutate(method = str_sub(method, 4))
tidy_results <- ATEhats %>%
  left_join(SEs, by=c("nbar", "J", "ICC", "tau", "tx_var", "runID", "sid", "method")) %>%
  mutate(pvalue_one = pnorm(-ATEhat/SE),
         method = ifelse(method == "single", "Single",
                         ifelse(method == "firc", "FIRC",
                                ifelse(method == "rirc", "RIRC", "Bayes"))))

SE1s <- hits %>%
  mutate(SE_rirc = SE_rirc_rand,
         SE_firc = SE_firc_rand) %>%
  # mutate(SE_rirc = sqrt(SE_rirc_fixed^2 + SE_rirc_rand^2),
  #        SE_firc = sqrt(SE_firc_fixed^2 + SE_firc_rand^2)) %>%   # using both SEs
  group_by(nbar, J, ICC, tau, tx_var, runID) %>%
  pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
  select(sid, method, SE) %>%
  mutate(method = str_sub(method, 4))
tidy_results1 <- ATEhats %>%
  left_join(SE1s, by=c("nbar", "J", "ICC", "tau", "tx_var", "runID", "sid", "method")) %>%
  mutate(pvalue_one = pnorm(-ATEhat/SE),
         method = ifelse(method == "single", "Single",
                         ifelse(method == "firc", "FIRC",
                                ifelse(method == "rirc", "RIRC", "Bayes"))))


#####
# Global settings
#####

ALPHA <- 0.1


#####
# checking coverage
#####

# plot coverage for a single trial
tidy_results %>% 
  filter(nbar == 25, J == 25, ICC == 0.6, tau == 0.01, runID == 7) %>%
  # filter(nbar == 25, J == 20, ICC == 0, tau == 0.01, runID == 2) %>%
  ggplot(aes(x = ATE, color = method)) +
  geom_point(aes(y = ATEhat)) +
  geom_errorbar(aes(ymin = ATEhat - qnorm(1-ALPHA/2)*SE, 
                    ymax = ATEhat + qnorm(1-ALPHA/2)*SE)) +
  facet_wrap(~method) +
  geom_abline(slope = 1)

# plot coverage vs. ATE size
tidy_results1 %>%
  filter(J == 50) %>%
  mutate(interval = qnorm(1-ALPHA/2)*SE,
         covered = (ATE <= ATEhat+interval) & (ATE >= ATEhat-interval)) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, ATE) %>%
  summarize(coverage = mean(covered)) %>%
  filter(ICC == 0,
         ATE >= as.numeric(as.character(tau)) - 1,
         ATE <= as.numeric(as.character(tau)) + 1) %>%
  ggplot() +
  geom_line(aes(x=ATE, y=coverage, color=method, group=method)) +
  facet_grid(tau ~ nbar, labeller=label_both) +
  geom_hline(yintercept=1-ALPHA, lty="dashed") +
  labs(title = "Conditional coverage rates",
       y = "Coverage",
       x = "True site ATE",
       color = "Method")
ggsave(glue("writeup/images/coverage_plot.png"), width=200, height=125, units="mm")

# plot EB coverage
tidy_results %>%
  mutate(interval = qnorm(1-ALPHA/2)*SE,
         covered = (ATE <= ATEhat+interval) & (ATE >= ATEhat-interval)) %>%
  group_by(nbar, J, ICC, tau, tx_var, method) %>%
  summarize(coverage = mean(covered)) %>%
  ggplot(aes(x=tau, y=coverage, color=method)) +
  geom_point(aes(pch=as.factor(J))) +
  geom_line(aes(group=interaction(method, J))) +
  facet_grid(ICC ~ nbar, labeller=label_both, scales="free") +
  geom_hline(yintercept=1-ALPHA, lty="dashed")


#####
# checking power
#####

# Plot power vs. ATE size, for a specified true tau and ICC
tau_plot <- 0.2

tidy_results %>%
  filter((ICC == 0 & nbar == 25) | 
           (ICC == 0.3 & nbar == 50) | 
           (ICC == 0.6 & nbar == 75)) %>%
  filter(tau == tau_plot) %>%
  mutate(reject = pvalue_one < ALPHA) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, ATE) %>%
  summarize(power = mean(reject)) %>%
  ggplot(aes(x = ATE, y = power, col = method, group=method ) ) +
  facet_grid(J ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = ALPHA, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true site effect"),
       subtitle = glue("Overall ATE = {tau_plot}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method")
ggsave(glue("writeup/images/power_plot_{tau_plot*10}.png"), width=200, height=150, units="mm")


#####
# visualizing power for sites with ATE=0.2
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- tidy_results %>% 
  filter( ATE == 0.2 ) %>%
  group_by( tau, nbar, J, ICC, tx_var, method ) %>%
  summarise( n = n(),
             power = mean( pvalue_one <= ALPHA ) )

ggplot( circ20, aes( x=tau, y=power, col=method, group=method ) ) +
  geom_point(aes(pch = J)) +
  geom_line(aes(group=interaction(J, method))) +
  geom_hline(yintercept=1-ALPHA, lty="dashed") +
  facet_grid(ICC~nbar, labeller=label_both) +
  labs( title = glue("Power (\u03B1 = {ALPHA}) to detect a site with ATE=0.2"), 
        x = "tau" )
# ggsave("writeup/images/power_plot_ATE02.png", width=200, height=75, units="mm")

# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
tau_plot <- 0.2
tidy_results %>%
  filter(ATE == 0.2) %>%
  filter(ICC == 0, tau == 0.5) %>%
  ggplot( aes( x=ATEhat ) ) +
  geom_density(aes(color=method, group=method), 
               position="identity", fill="transparent") +
  geom_vline( xintercept = 0.2, col="red" ) +
  facet_grid( J ~ nbar, labeller = label_both ) +
  labs(title = "Estimated ATEs for sites with ATE=0.2")

tidy_results %>%
  filter(ATE == 0.2) %>%
  filter(J == 50) %>%
  filter((ICC == 0 & nbar == 25) | 
           (ICC == 0.3 & nbar == 50) | 
           (ICC == 0.6 & nbar == 75)) %>%
  ggplot( aes( x=ATEhat ) ) +
  geom_density(aes(color=method, group=method), 
               position="identity", fill="transparent") +
  geom_vline( xintercept = 0.2, col="red" ) +
  facet_grid( tau ~ nbar, labeller = label_both ) +
  labs(title = "Estimated ATEs for sites with ATE=0.2",
       color = "Method",
       x = "Estimated site ATE",
       y = "Density")
ggsave("writeup/images/power_plot_ATE02_dens.png", width=200, height=150, units="mm")


#####
# RMSE plots
#####

# boxplots of RMSE values
tidy_results %>%
  filter((ICC == 0 & nbar == 25) | 
           (ICC == 0.3 & nbar == 50) | 
           (ICC == 0.6 & nbar == 75)) %>%
  group_by(nbar, J, ICC, tau, tx_var, method, runID) %>%
  summarize(rmse = sqrt(mean((ATE - ATEhat)^2))) %>%
  ggplot(aes(x=method, y=rmse, color=method)) +
  geom_boxplot() +
  facet_grid(tau ~ nbar, labeller = label_both) +
  labs(y = "RMSE",
       x = "Method") +
  guides(color = F)
ggsave("writeup/images/rmse_plot.png", width=200, height=150, units="mm")



#####
# check out standard errors
#####

SEs %>%
  filter(method %in% c("bayesnorm", "firc", "rirc", "single", "firc_rand", "rirc_rand")) %>%
  filter(ICC == 0, tau == 0.01) %>%
  group_by(nbar, J, ICC, tau, tx_var, method) %>%
  summarize(SE = mean(SE)) %>%
  ggplot(aes(x = nbar, y = SE, color = method)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~J)
