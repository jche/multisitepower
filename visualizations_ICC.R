
#####
# Visualize simulation results
#####

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
# Global settings
#####

ALPHA <- 0.1

hits2 <- hits %>%
  filter(n_bar %in% c(25, 50, 100),
         ICC %in% c(0, 0.3, 0.9),
         tau %in% c(0.01, 0.2, 0.8))


#####
# visualizing power vs. ATE, focused on MLMs
#####

# aggregated power results:
# per simulation setting & ATE value, how often do we reject the null using pvalue_one?
agg_hits =  hits %>%
  group_by( tau, n_bar, J, ICC, ATE ) %>%
  summarise( power = mean( pvalue_one <= ALPHA ) )

# plot power vs. ATE size
ggplot( agg_hits, aes( x=ATE, y=power, col=tau ) ) +
  facet_grid(ICC ~ n_bar, labeller=label_both ) +
  geom_line() +
  # geom_point(alpha=0.5) +
  geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = ALPHA, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true ATE"))
# ggsave("writeup/images/power_plot.png", width=200, height=150, units="mm")


#####
# visualizing power for sites with ATE=0.2, focused on MLMs
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- hits %>% 
  filter( ATE == 0.2 ) %>%
  group_by( tau, n_bar, J, ICC ) %>%
  summarise( n = n(),
             power = mean( pvalue_one <= ALPHA ) )

ggplot( circ20, aes( x=tau, y=power, col=ICC, group=interaction(n_bar, ICC) ) ) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept=0.8, lty="dashed") +
  facet_grid(~n_bar, labeller=label_both) +
  labs( title = glue("Power (\u03B1 = {ALPHA}) to detect a site with ATE=0.2"), 
        x = "tau" )
# ggsave("writeup/images/power_plot_ATE02.png", width=200, height=75, units="mm")


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
hits2 %>% 
  filter( ATE == 0.2 ) %>%
  ggplot( aes( x=ATE_hat ) ) +
  geom_density(aes(color=tau, group=tau), 
                 position="identity", fill="transparent") +
  geom_vline( xintercept = 0.2, col="red" ) +
  facet_grid( ICC ~ n_bar, labeller = label_both ) +
  labs(title = "Estimated ATEs for sites with ATE=0.2")
# ggsave("writeup/images/power_plot_ATE02_dens.png", width=200, height=150, units="mm")



#####
# visualizing power vs. ATE, comparing MLMs to single-site estimates
#####

agg_hits =  hits %>%
  group_by( tau, n_bar, J, ICC, ATE ) %>%
  summarise( power = mean( pvalue_one <= ALPHA ),
             power_single = mean(pvalue_single <= ALPHA) )

# power vs. ATE, comparing multi-site to single-site
power_single_key <- agg_hits %>%
  group_by(n_bar, J, ICC, ATE) %>%
  summarize(power_single = mean(power_single))

agg_hits %>%
  ggplot(aes(x=ATE)) +
  geom_line(data = power_single_key, aes(y=power_single)) +
  geom_line(aes(y=power, group=tau, color=tau)) +
  facet_grid(ICC ~ n_bar, labeller = label_both) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  geom_hline(aes(yintercept=0.8), lty="dashed") +
  geom_hline(aes(yintercept=ALPHA), lty="dashed") +
  geom_vline(aes(xintercept=0)) +
  labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true ATE"),
       subtitle = "Black line = single-site estimates")

# ggsave("writeup/images/power_plot_comp.png", width=200, height=150, units="mm")


# # bar chart of reject/no-reject for each site
# #  - note: ICC doesn't seem to have much of an effect here
# hits %>%
#   mutate(rej = pvalue_one <= 0.05,
#          rej_single = pvalue_single <= 0.05) %>%
#   group_by(ATE, tau, n_bar) %>%
#   summarize(both = sum(rej & rej_single),
#             neither = sum(!rej & !rej_single),
#             multi_only = sum(rej & !rej_single),
#             single_only = sum(!rej & rej_single)) %>%
#   pivot_longer(cols = c(both, neither, multi_only, single_only),
#                names_to = "rejections",
#                values_to = "n") %>%
#   # mutate(rejections = factor(rejections, levels = c("both", "neither", "multi_only", "single_only"))) %>%
#   ggplot(aes(x=ATE, y=n, fill=rejections)) +
#   geom_col(position = "identity", alpha=0.5) +
#   facet_grid(tau ~ n_bar, labeller=label_both)
# # ggsave("plots/reject_comp.png", width=9, height=4)


#####
# visualizing power for sites with ATE=0.2, comparing MLMs to single-site
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- hits %>% 
  filter( ATE == 0.2 ) %>%
  group_by( tau, n_bar, J, ICC ) %>%
  summarise( n = n(),
             power = mean( pvalue_one <= ALPHA ),
             power_single = mean(pvalue_single <= ALPHA))

circ20 %>%
  pivot_longer(cols = c(power, power_single),
               names_to = "estimator",
               values_to = "power") %>%
  mutate(estimator = ifelse(estimator == "power", "MLM", "single-site"))

power_single_key20 <- circ20 %>%
  group_by(n_bar, J, ICC) %>%
  summarize(power_single = mean(power_single))

circ20 %>%
ggplot( aes( col=ICC, group=interaction(n_bar, ICC) ) ) +
  geom_point(aes(x=tau, y=power)) +
  geom_line(aes(x=tau, y=power)) +
  geom_segment(data = power_single_key20, 
               aes(y=power_single, yend=power_single, x=1, xend=4),
               lty = "dashed") +
  geom_hline(yintercept=0.8, lty="dashed") +
  facet_grid(~n_bar, labeller=label_both) +
  labs( title = glue("Power (\u03B1 = {ALPHA}) to detect a site with ATE=0.2"),
        subtitle = "Solid lines = MLMs, dashed lines = single-site",
        x = "tau" )
# ggsave("writeup/images/power_plot_comp_ATE02.png", width=200, height=75, units="mm")


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
hits %>% 
  filter( ATE == 0.2, ICC==0 ) %>%
  pivot_longer(cols = c(ATE_hat, ATE_hat_single),
               names_to = "estimator",
               values_to = "ATE_hat") %>%
  mutate(estimator = ifelse(estimator == "ATE_hat", "MLM", "single-site")) %>%
  ggplot( aes( x=ATE_hat, group=estimator, color=estimator ) ) +
  facet_grid( tau ~ n_bar, labeller = label_both ) +
  geom_density(position = "identity", fill="transparent") +
  geom_vline( xintercept = 0.2, col="red" ) +
  labs(title = "Estimated ATEs for sites with ATE=0.2",
       subtitle = "(ICC = 0)")
# ggsave("writeup/images/power_plot_comp_ATE02_dens.png", width=200, height=150, units="mm")



#####
# RMSE plots
#####

# 
hits %>%
  group_by(n_bar, J, ICC, tau, tx_var, runID) %>%
  summarize(rmse = sqrt(mean((ATE - ATE_hat)^2)),
            rmse_single = sqrt(mean((ATE - ATE_hat_single)^2))) %>%
  ggplot(aes(x=rmse, y=rmse_single, col=ICC)) +
  geom_point() +
  facet_grid(tau ~ n_bar) +
  geom_abline()

# RMSE for tau_j=0.2 sites, across sims
hits %>%
  filter(ATE == 0.2) %>%
  group_by(n_bar, J, ICC, tau, tx_var) %>%
  summarize(rmse = sqrt(mean((ATE - ATE_hat)^2)),
            rmse_single = sqrt(mean((ATE - ATE_hat_single)^2))) %>%
  ggplot(aes(x=rmse, y=rmse_single, col=ICC)) +
  geom_point() +
  facet_grid(~ n_bar) +
  geom_abline()

