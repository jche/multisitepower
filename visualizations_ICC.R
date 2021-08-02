
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

FNAME <- "sim_results_ICC"

fname <- glue("results/", FNAME, ".csv")
hits <- read_csv(fname)


#####
# visualizing power vs. ATE, focused on MLMs
#####

# aggregated power results:
# per simulation setting & ATE value, how often do we reject the null using pvalue_one?
agg_hits =  hits %>%
  group_by( tau, n_bar, J, ICC, ATE ) %>%
  summarise( power = mean( pvalue_one <= 0.05 ) )

# plot power vs. ATE size
ggplot( agg_hits, aes( x=ATE, y=power, col = as.factor(tau) ) ) +
  facet_grid(ICC ~ n_bar, labeller=label_both ) +
  # geom_smooth(se=F) +
  geom_line() +
  geom_point(alpha=0.5) +
  geom_hline( yintercept = 0.8 ) +
  geom_vline( xintercept = 0 )
# ggsave("plots/power_plot.png")


#####
# visualizing power for sites with ATE=0.2, focused on MLMs
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- hits %>% 
  filter( ATE == 0.2 ) %>%
  group_by( tau, n_bar, J ) %>%
  summarise( n = n(),
             power = mean( pvalue_one <= 0.05 ) )

ggplot( circ20, aes( tau, power, col=as.factor( n_bar) ) ) +
  geom_point() +
  geom_line() +
  labs( title = "Power to detect a site with 0.20 ATE", 
        x = "Average ATE across sites" )
# ggsave("plots/power_plot_ATE02.png")


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
hits %>% 
  filter( ATE == 0.2 ) %>%
  ggplot( aes( x=ATE_hat ) ) +
  geom_histogram(aes(y = ..density.., fill=factor(ICC), group=ICC), 
                 position="identity", alpha=0.3) +
  geom_vline( xintercept = 0.2, col="red" ) +
  facet_grid( tau ~ n_bar, labeller = label_both )
# ggsave("plots/power_plot_ATE02_hist.png", width=9, height=4)


#####
# visualizing power vs. ATE, comparing MLMs to single-site estimates
#####

# aggregated power results:
# per simulation setting & ATE value, how often do we reject the null using pvalue_one?
agg_hits =  hits %>%
  group_by( tau, n_bar, J, ATE ) %>%
  summarise( power = mean( pvalue_one <= 0.05 ),
             power_single = mean(pvalue_single <= 0.05))

# power vs. ATE, comparing multi-site to single-site
agg_hits %>%
  pivot_longer(cols = c(power, power_single),
               names_to = "class",
               values_to = "power") %>%
  ggplot(aes(x=ATE, y=power, group=class, color=class)) +
  geom_point() +
  geom_line() +
  facet_grid(tau ~ n_bar, labeller = label_both) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  geom_hline(aes(yintercept=0.8)) +
  geom_vline(aes(xintercept=0))
# ggsave("plots/power_plot_comp.png")

# bar chart of reject/no-reject for each site
hits %>%
  mutate(rej = pvalue_one <= 0.05,
         rej_single = pvalue_single <= 0.05) %>%
  group_by(ATE, tau, n_bar) %>%
  summarize(both = sum(rej & rej_single),
            neither = sum(!rej & !rej_single),
            multi_only = sum(rej & !rej_single),
            single_only = sum(!rej & rej_single)) %>%
  pivot_longer(cols = c(both, neither, multi_only, single_only),
               names_to = "rejections",
               values_to = "n") %>%
  # mutate(rejections = factor(rejections, levels = c("both", "neither", "multi_only", "single_only"))) %>%
  ggplot(aes(x=ATE, y=n, fill=rejections)) +
  geom_col(position = "identity", alpha=0.5) +
  facet_grid(tau ~ n_bar, labeller=label_both)
# ggsave("plots/reject_comp.png", width=9, height=4)

# # same as above, with proportions
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
#   ggplot(aes(x=ATE, y=n, fill=rejections)) +
#   geom_col(position = "fill") +
#   facet_grid(tau ~ n_bar, labeller=label_both)


#####
# visualizing power for sites with ATE=0.2, comparing MLMs to single-site
#####

# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
hits %>% 
  filter( ATE == 0.2 ) %>%
  pivot_longer(cols = c(ATE_hat, ATE_hat_single),
               names_to = "class",
               values_to = "ATE_hat") %>%
  ggplot( aes( x=ATE_hat, group=class, color=class, fill=class ) ) +
  facet_grid( tau ~ n_bar, labeller = label_both ) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha=0.5) +
  geom_vline( xintercept = 0.2, col="red" )
# ggsave("plots/hist_comp.png", width=9, height=4)


