
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

FNAME <- "sim_results_full"

fname <- glue("results/", FNAME, ".csv")
hits <- read_csv(fname) %>%
  mutate(n_bar = as.factor(n_bar),
         ICC = as.factor(ICC),
         tau = as.factor(tau),
         tx_var = as.factor(tx_var))


#####
# Global settings
#####

ALPHA <- 0.1

hits_list <- list(
  hits %>% filter(tx_var == 0), 
  hits %>% filter(tx_var == 0.3), 
  hits %>% filter(tx_var == 0.6))


#####
# visualizing power vs. ATE, focused on MLMs
#####

# aggregated power results:
# per simulation setting & ATE value, how often do we reject the null using pvalue_one?
agg_hits <- map(hits_list, function(df) {
  TX_VAR <- df$tx_var[1]
  
  p <- df %>%
    group_by( tau, n_bar, J, ICC, ATE ) %>%
    summarise( power = mean( pvalue_one <= ALPHA )) %>%
    ggplot( aes( x=ATE, y=power, col=tau ) ) +
      facet_grid(ICC ~ n_bar, labeller=label_both ) +
      geom_line() +
      geom_hline( yintercept = 0.8, lty = "dashed" ) +
      geom_hline( yintercept = ALPHA, lty = "dashed" ) +
      geom_vline( xintercept = 0 ) +
      coord_cartesian(xlim = c(-0.5, 1), ylim = c(0, 1)) +
      labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true ATE"),
           subtitle = glue("Treatment variation: {TX_VAR}"))
  
  if (TX_VAR == 0) {
    p + geom_point(alpha=0.5)
  } else {
    p
  }
} )

# plot power vs. ATE size
agg_hits[[1]]
agg_hits[[2]]
agg_hits[[3]]

ggsave("writeup/images/power_plot6.png", width=200, height=150, units="mm")


#####
# visualizing power for sites with ATE=0.2, focused on MLMs
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- map(hits_list, function(df) {
  TX_VAR <- df$tx_var[1]
  
  df %>% 
    filter( ATE == 0.2 ) %>%
    group_by( tau, n_bar, J, ICC, tx_var ) %>%
    summarise( n = n(),
               power = mean( pvalue_one <= ALPHA ) ) %>%
    ggplot( aes( x=tau, y=power, col=ICC, group=interaction(n_bar, ICC) ) ) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept=0.8, lty="dashed") +
    facet_grid(~n_bar, labeller=label_both) +
    labs( title = glue("Power (\u03B1 = {ALPHA}) to detect a site with ATE=0.2"), 
          subtitle = glue("Treatment variation: {TX_VAR}"),
          x = "tau" )
})

circ20[[1]]
circ20[[2]]
circ20[[3]]

ggsave("writeup/images/power_plot_ATE02.png", width=200, height=75, units="mm")


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
ate_hists <- map(hits_list, function(df) {
  TX_VAR <- df$tx_var[1]
  
  df %>% 
    filter( ATE == 0.2 ) %>%
    ggplot( aes( x=ATE_hat ) ) +
    geom_density(aes(color=tau, group=tau), 
                 position="identity", fill="transparent") +
    geom_vline( xintercept = 0.2, col="red" ) +
    facet_grid( ICC ~ n_bar, labeller = label_both ) +
    labs(title = "Estimated ATEs for sites with ATE=0.2",
         subtitle = glue("Treatment variation: {TX_VAR}"))
}) 

ate_hists[[1]]
ate_hists[[2]]
ate_hists[[3]]

# ggsave("writeup/images/power_plot_ATE02_dens.png", width=200, height=150, units="mm")



#####
# visualizing power vs. ATE, comparing MLMs to single-site estimates
#####

agg_hits_comp <- map(hits_list, function(df) {
  TX_VAR <- df$tx_var[1]
  
  p <- df %>%
    group_by( tau, n_bar, J, ICC, ATE ) %>%
    summarise( power = mean( pvalue_one <= ALPHA ),
               power_single = mean(pvalue_single <= ALPHA)) %>%
    pivot_longer(cols = c(power, power_single),
                 names_to = "estimator",
                 values_to = "power") %>%
    mutate(estimator = ifelse(estimator == "power", "MLM", "single-site")) %>%
    ggplot(aes(x=ATE, y=power, group=interaction(estimator, tau), color=tau)) +
    # geom_point() +
    geom_line(aes(lty=estimator)) +
    facet_grid(ICC ~ n_bar, labeller = label_both) +
    coord_cartesian(xlim = c(-0.5, 1)) +
    geom_hline(aes(yintercept=0.8), lty="dashed") +
    geom_hline(aes(yintercept=ALPHA), lty="dashed") +
    geom_vline(aes(xintercept=0)) +
    labs(title = glue("Power (\u03B1 = {ALPHA}) vs. true ATE"),
         subtitle = glue("Treatment variation: {TX_VAR}"))
  
  if (TX_VAR == 0) {
    p + geom_point(alpha=0.5)
  } else {
    p
  }
} )

# plot power vs. ATE size
agg_hits_comp[[1]]
agg_hits_comp[[2]]
agg_hits_comp[[3]]


#####
# visualizing power for sites with ATE=0.2, comparing MLMs to single-site
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20_comp <- map(hits_list, function(df) {
  TX_VAR <- df$tx_var[1]
  
  df %>% 
    filter( ATE == 0.2 ) %>%
    group_by( tau, n_bar, J, ICC, tx_var ) %>%
    summarise( n = n(),
               power = mean( pvalue_one <= ALPHA ),
               power_single = mean(pvalue_single <= ALPHA)) %>%
    pivot_longer(cols = c(power, power_single),
                 names_to = "estimator",
                 values_to = "power") %>%
    mutate(estimator = ifelse(estimator == "power", "MLM", "single-site")) %>%
    ggplot( aes( x=tau, y=power, col=ICC, group=interaction(n_bar, ICC, estimator) ) ) +
    geom_point() +
    geom_line(aes(lty = estimator)) +
    geom_hline(yintercept=0.8, lty="dashed") +
    facet_grid(~n_bar, labeller=label_both) +
    labs( title = glue("Power (\u03B1 = {ALPHA}) to detect a site with ATE=0.2"), 
          subtitle = glue("Treatment variation: {TX_VAR}"),
          x = "tau" )
})

circ20_comp[[1]]
circ20_comp[[2]]
circ20_comp[[3]]

ggsave("writeup/images/power_plot_ATE02.png", width=200, height=75, units="mm")


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
ate_hists_comp <- map(hits_list, function(df) {
  TX_VAR <- df$tx_var[1]
  ICC_VAL <- 0
  
  df %>% 
    filter( ATE == 0.2, ICC==ICC_VAL ) %>%
    pivot_longer(cols = c(ATE_hat, ATE_hat_single),
                 names_to = "estimator",
                 values_to = "ATE_hat") %>%
    mutate(estimator = ifelse(estimator == "ATE_hat", "MLM", "single-site")) %>%
    ggplot( aes( x=ATE_hat, group=estimator, color=estimator ) ) +
    facet_grid( tau ~ n_bar, labeller = label_both ) +
    geom_density(position = "identity", fill="transparent") +
    geom_vline( xintercept = 0.2, col="red" ) +
    labs(title = "Estimated ATEs for sites with ATE=0.2",
         subtitle = glue("ICC = {ICC_VAL}, Treatment variation: {TX_VAR}"))
}) 

ate_hists_comp[[1]]
ate_hists_comp[[2]]
ate_hists_comp[[3]]

# ggsave("writeup/images/power_plot_ATE02_dens.png", width=200, height=150, units="mm")

