
#####
# Power analysis for individual site effects in multisite trial
# 
# uses blkvar package: devtools::install_github("https://github.com/lmiratrix/blkvar")
#####

require(tidyverse)
require(glue)
require(tictoc)
# devtools::install_github("lmiratrix/blkvar")
library( blkvar )
library( arm )

#####
# simulation functions
#####

# function: (# obs, effect size) => (reject null? T/F)
#' 
#'
#' @param n number of observations per site (constant, for now)
#' @param J number of sites
#' @param tau overall ATE
#' @param ICC 
#'
#' @return tibble with J rows, 
#' @export
#'
#' @examples
one_sim <- function(n, J, tau, round_sites = 0.05) {
  
  browser()
  
  sdat = blkvar::generate_multilevel_data( n.bar=n, J=J,
                                           variable.n = FALSE,
                                           tau.11.star = 0.3,      # cross-site tx var
                                           gamma.10 = tau,         # cross-site ATE
                                           # ICC = ICC,              # ICC
                                           return.sites = TRUE,
                                           verbose = FALSE) %>%
    mutate(sid = as.character(1:n()),
           beta.1 = round(beta.1/round_sites) * round_sites)       # round to nearest round_sites
  
  # Note: generate_individual_data() is not in CRAN version of package
  dat = blkvar::generate_individual_data( sdat )
  head( dat )
  
  # Switch to FIRC model someday.
  # Or our bayesian modeling to get posteriors.
  # Currently: RIRC model (so no site-by-treatment interactions)
  mod = lmer( Yobs ~ 1 + Z + (1+Z|sid), data=dat )
  
  res = as.data.frame( coef( mod )$sid )     # random effects
  ses = as.data.frame( se.ranef(mod)$sid )   # standard errors
  
  res = tibble( sid = rownames( ranef( mod )$sid ),
                est = res$Z,
                SE = ses$Z,
                t = est / SE,
                pvalue = 2*pnorm( -abs(t) ) )
  
  res = left_join( res, sdat, by="sid" ) %>%
    dplyr::select( -u0, -u1, -beta.0 )  %>%
    rename( ATE_hat = est,
            ATE = beta.1 )
  
  res    
}


if ( F ) {
  foo <- function(n, J, tau, ICC) {
    
    sdat = blkvar::generate_multilevel_data( n.bar=n, J=J,
                                             variable.n = FALSE,
                                             tau.11.star = 0.3,      # cross-site tx var
                                             gamma.10 = tau,         # cross-site ATE
                                             ICC = ICC,      
                                             # zero.corr=T prevents breaking. If ICC is too low, there's not a lot of variation on the intercept, so if we have treatment variation but not site variation, we can't undo it to get tx/co variation to be the same.
                                             zero.corr=T,# ICC
                                             return.sites = TRUE,
                                             verbose = TRUE)
  }
  browser()
  os <- foo( n=20, J=1000, tau=0.2, ICC=0.1)
  os <- foo( n=20, J=1000, tau=0.2, ICC=0.4)
  
  # issue: currently, low ICC values throw an error
  # Q: how does covariate W affect things? Should I use gen_multi_data_no_cov()?
  
  
  os <- one_sim( n=20, J=1000, tau=0.2)
  mean( os$ATE )
  mean( os$ATE_hat )
  ggplot( os, aes( ATE, ATE_hat ) ) +
    geom_point() +
    geom_abline( slope=1, intercept= 0 )
}

# function: (# obs, effect size) => (power)
#  - runs one_sim NUMSIM times, so we can aggregate across runIDs to get the power
power_sim <- function(n, J, tau, NUMSIM = 250) {
  cat(glue("Working on n = {n}, tau = {tau}\n\n"))
  rs = tibble( runID = 1:NUMSIM )
  rs$data = map( rs$runID, ~one_sim(n, J, tau))
  rs = unnest( rs, cols = data )
  rs
}

if ( FALSE ) {
  power_sim( 20, 20, 0.2, 3 )
}


#####
# run power simulation
#####

# run power simulation:
#  - n = number of observations
#  - tau = true effect size
#  - (real power sim would have more knobs: J, site.size, ICC/variances, etc.)
df_sim <- expand_grid(
  n = c(25, 50, 75, 100),
  J = 20,
  # ICC knob
  tau = c(0.01, 0.2, 0.5, 0.8)
)

# run simulation: store power_sim() results in df_sim as list column
df_sim <- df_sim %>%
  rowwise() %>%
  mutate(data = list(power_sim(n, J, tau, NUMSIM = 100)))
# df_sim$data = pmap( df_sim, power_sim, NUMSIM = 100 )

# raw results: unnest df_sim & record pvalue_one per site
hits = df_sim %>% 
  rename( n_bar = n ) %>%
  unnest( cols=data ) %>%
  mutate( pvalue_one = pnorm( -t ))


FNAME <- "temp_results"
fname <- glue("results/", FNAME, ".csv")

if (FALSE) {
  if (file.exists(fname)) {
    max_runID <- read_csv(fname) %>%
      pull(runID) %>%
      max()
    
    hits %>% 
      mutate(runID = runID + max_runID) %>%
      write_csv(fname, append=T)
  } else {
    write_csv(hits, fname)
  }
}

if (TRUE) {
  hits <- read_csv(fname)
}


#####
# visualizing power vs. ATE
#####

# aggregated power results:
# per simulation setting & ATE value, how often do we reject the null using pvalue_one?
agg_hits =  hits %>%
  group_by( tau, n_bar, J, ATE ) %>%
  summarise( power = mean( pvalue_one <= 0.05 ) )

# plot power vs. ATE size
ggplot( agg_hits, aes( x=ATE, y=power, col = as.factor(tau) ) ) +
  facet_wrap( ~ n_bar, labeller=label_both ) +
  # geom_smooth(se=F) +
  geom_line() +
  geom_point(alpha=0.5) +
  geom_hline( yintercept = 0.8 ) +
  geom_vline( xintercept = 0 )
# ggsave("plots/power_plot.png")


ggplot( agg_hits, aes( x=ATE, y=power, col = as.factor(n_bar) ) ) +
  facet_wrap( ~ tau, labeller=label_both ) +
  # geom_smooth(se=F) +
  geom_line() +
  geom_point(alpha=0.5) +
  geom_hline( yintercept = 0.8 ) +
  geom_vline(xintercept = 0)
# ggsave("plots/power_plot_v2.png")


#####
# visualizing power for sites with ATE=0.2
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- hits %>% 
  filter( ATE == 0.2 ) %>%
  group_by( tau, n_bar, J ) %>%
  summarise( n = n(),
             power = mean( pvalue_one <= 0.05 ) )
head( circ20 )

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
  facet_grid( tau ~ n_bar, labeller = label_both ) +
  geom_histogram(aes(y = ..density..)) +
  geom_vline( xintercept = 0.2, col="red" )
# ggsave("plots/power_plot_ATE02_hist.png", width=9, height=4)




