
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
# devtools::install_github("lmiratrix/blkvar")
library( blkvar )
library( arm )

#####
# simulation functions
#####

# function: (# obs, effect size) => (reject null? T/F)
one_sim <- function(n, J, tau) {
    
    sdat = blkvar::generate_multilevel_data( n.bar=n, J=J,
                                             variable.n = FALSE,
                                             tau.11.star = 0.3,      # tx var
                                             gamma.10 = tau,         # ATE
                                             return.sites = TRUE ) 
    head( sdat )
    sdat$sid = as.character( 1:nrow(sdat) )
    
    # Note: gen_indiv_data() is not in CRAN version of package
    dat = blkvar::generate_individual_data( sdat )
    head( dat )
    
    # Switch to FIRC model someday.
    # Or our bayesian modeling to get posteriors.
    # Currently: RIRC model (so no site-by-treatment interactions)
    mod = lmer( Yobs ~ 1 + Z + (1+Z|sid), data=dat )
    
    # Q: why are standard errors of the random effects constant?
    # A: site sizes are constant
    res = as.data.frame( coef( mod )$sid )     # random effects
    ses = as.data.frame( se.ranef(mod)$sid )   # standard errors
    
    # # Q: does SE match the formula?
    # #  - Q: do we just plug in the residual variance for sigma_m?
    # # A: no! formula is for something global, SEs are for site random effects...?
    # # A: no! formula doesn't match standard error of fixed effect of Z either...
    # # TODO: match up SE with the formula
    # sigma_hat <- mod %>% resid() %>% sd()
    # sigma_hat <- resid(mod)^2 %>% sum() / (J*n - 2)
    # sigma_hat * sqrt(1 / (0.5 * 0.5 * J * n))

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


if ( FALSE ) {
    os <- one_sim( 20, 20, 0.8 )
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

# store power_sim() results in df_sim as list column
df_sim$data = pmap( df_sim, power_sim, NUMSIM = 100 )


### Q: what is pvalue_one?
# - idea: df is high, so pvalue = pnorm(t)
# - A: so pvalue_one is the p-value for the one-tailed test of whether the true effect is greater than zero

# unnest df_sim & record pvalue_one per site
hits = df_sim %>% 
    rename( n_bar = n ) %>%
    unnest( cols=data ) %>%
    mutate( ATE_bin = cut( ATE, 10 ),
            pvalue_one = pnorm( -t ))

# per simulation setting & ATE_bin: 
#  how often do we reject the null, based on pvalue_one?
agg_hits =  hits %>%
    group_by( tau, n_bar, J, ATE_bin ) %>%
    summarise( ATE = mean( ATE ),
               power = mean( pvalue_one <= 0.05 ) )
agg_hits

### Q: why is "power" low for strongly negative ATEs?
# A: we're only looking at the one-tailed test

# plot power vs. ATE size
ggplot( agg_hits, aes( x=ATE, y=power, col = as.factor(tau) ) ) +
    facet_wrap( ~ n_bar, labeller=label_both ) +
    geom_line() + 
    geom_point() +
    geom_hline( yintercept = 0.8 )
ggsave("plots/power_plot.png")

# issue: for the "small effect" bin, 
#  increasing tau also happens to increase the true average ATE...
agg_hits %>%
    filter(ATE_bin == levels(agg_hits$ATE_bin)[5]) %>%
    arrange(n_bar)


#####
# visualizing power for sites with ATE=0.2
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- hits %>% filter( abs( ATE - 0.20 ) < 0.02 ) %>%
    group_by( tau, n_bar, J ) %>%
    summarise( n = n(),
               ATE_bar = mean(ATE),
               power = mean( pvalue_one <= 0.05 ) )
head( circ20 )

ggplot( circ20, aes( tau, power, col=as.factor( n_bar) ) ) +
    geom_point() +
    geom_line() +
    labs( title = "Power to detect a site with 0.20 ATE", 
          x = "Average ATE across sites" )
ggsave("plots/power_plot_ATE02.png")


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
hits %>% 
    filter( abs( ATE - 0.20 ) < 0.02 ) %>%
    ggplot( aes( x=ATE_hat ) ) +
        facet_grid( tau ~ n_bar, labeller = label_both ) +
        geom_histogram() +
        geom_vline( xintercept = 0.2, col="red" )
ggsave("plots/power_plot_ATE02_hist.png", width=9, height=4)


if ( FALSE ) {
df_sim_res <- df_sim %>%
    rowwise() %>%
    mutate(power = power_sim(n, tau))

# visualize results
ggplot(df_sim_res, aes(x=n, y=power, color=tau)) +
    geom_point() +
    geom_line(aes(group=tau)) +
    geom_hline(aes(yintercept=0.8), color="red", lty="dashed") +
    coord_cartesian(ylim=c(0,1))
}


#####
# visualizing power 
#####




# Note on "effect size" when we have an effect distribution:
#  - each site has a (random) effect size
#     - the simple idea is to aggregate "effect sizes" across simulations
#       (within effect-size bins) to compute power for that effect-size bin
#
#  - wrinkle: detection of an effect for each site *technically* depends on all of
#    the other sites in a given simulation (via the s.e. estimate)
#     - so just averaging all sites within the same effect-size bin (across simulations)
#       doesn't tell the whole story, but it's close enough, I think




