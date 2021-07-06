
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
                                            tau.11.star = 0.3,  # tx var
                                            gamma.10 = tau,  # ATE
                                            return.sites = TRUE) 
    head( sdat )
    sdat$sid = as.character( 1:nrow(sdat) )
    
    # Note: gen_indiv_data() is not in CRAN version of package
    dat = blkvar::generate_individual_data( sdat )
    head( dat )
    
    # Switch to FIRC model someday.
    # Or our bayesian modeling to get posteriors.
    # Currently: RIRC model
    mod = lmer( Yobs ~ 1 + Z + (1+Z|sid), data=dat )

    res = as.data.frame( coef( mod )$sid )
    ses = as.data.frame( se.ranef(mod)$sid )

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
    # n = c(25, 50, 75, 100),
    n = c(25, 100),
    J = 20,
    tau = c(0.01, 0.8)
    # tau = c(0.01, 0.2, 0.5, 0.8)
)

# store power_sim() results in df_sim as list column
df_sim$data = pmap( df_sim, power_sim, NUMSIM = 100 )

# unnest df_sim & record pvalue_one per site
#  - Q: what is pvalue_one???
hits = df_sim %>% 
    rename( n_bar = n ) %>%
    unnest( cols=data ) %>%
    mutate( ATE_bin = cut( ATE, 10 ),
            # reject = pvalue <= 0.05,   # doesn't this make more sense?
            pvalue_one = pnorm( -t ))

# per simulation setting & ATE_bin: 
#  how often do we reject the null, based on pvalue_one?
agg_hits =  hits %>%
    group_by( tau, n_bar, ATE_bin ) %>%
    summarise( ATE = mean( ATE ),
               power = mean( pvalue_one <= 0.05 ) )
# # doesn't this make more sense?
# agg_hits <- hits %>%
#     group_by(tau, n_bar, ATE_bin) %>%
#     summarize(ATE = mean(ATE),
#               power = mean(reject))
agg_hits

# Q: why is "power" low for strongly negative ATEs?
ggplot( agg_hits, aes( x=ATE, y=power, col = as.factor(tau) ) ) +
    facet_wrap( ~ n_bar ) +
    geom_line() + geom_point() +
    geom_hline( yintercept = 0.8 )


#####
# visualizing power for sites with ATE=0.2
#####

# for sites with ATE=0.2, plot power (across sim settings)
circ20 <- hits %>% filter( abs( ATE - 0.20 ) < 0.02 ) %>%
    group_by( tau, n_bar ) %>%
    summarise( power = mean( pvalue_one <= 0.05 ) )
head( circ20 )

ggplot( circ20, aes( tau, power, col=as.factor( n_bar) ) ) +
    geom_point() +
    geom_line() +
    labs( title = "Power to detect a site with 0.20 ATE", 
          x = "Average ATE across sites" )


# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
hits %>% 
    filter( abs( ATE - 0.20 ) < 0.02 ) %>%
    ggplot( aes( x=ATE_hat ) ) +
        facet_wrap( tau ~ n_bar, labeller = label_both ) +
        geom_histogram() +
        geom_vline( xintercept = 0.2, col="red" )


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

# Note on "effect size" when we have an effect distribution:
#  - each site has a (random) effect size
#     - the simple idea is to aggregate "effect sizes" across simulations
#       (within effect-size bins) to compute power for that effect-size bin
#
#  - wrinkle: detection of an effect for each site *technically* depends on all of
#    the other sites in a given simulation (via the s.e. estimate)
#     - so just averaging all sites within the same effect-size bin (across simulations)
#       doesn't tell the whole story, but it's close enough, I think




