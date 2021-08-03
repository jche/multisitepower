
#####
# Power analysis for individual site effects in multisite trial,
# comparing MLMs to just using each single site
# 
# uses blkvar package: devtools::install_github("https://github.com/lmiratrix/blkvar")
#####

library( arm )   # note: loads MASS package
# devtools::install_github("lmiratrix/blkvar")
library( blkvar )

require(tidyverse)
require(glue)
require(tictoc)


#####
# simulation functions
#####

# given a df for a single site, run a one-sided t-test
run_t_test <- function(df) {
  
  df %>%
    summarize(t = list(t.test(Y1[Z==1], Y0[Z==0], alternative = "greater"))) %>%
    mutate(ATE_hat_single = t[[1]]$estimate["mean of x"] - t[[1]]$estimate["mean of y"],
           SE_single = t[[1]]$stderr,
           t_single = t[[1]]$statistic,
           pvalue_single = t[[1]]$p.value) %>%
    select(-t)
}

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
one_sim <- function(n, J, tau, ICC, round_sites = 0.05) {
  
  # browser()
  
  sdat = blkvar::generate_multilevel_data( n.bar=n, J=J,
                                           variable.n = FALSE,
                                           tau.11.star = 0.3,   # cross-site tx var
                                           gamma.10 = tau,      # cross-site ATE
                                           ICC = ICC,           # ICC [really var(intercepts)!]
                                           # ICC = var_alpha,
                                           # sigma2.e = var_eps,
                                           rho2.0W = 0,         # covariate has no explanatory power 
                                           rho2.1W = 0,
                                           zero.corr = T,       # don't correlate site intercepts and treatment effects (so treatment group is higher variance)
                                           return.sites = TRUE,
                                           verbose = FALSE) %>%
    mutate(sid = as.character(1:n()),
           beta.1 = round(beta.1/round_sites) * round_sites)    # round to nearest round_sites
  head(sdat)
  
  # Note: generate_individual_data() is not in CRAN version of package
  # dat = blkvar::generate_individual_data( sdat )
  dat = generate_individual_data(sdat,
                                 sigma2.e = 1-ICC)   # need to specify this! or else sigma2.e = 1...
  head( dat )
  
  # conclusion: "ICC" argument is just Var(\alpha_i), NOT the ICC itself!
  if (FALSE) {
    var(sdat$beta.0)
    var(sdat$u0)   # what is u0? has to do with covariate, I think...?
    
    # cross-site variation: 1 for control group, 1.3 for treatment group
    dat %>%
      group_by(Z) %>%
      summarize(var = var(Yobs))
    
    # within-site variation: always 1 if we don't specify sigma2.e in gen_indiv_data()
    #  - 1-ICC if we correctly specify it!
    dat %>%
      group_by(Z, sid) %>%
      summarize(var = var(Yobs)) %>%
      summarize(mn_var = mean(var))
  }
  
  # run single-site models!
  res_single <- dat %>%
    group_by(sid) %>%
    dplyr::group_modify(~run_t_test(.))
    # summarize(ttest = list(run_t_test()))
  
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
  
  res %>%
    left_join(res_single, by="sid")
}


if ( F ) {
  os <- one_sim( n=200, J=200, tau=0.2, ICC=0.3)
  mean( os$ATE )
  mean( os$ATE_hat )
  ggplot( os, aes( ATE, ATE_hat ) ) +
    geom_point() +
    geom_abline( slope=1, intercept= 0 )
}

# function: (# obs, effect size) => (power)
#  - runs one_sim NUMSIM times, so we can aggregate across runIDs to get the power
power_sim <- function(n, J, tau, ICC, NUMSIM = 250) {
  cat(glue("Working on n = {n}, tau = {tau}, ICC = {ICC}\n\n"))
  rs = tibble( runID = 1:NUMSIM )
  rs$data = map( rs$runID, ~one_sim(n, J, tau, ICC))
  rs = unnest( rs, cols = data )
  rs
}

if ( FALSE ) {
  power_sim( 20, 20, 0.2, 0, 3 )
}


#####
# run power simulation
#####

# run power simulation:
#  - n = number of observations
#  - tau = true effect size
#  - (real power sim would have more knobs: J, site.size, ICC/variances, etc.)
df_sim <- expand_grid(
  n   = c(25, 50, 75, 100),
  J   = c(20),
  # ICC = c(0, 0.5, 1, 2),   # equiv to ICC = c(0, 0.33, 0.5, 0.67)
  ICC = c(0, 0.3, 0.6, 0.9),
  tau = c(0.01, 0.2, 0.5, 0.8)
)

# run simulation: store power_sim() results in df_sim as list column
tic()
df_sim <- df_sim %>%
  rowwise() %>%
  mutate(data = list(power_sim(n, J, tau, ICC, NUMSIM = 400)))
toc()

# raw results: unnest df_sim & record pvalue_one per site
hits = df_sim %>% 
  rename( n_bar = n ) %>%
  unnest( cols=data ) %>%
  mutate( pvalue_one = pnorm( -t ))


#####
# save results
#####

FNAME <- "sim_results_ICC_fixed"
fname <- glue("results/", FNAME, ".csv")

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


