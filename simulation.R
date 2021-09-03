
#####
# Power analysis for individual site effects in multisite trial,
# comparing MLMs to just using each single site
# 
# uses blkvar package: devtools::install_github("https://github.com/lmiratrix/blkvar")
#####

library( arm )   # bayesian-type functions for lmer, note this loads MASS package
# devtools::install_github("lmiratrix/blkvar")
library( blkvar )
library(rstan)   # bayesian models!

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
    mutate(ATEhat_single = t[[1]]$estimate["mean of x"] - t[[1]]$estimate["mean of y"],
           SE_single = t[[1]]$stderr
           # t_single = t[[1]]$statistic,
           # pvalue_single = t[[1]]$p.value
           ) %>%
    select(-t)
}

# given the full df with individual observations,
# make df with site-level summaries for rstan functions
make_site_summaries <- function( df ) {
  df %>%
    group_by(sid, Z) %>%
    summarize(ybar = mean(Yobs),
              n = n(),
              V = var(Yobs)) %>%
    pivot_wider(names_from = "Z", 
                values_from = ybar:V, 
                names_sep = ".") %>%
    ungroup() %>%
    mutate(pool.se2 = sum( V.1 * (n.1-1) + V.0 * (n.0-1) ) / sum( n.1 + n.0 - 2 ),
           tau.hat = ybar.1 - ybar.0,
           SE = sqrt(pool.se2 / n.1 + pool.se2 / n.0))
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
one_sim <- function(n, J, tau, ICC, tx_var, round_sites = 0.05) {
  
  ##### simulate data #####
  
  sdat = blkvar::generate_multilevel_data( n.bar=n, J=J,
                                           # variable.n = FALSE,
                                           variable.n = TRUE,
                                           tau.11.star = tx_var,   # cross-site tx var
                                           gamma.10 = tau,      # cross-site ATE
                                           ICC = ICC,           # ICC [really var(intercepts)]
                                           rho2.0W = 0,         # covariate has no explanatory power 
                                           rho2.1W = 0,
                                           zero.corr = T,       # don't correlate site intercepts and treatment effects (so treatment group is higher variance)
                                           return.sites = TRUE,
                                           verbose = FALSE) %>%
    mutate(sid = as.character(1:n()),
           beta.1 = round(beta.1/round_sites) * round_sites)    # round to nearest round_sites
  
  # manually add in a site with ATE = 0.2
  if (tx_var == 0) {
    sdat$beta.1[1] <- 0.2   # idea: set treatment effect for site 1 to 0.2, so there is at least one site with ATE=0.2 for our visualizations
  }
  head(sdat)
  
  # Note: generate_individual_data() is not in CRAN version of package
  # dat = blkvar::generate_individual_data( sdat )
  dat = generate_individual_data(sdat,
                                 sigma2.e = 1-ICC)   # need to specify this! or else sigma2.e = 1...
  head( dat )
  
  
  ##### run models #####
  
  ### run single-site models
  
  res_single <- dat %>%
    group_by(sid) %>%
    dplyr::group_modify(~run_t_test(.)) %>%
    ungroup() %>%
    mutate(sid = as.character(sid))
  
  ### run frequentist multilevel models (FIRC and RIRC)

  # FIRC model
  #  - note: no intercept, so site #1 is the "intercept"
  mod_firc = lmer( Yobs ~ 0 + as.factor(sid) + Z + (0+Z|sid), data=dat)
  res_firc <- tibble(
    sid = as.character(1:J),
    ATEhat_firc = coef(mod_firc)$sid$Z,
    SE_firc_fixed = se.fixef(mod_firc)["Z"],
    SE_firc_rand = se.ranef(mod_firc)$sid[,"Z"]
  )

  # RIRC model
  mod_rirc = lmer( Yobs ~ 1 + Z + (1+Z|sid), data=dat )
  res_rirc <- tibble(
    sid = as.character(1:J),
    ATEhat_rirc = coef(mod_rirc)$sid$Z,
    SE_rirc_fixed = se.fixef(mod_rirc)["Z"],
    SE_rirc_rand = se.ranef(mod_rirc)$sid[,"Z"]
  )
  
  
  
  sim_firc <- sim(mod_firc, n.sims = 5000)
  
  browser()
  
  if (FALSE) {
    
    # se.ranef uses:
    object <- mod_firc
    sqrt(attr( ranef( object, condVar = TRUE )[[1]], "postVar" ))
    se.ranef(mod_firc)
    #  this is square root of the "conditional variance-covariance matrices of the random effects"
    #   - random effects shouldn't have covariances, right? so it's a diagonal matrix
    #   - what does the "conditional" piece mean?
    
    # using the model built-in function gives lower standard errors...?
    #  - worse for lower standard errors (larger sites)...
    tibble(
      # using_samples = apply(coef(sim_firc)[["ranef"]]$sid[,,"Z"], 2, function(x) sd(x + coef(sim_firc)[["fixef"]][,"Z"])),   # our actual standard errors that we want!
      using_samples = apply(coef(sim_firc)[["ranef"]]$sid[,,"Z"], 2, sd),
      # using_model = se.ranef(mod_firc)$sid[,"Z"]
      using_model = se.ranef(mod_firc)$sid[,"Z"] + se.fixef(mod_firc)["Z"]
    ) %>%
      ggplot() +
      geom_point(aes(x=using_model, y=using_samples)) +
      geom_abline() +
      coord_cartesian(xlim = c(0.05, 0.25), ylim = c(0.1, 0.3))
    
    # Q: how to recover model-function's SEs for random effects?
    #  - we match exactly for fixed effects
    #  - we can't recover what se.ranef outputs by using our samples...
    
    # the fixed effects have the same standard errors...
    tibble(using_model = se.fixef(mod_firc),
           using_samples =   coef(sim_firc)[["fixef"]] %>% apply(., 2, sd)) %>%
      ggplot() +
      geom_point(aes(x=using_model, y=using_samples)) +
      geom_abline()
  }
  
  
  ### run bayesian multilevel models
  
  # make dataset for bayesian models
  stan_df <- make_site_summaries(dat)
  stan_list <- list(
    J = J,
    site_mn_obs = stan_df$tau.hat,
    site_sd_obs = stan_df$SE)
  
  ## normal bayesian model

  # load model, if not already loaded (loading models takes time)
  if (!exists("mod_norm")) {
    # global assignment to avoid garbage collection and subsequent recompiling
    mod_norm <<- stan_model("Stan/dp_normal_reparam.stan",
                            model_name = "norm",
                            auto_write = T)
  }
  
  # sample!
  fit_norm <- sampling(mod_norm,
                       data = stan_list,
                       iter = 2000,
                       chains = 4,
                       control = list(max_treedepth = 12,
                                      adapt_delta = 0.95),
                       verbose = F, 
                       show_messages = F, 
                       refresh = 0)
  samples_norm <- rstan::extract(fit_norm)
  names(samples_norm)
  
  # get site-effect estimates from samples
  site_effects_norm <- samples_norm$site_mn
  if (FALSE) {
    site_effects_norm %>%
      as_tibble() %>%
      pivot_longer(everything()) %>%
    ggplot() +
      geom_histogram(aes(x=value)) +
      facet_wrap(~name)
  }
  res_bayesnorm <- tibble(
    sid = as.character(1:J),
    ATEhat_bayesnorm = apply(site_effects_norm, 2, mean),
    SE_bayesnorm = apply(site_effects_norm, 2, sd))
  
  
  ##### compile results #####
  
  as_tibble(sdat) %>%
    select(sid, n, ATE = beta.1) %>%
    left_join(res_single, by="sid") %>%
    left_join(res_firc, by="sid") %>%
    left_join(res_rirc, by="sid") %>%
    left_join(res_bayesnorm, by="sid") %>%
    mutate(is_singular = isSingular(mod_firc))   # indicate if FIRC model fit was singular
}


if ( T ) {
  os <- one_sim( n=100, J=100, tau=5, ICC=0.6, tx_var=0.3)
  
  ATEhats <- os %>%
    pivot_longer(contains("ATEhat"), names_to = "method", values_to = "ATEhat") %>%
    select(sid, method, ATEhat) %>%
    mutate(method = str_sub(method, 8))
  SEs <- os %>%
    mutate(SE_firc = sqrt(SE_firc_fixed^2 + SE_firc_rand^2),
           SE_rirc = sqrt(SE_rirc_fixed^2 + SE_rirc_rand^2)) %>%
    pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
    select(sid, method, SE) %>%
    mutate(method = str_sub(method, 4))
  
  tidy_results <- os %>%
    select(sid, n, ATE) %>%
    left_join(ATEhats, by=c("sid")) %>%
    left_join(SEs, by=c("sid", "method"))
  
  ggplot(tidy_results, aes(x = ATE, color = method)) +
    geom_point(aes(y = ATEhat)) +
    geom_errorbar(aes(ymin = ATEhat - 1.96*SE, 
                      ymax = ATEhat + 1.96*SE)) +
    facet_wrap(~method) +
    geom_abline(slope = 1)
  
  # check one-tailed test
  tidy_results %>%
    mutate(pvalue_one = pnorm(-ATEhat/SE),
           reject = pvalue_one < 0.1)
}

# function: (# obs, effect size) => (power)
#  - runs one_sim NUMSIM times, so we can aggregate across runIDs to get the power
power_sim <- function(n, J, tau, ICC, tx_var, NUMSIM = 250) {
  cat(glue("Working on n = {n}, tau = {tau}, ICC = {ICC}, tx_var = {tx_var}\n\n"))
  rs = tibble( runID = 1:NUMSIM )
  rs$data = map( rs$runID, ~one_sim(n, J, tau, ICC, tx_var))
  rs = unnest( rs, cols = data )
  rs
}

if ( F ) {
  power_sim( n=20, J=20, tau=0.2, ICC=0, tx_var=0, NUMSIM=3 )
}


#####
# run power simulation
#####

# simulation settings
df_sim <- expand_grid(
  n   = c(25, 50, 75),
  J   = c(25, 50, 75),
  ICC = c(0, 0.3, 0.6),
  tau = c(0.01, 0.2, 0.5),
  tx_var = c(0.3)
)

# run simulation: store power_sim() results in df_sim as list column
tic()
df_sim <- df_sim %>%
  rowwise() %>%
  mutate(data = list(power_sim(n, J, tau, ICC, tx_var, NUMSIM = 50)))
toc()

# raw results: unnest df_sim & record pvalue_one per site
hits = df_sim %>% 
  rename( n_bar = n ) %>%
  unnest( cols=data )


#####
# save results
#####

FNAME <- "sim_results_bayes"
fname <- glue("results/", FNAME, ".csv")

if (file.exists(fname)) {
  # ASSUMING that all sim settings are run equally
  max_runID <- read_csv(fname) %>%
    pull(runID) %>%
    max()
  
  hits %>% 
    mutate(runID = runID + max_runID) %>%
    write_csv(fname, append=T)
} else {
  write_csv(hits, fname)
}


