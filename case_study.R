
# power study

# idea: we have actual sample sizes, and number of sites.
#  - we want to know: how precise can we get for individual site-level tx effects?

# note: model is wrong, since p_ij is not bounded in (0,1)
#  - but it's much easier to follow along with, so let's stick with it


require(tidyverse)
require(mvtnorm)
require(rstan)


# simulate data -----------------------------------------------------------

sim_case_study <- function(site_sizes,
                           tau, sig_tau,
                           alpha, sig_alpha,
                           rho) {
  # simulate site-level parameters
  site_pars <- rmvnorm(n = length(site_sizes),
                       mean = c(tau, alpha),
                       sigma = matrix(c(sig_tau^2, rho*sig_tau*sig_alpha, 
                                        rho*sig_tau*sig_alpha, sig_alpha^2),
                                      ncol = 2)) %>% 
    as_tibble(.name_repair = "unique") %>%
    rename(tau_j = ...1,
           alpha_j = ...2) %>% 
    mutate(sid = 1:length(site_sizes), 
           .before=tau_j)
  
  # confirm that all probabilities are within (0,1)
  site_pars %>% 
    mutate(valid = tau_j+alpha_j >= 0 & tau_j+alpha_j <= 1) %>% 
    summarize(valid = min(valid)) %>% 
    pull(valid) %>% 
    if(!.) {stop("Bad simulation (p not in (0,1)): try again")}
  
  # simulate individual-level data
  sim <- tibble(
    sid = rep(1:length(site_sizes), site_sizes)
  ) %>% 
    left_join(site_pars, by="sid") %>% 
    mutate(z = purrr::rbernoulli(n(), p = 0.5),
           y = purrr::rbernoulli(n(), p = alpha_j + tau_j*z))
  
  # make site-level summaries
  sim_sum <- sim %>% 
    group_by(sid) %>% 
    summarize(n1 = sum(z),
              n0 = sum(1-z),
              p1 = sum(z*y)/n1,
              p0 = sum((1-z)*y)/n0,
              tau_j = first(tau_j),
              tau_j_hat = p1-p0,
              se_j = sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0))
  
  return(sim_sum)
}

if (F) {
  # set site sample sizes
  site_sizes <- c(551, 412, 343, 173, 464, 544, 499, 396, 197, 116)
  
  # set data-generating parameters (note: no ICC)
  tau <- 0.2
  sig_tau <- 0.1
  alpha <- 0.175
  sig_alpha <- 0.01
  rho <- 0
  
  sim_case_study(site_sizes,
                 tau, sig_tau,
                 alpha, sig_alpha,
                 rho)
}


# run models -----------------------------------------------------------

run_t_test <- function(sdat) {
  sdat %>% 
    select(sid, tau_j, tau_j_hat, se_j) %>% 
    mutate(q5 = tau_j_hat + qnorm(0.05) * se_j,
           q95 = tau_j_hat + qnorm(0.95) * se_j,
           method = "single")
}

run_mlm <- function(sdat) {
  # make dataset for bayesian models
  stan_list <- list(
    J = length(site_sizes),
    site_mn_obs = sdat$tau_j_hat,
    site_sd_obs = sdat$se_j)
  
  # load model, if not already loaded (loading models takes time)
  if (!exists("mod_norm")) {
    # global assignment to avoid garbage collection and subsequent recompiling
    mod_norm <<- stan_model("Stan/dp_normal_reparam.stan",
                            model_name = "norm",
                            auto_write = T)
  }
  
  # fit model
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
  site_effects_norm <- samples_norm$site_mn
  if (FALSE) {
    site_effects_norm %>%
      as_tibble() %>%
      pivot_longer(everything()) %>%
      ggplot() +
      geom_histogram(aes(x=value)) +
      facet_wrap(~name)
  }
  
  sdat %>% 
    select(sid, tau_j) %>% 
    mutate(tau_j_hat = apply(site_effects_norm, 2, mean),
           se_j = apply(site_effects_norm, 2, sd),
           q5   = apply(site_effects_norm, 2, function(x) quantile(x, 0.05)),
           q95  = apply(site_effects_norm, 2, function(x) quantile(x, 0.95)),
           method = "MLM")
}

if (F) {
  sdat <- sim_case_study(
    site_sizes,
    tau, sig_tau,
    alpha, sig_alpha,
    rho)
  
  run_t_test(sdat)
  run_mlm(sdat)
}












# attempting to simulate data the right way
#  - sorta following https://www.barelysignificant.com/post/icc/
if (F) {
  logit <- function(x) {
    log(x/(1-x))
  }
  dlogit <- function(x) {
    1/x + 1/(1-x)
  }
  ddlogit <- function(x) {
    1/(1-x)^2 - 1/x^2
  }
  
  inv_logit <- function(x) {
    1/(1+exp(-x))
  }
  
  # simulate alpha_j: 
  #  - GOAL: p is "...base rates will likely be in the 15-20% range..."
  #  - IDEA: set mu_alpha, sig_alpha such that p_0j ~ N(0.175, 0.01^2),
  #     i.e., approximate logit(p_0j) with a N(mu_alpha, sig_alpha^2)
  #  - METHOD: https://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables
  mu_p <- 0.175
  sig_p <- 0.01
  mu_alpha <- logit(mu_p) + ddlogit(mu_p)*sig_p^2/2
  sig_alpha <- sqrt(dlogit(mu_p)^2*sig_p^2 - ddlogit(mu_p)^2*sig_p^4/4)
  
  # simulate tau_j: 
  #  - GOAL: txeff is "..something in the 2-4 pp range..."
  mu_txeff <- 0.2
  sig_txeff <- 0.1
  rho <- 0
  
  # check that approximations look okay
  if (F) {
    # alpha dist approximation
    x <- rnorm(10000, mean=mu_p, sd=sig_p)
    x_approx <- rnorm(10000, mean=mu_alpha, sd=sig_alpha)
    tibble(
      x = logit(x),
      x_approx = x_approx
    ) %>% 
      pivot_longer(everything()) %>% 
      ggplot(aes(x=value, color=name)) +
      geom_density()
    
    # tau dist approximation
  }
  
  # from method 1 in: https://www.barelysignificant.com/post/icc/
  #  - note: icc is function of sig_alpha, so it's not user-specified
  icc <- sig_alpha^2/(sig_alpha^2 + pi^2/3)
}