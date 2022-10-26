
# same as case_study, but with a gamma tx effect distribution

require(tidyverse)
require(mvtnorm)
require(rstan)

ON_CLUSTER <- T


# simulate data -----------------------------------------------------------

# simulate treatment effects from a gamma distribution
sim_case_study2 <- function(site_sizes,
                            tau, b,
                            alpha, sig_alpha) {
  
  site_pars <- tibble(
    sid = 1:length(site_sizes),
    tau_j = rgamma(length(site_sizes), shape=tau*b, rate=b),
    alpha_j = rnorm(length(site_sizes), mean=alpha, sd=sig_alpha),
  )
  
  # simulate individual-level data
  sim <- tibble(
    sid = rep(1:length(site_sizes), site_sizes)
  ) %>% 
    left_join(site_pars, by="sid") %>% 
    mutate(z = purrr::rbernoulli(n(), p = 0.5),
           p = alpha_j + tau_j*z,
           p = ifelse(p <= 0, 0, p),
           p = ifelse(p >= 1, 1, p),
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
  site_sizes <- c(551, 928, 895, 1008, 309)
  
  
  # set data-generating parameters (note: no ICC)
  tau <- 0.03
  b <- 50
  alpha <- 0.175
  sig_alpha <- 0.01
  
  sim_case_study2(site_sizes,
                 tau, b,
                 alpha, sig_alpha)
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
    tau_j_hat = sdat$tau_j_hat,
    se_j = sdat$se_j)
  
  # load model, if not already loaded (loading models takes time)
  if (!exists("mod_norm")) {
    # global assignment to avoid garbage collection and subsequent recompiling
    mod_norm <<- stan_model(# "Stan/dp_normal_reparam.stan",
      "Stan/case_study_model.stan",
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
  
  # diagnose model
  if (F) {
    # # shinythemes dependency is broken for whatever reason....
    # require(shinystan)
    # launch_shinystan(fit_norm)
    
    stan_diag(fit_norm)
    print(fit_norm)
    
    require(bayesplot)
    mcmc_pairs(fit_norm, pars=c("tau", "sig_tau"))
  }
  
  samples_norm <- rstan::extract(fit_norm)
  site_effects_norm <- samples_norm$tau_j
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

run_case_study2 <- function(site_sizes,
                            tau, b,
                            alpha, sig_alpha) {
  sdat <- sim_case_study2(
    site_sizes,
    tau, b,
    alpha, sig_alpha)
  
  run_t_test(sdat) %>% 
    bind_rows(run_mlm(sdat))
}

if (F) {
  # set site sample sizes
  site_sizes <- c(551, 412, 343, 173, 464, 544, 499, 396, 197, 116)
  site_sizes <- c(551, 928, 895, 1008, 309)
  
  
  # set data-generating parameters (note: no ICC)
  tau <- 0.02
  b <- 50
  alpha <- 0.175
  sig_alpha <- 0.01
  
  sdat <- sim_case_study2(
    site_sizes,
    tau, b,
    alpha, sig_alpha)
  run_t_test(sdat)
  
  set.seed(90210)
  run_mlm(sdat)
  
  run_case_study(site_sizes,
                 tau, sig_tau,
                 alpha, sig_alpha,
                 rho)
}



# run simulation ----------------------------------------------------------

# set site sample sizes
site_sizes <- c(551, 412, 343, 173, 464, 544, 499, 396, 197, 116)
# site_sizes <- c(551, 928, 895, 1008, 309)

# expand parameter grid
df_sim <- expand_grid(
  tau   = c(0.03),
  b     = c(50, 150, 250),
  alpha = c(0.175),
  sig_alpha = c(0.01)
)

if(ON_CLUSTER) {
  require(uuid)
  require(glue)
  
  for (i in 1:10) {
    FNAME <- glue("case_study_results/res2_{str_sub(UUIDgenerate(), 1, 7)}.csv")
    # FNAME <- glue("case_study_results/stateres2_{str_sub(UUIDgenerate(), 1, 7)}.csv")
    res <- df_sim %>% 
      rowwise() %>% 
      mutate(res = list(run_case_study2(site_sizes,
                                        tau, b,
                                        alpha, sig_alpha))) %>% 
      unnest(res)
    res %>% 
      mutate(runID = i, .before=tau) %>% 
      write_csv(FNAME)
  }
} else { # running locally
  for (i in 1:1000) {
    res <- df_sim %>% 
      rowwise() %>% 
      mutate(res = list(run_case_study(site_sizes,
                                       tau, sig_tau,
                                       alpha, sig_alpha,
                                       rho))) %>% 
      unnest(res)
    
    # store results
    FNAME <- "case_study/case_study_results.csv"
    if (file.exists(FNAME)) {
      res %>%
        mutate(runID = counter, .before=tau) %>%
        write_csv(FNAME, append=T)
      
      counter <- counter+1
    } else {
      counter <- 1
      res %>% 
        mutate(runID = counter, .before=tau) %>% 
        write_csv(FNAME)
    }
  }
}



