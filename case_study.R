
# power study

# idea: we have actual sample sizes, and number of sites.
#  - we want to know: how precise can we get for individual site-level tx effects?


require(tidyverse)
require(mvtnorm)


# simulate data -----------------------------------------------------------

# set site sample sizes
site_sizes <- c(551, 412, 343, 173, 464, 544, 499, 396, 197, 116)







site_pars <- rmvnorm(n = length(site_sizes),
                     mean = c(tau, mu_alpha),
                     sigma = matrix(c(tx_sd^2, rho*tx_sd*sig_alpha, 
                                      rho*tx_sd*sig_alpha, sig_alpha^2),
                                    ncol = 2)) %>% 
  as_tibble(.name_repair = "unique") %>%
  rename(tau_j = ...1,
         alpha_j = ...2) %>% 
  mutate(sid = 1:length(site_sizes), 
         .before=tau_j) 

# simulate individual-level data
sim <- tibble(
  sid = rep(1:length(site_sizes), site_sizes)
) %>% 
  left_join(site_pars, by="sid") %>% 
  mutate(z = purrr::rbernoulli(n()),
         y = purrr::rbernoulli(n(), p = alpha_j + tau_j*z))

# make site-level summaries
sim_sum <- sim %>% 
  group_by(sid) %>% 
  summarize(n1 = sum(z),
            n0 = sum(1-z),
            p1 = sum(z*y)/n1,
            p0 = sum((1-z)*y)/n0,
            tau_j_hat = p1-p0,
            se_j = sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0))



# ah shoot it's bernoulli.



one_sim <- function(nbar, J, tau, ICC, tx_sd, 
                    variable.n = F,
                    round_sites = 0.05,
                    # round_sites = round(tx_sd*4/20, 3),   # 20 bins for +/- 2 sd
                    NUMSAMP = 2000) {
  
  # if tx_sd is very small, set smaller bins!
  if (tx_sd <= 0.1) {
    round_sites <- 0.02
    if (tx_sd <= 0.05) {
      round_sites <- 0.01
    }
  }
  print(glue("Using bin size of {round_sites}"))
  
  ##### simulate data #####
  
  sdat = blkvar::generate_multilevel_data( n.bar = nbar, 
                                           J = J,
                                           variable.n = variable.n,
                                           tau.11.star = tx_sd^2,   # cross-site tx var
                                           gamma.10 = tau,      # cross-site ATE
                                           ICC = ICC,           # ICC [really var(intercepts)]
                                           rho2.0W = 0,         # covariate has no explanatory power 
                                           rho2.1W = 0,
                                           zero.corr = T,       # don't correlate site intercepts and treatment effects (so treatment group is higher variance)
                                           return.sites = TRUE,
                                           verbose = FALSE) %>%
    mutate(sid = as.character(1:n()),
           beta.1 = round(beta.1/round_sites) * round_sites)
  head(sdat)
  
  # Note: generate_individual_data() is not in CRAN version of package
  # dat = blkvar::generate_individual_data( sdat )
  dat = generate_individual_data(sdat,
                                 sigma2.e = 1-ICC)   # need to specify this! or else sigma2.e = 1...
  head( dat )
  
  ##### run models #####
  
  NAME_VEC <- c("ATEhat", "SE", paste0("q", QUANTILES))
  
  ### run single-site models
  
  res_single <- dat %>%
    group_by(sid) %>%
    dplyr::group_modify(~run_t_test(., NAME_VEC)) %>%
    ungroup() %>%
    mutate(sid = as.character(sid))
  
  if (F) {
    ### run frequentist multilevel models (FIRC and RIRC)
    
    # FIRC model #1 (via generic arm package functions)
    mod_firc = lmer( Yobs ~ 0 + as.factor(sid) + Z + (0+Z|sid), data=dat)
    
    ATEhat_firc <- coef(mod_firc)$sid$Z
    SE_firc1 <- sqrt(se.fixef(mod_firc)["Z"]^2 + se.ranef(mod_firc)$sid[,"Z"]^2)
    
    res_firc1 <- tibble(sid = as.character(1:J),
                        ATEhat = ATEhat_firc,
                        SE = SE_firc1) %>%
      rowwise() %>%
      mutate(q = list(ATEhat + qnorm(QUANTILES) * SE)) %>%
      unnest(cols = q) %>%
      group_by(sid) %>%
      mutate(quantile = QUANTILES) %>%
      ungroup() %>%
      pivot_wider(names_from = quantile, values_from = q, names_prefix = "q")
    names(res_firc1) <- c("sid", paste0(NAME_VEC, "_firc1"))
    
    # FIRC model #2 (via arm package sampling functions)
    
    # grabbing arm package samples
    sim_firc <- sim(mod_firc, n.sims = NUMSAMP)
    fixef_samps_firc <- coef(sim_firc)[["fixef"]][,"Z"]
    ranef_samps_firc <- coef(sim_firc)[["ranef"]]$sid[,,"Z"]
    
    res_firc2 <- tibble(
      sid = as.character(1:J),
      ATEhat = ATEhat_firc,
      SE = apply(ranef_samps_firc, 2, function(x) sd(x + fixef_samps_firc))
    ) %>%
      cbind(t(apply(ranef_samps_firc, 2, 
                    function(x) quantile(x + fixef_samps_firc, probs=QUANTILES)))) %>%
      as_tibble(.name_repair = "minimal")
    names(res_firc2) <- c("sid", paste0(NAME_VEC, "_firc2"))
    
    # RIRC model #1 (via generic arm package functions)
    mod_rirc = lmer( Yobs ~ 1 + Z + (1+Z|sid), data=dat )
    
    ATEhat_rirc <- coef(mod_rirc)$sid$Z
    SE_rirc1 <- sqrt(se.fixef(mod_rirc)["Z"]^2 + se.ranef(mod_rirc)$sid[,"Z"]^2)
    
    res_rirc1 <- tibble(sid = as.character(1:J),
                        ATEhat = ATEhat_rirc,
                        SE = SE_rirc1) %>%
      rowwise() %>%
      mutate(q = list(ATEhat + qnorm(QUANTILES) * SE)) %>%
      unnest(cols = q) %>%
      group_by(sid) %>%
      mutate(quantile = QUANTILES) %>%
      ungroup() %>%
      pivot_wider(names_from = quantile, values_from = q, 
                  names_prefix = "q")
    names(res_rirc1) <- c("sid", paste0(NAME_VEC, "_rirc1"))
    
    # RIRC model #2 (via arm package sampling functions)
    
    # grabbing arm package samples
    sim_rirc <- sim(mod_rirc, n.sims = NUMSAMP)
    fixef_samps_rirc <- coef(sim_rirc)[["fixef"]][,"Z"]
    ranef_samps_rirc <- coef(sim_rirc)[["ranef"]]$sid[,,"Z"]
    
    res_rirc2 <- tibble(
      sid = as.character(1:J),
      ATEhat_rirc2 = ATEhat_rirc,
      SE_rirc2 = apply(ranef_samps_rirc, 2, function(x) sd(x + fixef_samps_rirc))
    ) %>%
      cbind(t(apply(ranef_samps_rirc, 2, 
                    function(x) quantile(x + fixef_samps_rirc, probs=QUANTILES)))) %>%
      as_tibble(.name_repair = "minimal")
    names(res_rirc2) <- c("sid", paste0(NAME_VEC, "_rirc2"))
    
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
  
  fit_norm <- sampling(mod_norm,
                       data = stan_list,
                       iter = NUMSAMP,
                       chains = 4,
                       control = list(max_treedepth = 12,
                                      adapt_delta = 0.95),
                       verbose = F, 
                       show_messages = F, 
                       refresh = 0)
  
  # wrap in error detection
  num_divergences <- 0
  more_samples <- F
  withCallingHandlers(
    warning = function(cnd) {
      # browser()
      
      if (str_detect(cnd$message, "Bulk Effective Samples Size")) {
        print("Bulk ESS too low: doubling number of samples")
        more_samples <<- T
        fit_norm <- sampling(mod_norm,
                             data = stan_list,
                             iter = NUMSAMP,
                             chains = 4,
                             control = list(max_treedepth = 12,
                                            adapt_delta = 0.95),
                             verbose = F, 
                             show_messages = F, 
                             refresh = 0)
      }
    },
    { # sample!
      fit_norm <- sampling(mod_norm,
                           data = stan_list,
                           iter = NUMSAMP,
                           chains = 4,
                           control = list(max_treedepth = 12,
                                          adapt_delta = 0.95),
                           verbose = F, 
                           show_messages = F, 
                           refresh = 0) }
  )
  
  samples_norm <- rstan::extract(fit_norm)
  names(samples_norm)
  
  # browser()
  
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
    SE_bayesnorm = apply(site_effects_norm, 2, sd)) %>%
    cbind(t(apply(site_effects_norm, 2, 
                  function(x) quantile(x, probs=QUANTILES)))) %>%
    as_tibble(.name_repair = "minimal")
  names(res_bayesnorm) <- c("sid", paste0(NAME_VEC, "_bayesnorm"))
  
  ##### compile results #####
  
  if (F) {
    # compile results for overall tau
    mod_pooled <- lm(Yobs ~ 0 + as.factor(sid) + Z, data=dat)   # use pooled model
    res_overall <- tibble(
      method = c("FIRC", "RIRC", "Bayes", "Single"),
      ATEhat = c(fixef(mod_firc)["Z"], fixef(mod_rirc)["Z"],
                 mean(samples_norm$pop_mn), coef(mod_pooled)["Z"]),
      SE     = c(se.fixef(mod_firc)["Z"], se.fixef(mod_rirc)["Z"],
                 sd(samples_norm$pop_mn), summary(mod_pooled)$coef["Z",2]))
  }
  
  # browser()
  
  # compile results for single-site tau_js
  res_sites <- as_tibble(sdat) %>%
    dplyr::select(sid, n, ATE = beta.1) %>%
    left_join(res_single, by="sid") %>%
    # left_join(res_firc1, by="sid") %>%
    # left_join(res_firc2, by="sid") %>%
    # left_join(res_rirc1, by="sid") %>%
    # left_join(res_rirc2, by="sid") %>%
    left_join(res_bayesnorm, by="sid") %>%
    mutate(# is_singular_firc = isSingular(mod_firc), # is FIRC model singular?
      # is_singular_rirc = isSingular(mod_rirc), # is RIRC model singular?
      ESS_low = more_samples)                  # low ESS warning from rstan?
  
  # return(list(overall = res_overall, sites = res_sites))
  return(list(sites = res_sites))
}















# attempting to simulate data the right way
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