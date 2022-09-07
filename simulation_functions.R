

QUANTILES <- c(0.05, 0.1, 0.95)


#####
# main simulation functions
#####

# function: (# obs, effect size) => (reject null? T/F)
#' 
#'
#' @param n number of observations per site (constant, for now)
#' @param J number of sites
#' @param tau overall ATE
#' @param ICC 
#'
#' @return tibble with J rows
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


if ( F ) {
  os <- one_sim( n=25, J=50, tau=0.01, ICC=0, tx_sd=0.3)
  browser()
  
  ATEhats <- os %>%
    pivot_longer(contains("ATEhat"), names_to = "method", values_to = "ATEhat") %>%
    dplyr::select(sid, method, ATEhat) %>%
    mutate(method = str_sub(method, 8))
  SEs <- os %>%
    # mutate(SE_firc = sqrt(SE_firc_fixed^2 + SE_firc_rand^2),
    #        SE_rirc = sqrt(SE_rirc_fixed^2 + SE_rirc_rand^2)) %>%
    pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
    dplyr::select(sid, method, SE) %>%
    mutate(method = str_sub(method, 4))
  
  tidy_results <- os %>%
    dplyr::select(sid, n, ATE) %>%
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
power_sim <- function(nbar, J, tau, ICC, tx_sd, 
                      variable.n = F,
                      NUMSIM = 250, 
                      WRITE_CSV = F, 
                      fname = NULL) {
  cat(glue("Working on nbar = {nbar}, J = {J}, ICC = {ICC}, tau = {tau}, tx_sd = {tx_sd}\n\n"))
  rs = tibble( runID = 1:NUMSIM )
  rs$data = map( rs$runID, ~one_sim(nbar, J, tau, ICC, tx_sd, variable.n))
  rs = unnest( rs, cols = data )
  
  if(WRITE_CSV) {
    hits <- rs %>%
      mutate(nbar = nbar, J=J, ICC=ICC, tau=tau, tx_sd=tx_sd,
             .before = 1)
    
    if (!file.exists(fname)) {
      write_csv(hits, fname)
    } else {
      write_csv(hits, fname, append=T)
    }
  }
  
  rs
}

if ( F ) {
  power_sim( nbar=20, J=20, tau=0.2, ICC=0, tx_sd=0, NUMSIM=3 )
}


#####
# simulation helper functions
#####

# given a df for a single site, run a one-sided t-test
run_t_test <- function(df, NAME_VEC = NAME_VEC) {

  # df %>%
  #   summarize(t = list(t.test(Y1[Z==1], Y0[Z==0]))) %>%
  #   mutate(ATEhat_single = t[[1]]$estimate["mean of x"] - t[[1]]$estimate["mean of y"],
  #          SE_single = t[[1]]$stderr,     # standard SE, but of a t dist.
  #          df_single = t[[1]]$parameter   # df of the t distribution
  #   ) %>%
  #   dplyr::select(-t)
  
  t_test <- with(df, t.test(Y1[Z==1], Y0[Z==0]))
  
  ATEhat_single <- t_test$estimate["mean of x"] - t_test$estimate["mean of y"]
  SE_single     <- t_test$stderr
  df_single     <- t_test$parameter
  quants        <- ATEhat_single + qt(QUANTILES, df=df_single) * SE_single
  
  res <- c(ATEhat_single, SE_single, quants) %>%
    t() %>%
    as_tibble(.name_repair = "minimal")
  names(res) <- paste0(NAME_VEC, "_single")
  
  return(res)
}

# given the full df with individual observations,
# make df with site-level summaries for rstan functions
make_site_summaries <- function( df ) {

  df %>%
    group_by(sid, Z) %>%
    summarize(ybar = mean(Yobs),
              n = n(),
              V = var(Yobs),
              .groups = "drop_last") %>%
    pivot_wider(names_from = "Z", 
                values_from = ybar:V, 
                names_sep = ".") %>%
    ungroup() %>%
    mutate(pool.se2 = sum( V.1 * (n.1-1) + V.0 * (n.0-1) ) / sum( n.1 + n.0 - 2 ),
           tau.hat = ybar.1 - ybar.0,
           SE = sqrt(pool.se2 / n.1 + pool.se2 / n.0))
}

