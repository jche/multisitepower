
# small script to demonstrate bayesian modeling of impacts

require(tidyverse)
require(rstan)
require(blkvar)

set.seed(90210)


# generate sample dataset -------------------------------------------------

# input: individual-level data
# output: site-level data
make_site_summaries <- function(df) {
  
  # compute estimated pooled sds for treated and control units
  #  - idea: remove site intercepts and take sd across all sites
  sds <- df %>%
    group_by(sid, Z) %>%
    mutate(Yobs = Yobs-mean(Yobs)) %>%
    ungroup() %>%
    group_by(Z) %>%
    summarize(sigma = sd(Yobs))
  sd1 <- sds %>% 
    filter(Z==1) %>% 
    pull(sigma)
  sd0 <- sds %>% 
    filter(Z==0) %>% 
    pull(sigma)
  
  df %>% 
    group_by(sid, Z) %>% 
    summarize(ybar = mean(Yobs),
              n = n(),
              .groups = "drop_last") %>% 
    pivot_wider(names_from = "Z",
                values_from = c("ybar", "n"),
                names_sep = ".") %>% 
    ungroup() %>% 
    mutate(tau.hat = ybar.1 - ybar.0,
           SE = sqrt(sd1^2 / n.1 + sd0^2 / n.0))   # use pooled sd
}


# generate individual-level data
dat <- blkvar::generate_multilevel_data_no_cov( 
  n.bar = 50, 
  J = 100,
  variable.n = T,
  tau.11.star = 0.2^2, # cross-site tx var
  gamma.10 = 0.2,      # cross-site ATE
  ICC = 0.2,           # ICC, i.e., var(intercepts)
  # rho2.0W = 0,         # covariate has no explanatory power 
  # rho2.1W = 0,
  # zero.corr = T,       # don't correlate site intercepts and treatment effects (so treatment group is higher variance)
  return.sites = FALSE,
  verbose = FALSE) %>%
  tibble()

# make site-level data
sdat <- make_site_summaries(dat)

# store true site-level effects
sdat_true <- dat %>% 
  group_by(sid) %>% 
  summarize(tau_j = mean(Y1) - mean(Y0))


# run rstan model ---------------------------------------------------------

# prepare list input for stan stuff
stan_dat <- list(
  J = nrow(sdat),
  tau_j_hat = sdat$tau.hat,
  se_j = sdat$SE
)

# compile rstan file & load model
mod <- stan_model(
  "Stan/normal_reparam.stan",
  model_name = "norm",
  auto_write = T)

# fit bayesian model
fit <- sampling(mod, data = stan_dat)

# get parameter samples
samples <- rstan::extract(fit)
names(samples)

# summarize posterior mean, sd, quantiles for tau_j parameters
ests_rstan <- sdat %>% 
  select(sid) %>% 
  mutate(tau_j_hat = apply(samples$tau_j, 2, mean),
         se_j = apply(samples$tau_j, 2, sd),
         q5   = apply(samples$tau_j, 2, function(x) quantile(x, 0.05)),
         q95  = apply(samples$tau_j, 2, function(x) quantile(x, 0.95))) %>% 
  left_join(sdat_true, by="sid")


# plot histogram of unadjusted posterior means of sites
ests_rstan %>% 
  pivot_longer(cols = c(tau_j_hat, tau_j)) %>% 
ggplot() +
  geom_density(aes(x=value, color=name))


# plot histogram of unadjusted samples of site effects
ggplot() +
  geom_density(aes(x=tau_j),
               data=sdat_true) +
  geom_density(aes(x=as.numeric(samples$tau_j)),   # all tau_j samples in red
               color = "red")



# rstanarm? ---------------------------------------------------------------

# NOTE: could also use rstanarm package, which is much more accessible
#  - i.e., much less complicated setup, familiar lmer syntax
#  - issue is that I'm not sure how to use it with the site-level summaries,
#    so the fitting is quite slow.

# require(rstanarm)
# fit_rstanarm <- stan_lmer(Yobs ~ 0 + sid + (Z | sid), data=dat)



