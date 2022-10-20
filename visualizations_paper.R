
require(tidyverse)
require(ggthemes)
require(wesanderson)
require(latex2exp)
require(rstan)

pal <- wes_palette("Zissou1", 5, type="continuous")



# analyzing simulation results --------------------------------------------

tidy_results <- read_csv("results_sree/sree_sims.csv") %>% 
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= 0 | q95 <= 0,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05),
         method = case_when(
           method == "bayesnorm" ~ "MLM",
           method == "single" ~ "Single"))


### MDES

tidy_results %>%
  filter(tx_sd == 0.2) %>%
  group_by(nbar, J, ICC, tau, tx_sd, method, ATE) %>% 
  summarize(n = n(),
            power = mean(reject_two)) %>%
  filter(n >= 50, power >= 0.8) %>%
  summarize(MDES = min(ATE)) %>%   # minimum ATE s.t. power >= 0.8
  ggplot(aes(x=nbar, y=MDES, color=method)) +
  geom_point() +
  geom_line()

### average MDES

# not defined


### average power

# in theory, maximum power should be 1-qnorm(0, 0.2, 0.2) = 0.84
tidy_results %>% 
  filter(tx_sd == 0.2) %>%
  # filter(method == "MLM") %>% 
  group_by(nbar, J, ICC, tau, tx_sd, method) %>% 
  summarize(power = mean(reject_two)) %>% 
  ggplot(aes(x=nbar, y=power, color=method)) +
  geom_point(alpha=0.5) +
  geom_line(aes(group=method), size=1) +
  scale_color_manual(values = c(pal[5], pal[1])) +
  coord_cartesian(ylim = c(0,1)) +
  theme_minimal() +
  theme(axis.line = element_line()) +
  labs(y = "Average power",
       x = "Average site size",
       color = "Method",
       pch = TeX("$\\sigma^2$ value"))

# "true" power, only for tau_j != 0
tidy_results %>% 
  filter(tx_sd == 0.2, ATE != 0) %>%
  # filter(method == "MLM") %>% 
  group_by(nbar, J, ICC, tau, tx_sd, method) %>% 
  summarize(power = mean(reject_two)) %>% 
  ggplot(aes(x=nbar, y=power, color=method)) +
  geom_point(alpha=0.5) +
  geom_line(aes(group=method), size=1) +
  scale_color_manual(values = c(pal[5], pal[1])) +
  coord_cartesian(ylim = c(0,1)) +
  theme_minimal() +
  theme(axis.line = element_line()) +
  labs(y = "Average power",
       x = "Average site size",
       color = "Method",
       pch = TeX("$\\sigma^2$ value"))






# illustrative example ----------------------------------------------------

source("simulation_functions.R")

set.seed(123)   # use this with J=100
if (T) {
  QUANTILES <- c(0.05, 0.1, 0.95)
  NAME_VEC <- c("ATEhat", "SE", paste0("q", QUANTILES))
  J <- 100
  
  sdat = blkvar::generate_multilevel_data( n.bar = 50, 
                                           J = J,
                                           variable.n = T,
                                           tau.11.star = 0.2^2,   # cross-site tx var
                                           gamma.10 = 0.2,      # cross-site ATE
                                           ICC = 0.2,           # ICC [really var(intercepts)]
                                           rho2.0W = 0,         # covariate has no explanatory power 
                                           rho2.1W = 0,
                                           zero.corr = T,       # don't correlate site intercepts and treatment effects (so treatment group is higher variance)
                                           return.sites = TRUE,
                                           verbose = FALSE) %>%
    mutate(sid = as.character(1:n()))
  dat = blkvar::generate_individual_data(sdat,
                                         sigma2.e = 1-0.2)
  
  # single-site estimates
  res_single <- dat %>%
    group_by(sid) %>%
    dplyr::group_modify(~run_t_test(., NAME_VEC)) %>%
    ungroup() %>%
    mutate(sid = as.character(sid))
  
  # make dataset for bayesian models
  stan_df <- make_site_summaries(dat)
  stan_list <- list(
    J = J,
    site_mn_obs = stan_df$tau.hat,
    site_sd_obs = stan_df$SE)
  # run stan
  mod_norm <- stan_model("Stan/dp_normal_reparam.stan",
                         model_name = "norm",
                         auto_write = T)
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
  res_bayesnorm <- tibble(
    sid = as.character(1:J),
    ATEhat_bayesnorm = apply(site_effects_norm, 2, mean),
    SE_bayesnorm = apply(site_effects_norm, 2, sd)) %>%
    cbind(t(apply(site_effects_norm, 2, 
                  function(x) quantile(x, probs=QUANTILES)))) %>%
    as_tibble(.name_repair = "minimal")
  
  names(res_bayesnorm) <- c("sid", paste0(NAME_VEC, "_bayesnorm"))
}



# plot --------------------------------------------------------------------

spdf <- res_single %>% 
  mutate(method = "Single",
         tau = sdat$beta.1,
         tau_hat = ATEhat_single,
         q5 = q0.05_single,
         q95 = q0.95_single) %>% 
  select(sid, method, tau, tau_hat, q5, q95) %>% 
  bind_rows(
    res_bayesnorm %>% 
      mutate(method = "MLM",
             tau = sdat$beta.1,
             tau_hat = ATEhat_bayesnorm,
             q5 = q0.05_bayesnorm,
             q95 = q0.95_bayesnorm) %>% 
      select(sid, method, tau, tau_hat, q5, q95)
  )

# coverage percentages:
spdf %>% 
  group_by(method) %>% 
  summarize(coverage = mean(tau <= q95 & tau >= q5))

# tau_hat vs. tau
YMIN1 <- 0.55
YMAX1 <- 0.65

spdf %>% 
  ggplot(aes(x=tau, y=tau_hat)) +
  geom_linerange(aes(ymin=q5, ymax=q95), alpha=0.3) +
  geom_point(alpha=0.5) +
  geom_abline(lty = "dotted") +
  facet_wrap(~method, scales="free_y") +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines")) +
  coord_cartesian(ylim = c(-1.5, 2.5)) +
  labs(y = TeX("Estimated $\\hat{\\tau}_j$"),
       x = TeX("True $\\tau_j$"))
ggsave("writeup/images/shrinkageplot.png")

spdf %>% 
  mutate(in_range = tau >= YMIN1 & tau <= YMAX1) %>%
  ggplot(aes(x=tau, y=tau_hat)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=YMIN1, xmax=YMAX1),
            fill = "gray",
            alpha=0.5) +
  geom_linerange(aes(ymin=q5, ymax=q95,
                     alpha = ifelse(in_range, 1, 0.01))) +
  geom_point(aes(alpha = ifelse(in_range, 1, 0.01))) +
  geom_abline(lty = "dotted") +
  facet_wrap(~method, scales="free_y") +
  
  guides(alpha = "none",
         color = "none") +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines")) +
  coord_cartesian(ylim = c(-1.5, 2.5)) +
  labs(y = TeX("Estimated $\\hat{\\tau}_j$"),
       x = TeX("True $\\tau_j$"))
ggsave("writeup/images/shrinkageplot_slice1.png")

# YMIN2 <- 1
# YMAX2 <- 1.5
# spdf %>% 
#   mutate(in_range = tau_hat >= YMIN2 & tau_hat <= YMAX2) %>%
#   ggplot(aes(x=tau, y=tau_hat)) +
#   geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=YMIN2, ymax=YMAX2),
#             fill = pal[1],
#             alpha=0.01) +
#   geom_linerange(aes(ymin=q5, ymax=q95,
#                      alpha = ifelse(in_range, 1, 0.01))) +
#   geom_point(aes(alpha = ifelse(in_range, 1, 0.01))) +
#   geom_abline(lty = "dotted") +
#   facet_wrap(~method, scales="free_y") +
#   
#   guides(alpha = "none",
#          color = "none") +
#   theme_minimal() +
#   theme(axis.line = element_line(),
#         panel.spacing = unit(2, "lines")) +
#   coord_cartesian(ylim = c(-1.5, 2.5)) +
#   labs(y = TeX("Estimated $\\hat{\\tau}_j$"),
#        x = TeX("True $\\tau_j$"))
# ggsave("results_sree/shrinkageplot_slice2.png")



YMIN2 <- 0.4
YMAX2 <- 0.7
spdf %>% 
  mutate(in_range = tau_hat >= YMIN2 & tau_hat <= YMAX2) %>%
  ggplot(aes(x=tau, y=tau_hat)) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=YMIN2, ymax=YMAX2),
            fill = "gray",
            alpha=0.5) +
  geom_linerange(aes(ymin=q5, ymax=q95,
                     alpha = ifelse(in_range, 1, 0.01))) +
  geom_point(aes(alpha = ifelse(in_range, 1, 0.01))) +
  geom_abline(lty = "dotted") +
  facet_wrap(~method, scales="free_y") +
  
  guides(alpha = "none",
         color = "none") +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines")) +
  coord_cartesian(ylim = c(-1.5, 2.5)) +
  labs(y = TeX("Estimated $\\hat{\\tau}_j$"),
       x = TeX("True $\\tau_j$"))
ggsave("writeup/images/shrinkageplot_slice3.png")



