
require(tidyverse)
require(ggthemes)
require(wesanderson)
require(latex2exp)
require(rstan)

pal <- wes_palette("Zissou1", 5, type="continuous")

MIN_SIZE <- 50

# analyzing simulation results --------------------------------------------

# tidy_results <- read_csv("final_sims/example_full.csv") %>% 
tidy_results <- read_csv("final_sims/example2_full.csv") %>%
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
  group_by(nbar, J, ICC, tau, tx_sd, method, ATE) %>% 
  summarize(n = n(),
            power = mean(reject_two)) %>%
  filter(n >= 50, power >= 0.8) %>%
  summarize(MDES = min(abs(ATE))) %>%   # minimum ATE s.t. power >= 0.8
  ggplot(aes(x=nbar, y=MDES, color=method)) +
  geom_point() +
  geom_line()

### average MDES

# not defined


### average power

# in theory, maximum power should be 1-qnorm(0, 0.2, 0.2) = 0.84
tidy_results %>% 
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
  filter(ATE != 0) %>%
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




# SREE plots --------------------------------------------------------------

### power vs. sample size 

# conditioned on tau_j
tidy_results %>% 
  filter(ATE %in% c(0, 0.2, 0.4)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
  summarize(n = n(),
            power = mean(reject_two)) %>% 
  filter(n >= MIN_SIZE) %>% 
  mutate(ATE = as.factor(ATE)) %>%
  ggplot(aes(x=nbar, y=power, color = ATE)) +
  geom_point(aes(size=n), alpha=0.1) +
  geom_smooth(aes(group=ATE), se=F) +
  geom_hline(aes(yintercept=0.1), lty="dashed") +
  geom_hline(aes(yintercept=0.8), lty="dashed") +
  facet_wrap(~ method, scales="free") +   # scales=free is hack to make axes appear
  scale_color_manual(values=pal[c(1,3,5)],
                     guide = guide_legend(reverse=T)) +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines")) +
  scale_x_continuous(limits=c(0, 200)) + 
  scale_y_continuous(limits=c(0, 1)) +
  labs(y = "Power",
       x = "Average site size",
       color = TeX("$\\tau_j$ value")) +
  guides(size="none")
ggsave("writeup/images/pp1.png")


# conditioned on tau_hat_j
#  - note: hard to compare between Bayes and Single, since
#    range of tau_hat_j values is very different between them
tidy_results %>% 
  filter(ATEhat_bin %in% c(0, 0.2, 0.4)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATEhat_bin, method) %>% 
  summarize(n = n(),
            power = mean(reject_two)) %>% 
  filter(n >= MIN_SIZE) %>% 
  mutate(ATEhat_bin = as.factor(ATEhat_bin)) %>%
  # # ensure that only tau-hats that occur for many nbars
  # #  get geom_smoothed (to avoid errors)
  # group_by(J, ICC, tau, tx_sd, ATEhat_bin, method) %>% 
  # filter(n() >= 25) %>% 
  ggplot(aes(x=nbar, y=power, color = ATEhat_bin)) +
  geom_point(aes(size=n), alpha=0.2) +
  # geom_line(aes(group=ATE)) +
  geom_smooth(aes(group=ATEhat_bin), se=F) +
  facet_wrap(~ method, scales="free") +
  scale_color_manual(values=pal[c(1,3,5)],
                     guide = guide_legend(reverse=T)) +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines")) +
  scale_x_continuous(limits=c(0, 200)) + 
  scale_y_continuous(limits=c(0, 1)) +
  labs(y = "Power",
       x = "Average site size",
       color = TeX("$\\hat{\\tau}_j$ value")) +
  guides(size="none")
ggsave("writeup/images/pp2.png")



### coverage vs. sample size

# conditioned on tau_j
tidy_results %>% 
  filter(ATE %in% c(0, 0.2, 0.4)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
  summarize(n = n(),
            coverage = mean(covered_two)) %>% 
  filter(n >= MIN_SIZE) %>% 
  mutate(ATE = as.factor(ATE)) %>%
  ggplot(aes(x=nbar, y=coverage, color = ATE)) +
  geom_point(aes(size=n), alpha=0.1) +
  # geom_line(aes(group=ATE)) +
  geom_smooth(aes(group=ATE), se=F) +
  geom_hline(aes(yintercept=0.9), lty="dashed", color="black") +
  facet_wrap(~ method, scales="free") +
  # scale_color_continuous(low="blue", high="orange") +
  scale_color_manual(values=pal[c(1,3,5)],
                     guide = guide_legend(reverse=T)) +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines")) +
  scale_x_continuous(limits=c(0, 200)) + 
  scale_y_continuous(limits=c(0.5, 1)) +
  labs(y = "Coverage",
       x = "Average site size",
       color = TeX("$\\tau_j$ value")) +
  guides(size="none")
ggsave("writeup/images/cp1.png")

# conditioned on tau_hat_j
tidy_results %>% 
  filter(ATEhat_bin %in% c(-0.4, -0.1, 0.2, 0.5, 0.8)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATEhat_bin, method) %>%
  summarize(n = n(),
            coverage = mean(covered_two)) %>% 
  filter(n >= MIN_SIZE) %>% 
  mutate(ATEhat_bin = factor(ATEhat_bin)) %>%
  
  # ensure that enough points exist for geom_smooth() to work nicely
  group_by(method, ATEhat_bin) %>% 
  filter(n() >= 10) %>% 
  
  ggplot(aes(x=nbar, y=coverage, color = ATEhat_bin)) +
  geom_point(aes(size=n), alpha=0.1) +
  geom_smooth(aes(group=ATEhat_bin), se=F) +
  geom_hline(aes(yintercept=0.9), lty="dashed", color="black") +
  
  facet_wrap(~ method, scales="free") +
  scale_y_continuous(limits=c(0.5, 1)) +
  scale_color_manual(# values=pal[c(1,3,5)],
    values = pal,
    guide = guide_legend(reverse=T)) +
  theme_minimal() +
  theme(axis.line = element_line()) +
  labs(y = "Coverage",
       x = "Average site size",
       color = TeX("$\\hat{\\tau}_j$ value")) +
  guides(size="none")
ggsave("writeup/images/cp2.png")




# illustrative example ----------------------------------------------------

# run single simulation
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

# aggregate results
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

# vertical slice
YMIN1 <- 0.55
YMAX1 <- 0.65
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

# horizontal slice
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





# full simulation study ---------------------------------------------------

### wrapper functions
# IDEA: For J, tau, tx_sd, ICC as [x]...
#  - plot [outcome] vs. nbar, colored by [x] levels
#  - facet by conditioning variable
# Only keep plots with differences
if (T) {
  # function to plot faceted power curves (almost lol)
  power_plot <- function(d, color_var, cond_var, 
                         color_caption, 
                         cond_var_values=c(0, 0.2, 0.4)) {
    d %>% 
      filter({{cond_var}} %in% cond_var_values) %>%
      mutate(across(c({{color_var}}, {{cond_var}}),
                    ~as.factor(.))) %>%
      group_by(nbar, J, ICC, tau, tx_sd, {{cond_var}}, method) %>% 
      summarize(n = n(),
                power = mean(reject_two)) %>% 
      filter(n >= MIN_SIZE) %>% 
      
      ggplot(aes(x=nbar, y=power, color = {{color_var}})) +
      geom_point(aes(size=n), alpha=0.1) +
      geom_smooth(aes(group={{color_var}}), se=F) +
      geom_hline(aes(yintercept=0.1), lty="dashed") +
      geom_hline(aes(yintercept=0.8), lty="dashed") +    
      scale_color_manual(values=pal[c(1,3,5)],
                         guide = guide_legend(reverse=T)) +
      theme_minimal() +
      theme(axis.line = element_line(),
            panel.spacing = unit(2, "lines")) +
      scale_x_continuous(limits=c(0, 200)) + 
      scale_y_continuous(limits=c(0, 1)) +
      labs(y = "Power",
           x = "Average site size",
           color = color_caption) +
      guides(size="none")
  }
  
  # function to plot faceted coverage curves (almost lol)
  coverage_plot <- function(d, color_var, cond_var, 
                            color_caption, 
                            cond_var_values=c(0, 0.2, 0.4)) {
    d %>% 
      filter({{cond_var}} %in% cond_var_values) %>%
      mutate(across(c({{color_var}}, {{cond_var}}),
                    ~as.factor(.))) %>%
      group_by(nbar, J, ICC, tau, tx_sd, {{cond_var}}, method) %>% 
      summarize(n = n(),
                coverage = mean(covered_two)) %>% 
      filter(n >= MIN_SIZE) %>% 
      
      ggplot(aes(x=nbar, y=coverage, color = {{color_var}})) +
      geom_point(aes(size=n), alpha=0.1) +
      geom_smooth(aes(group={{color_var}}), se=F) +
      geom_hline(aes(yintercept=0.9), lty="dashed", color="black") +
      
      scale_color_manual(values=pal[c(1,3,5)],
                         guide = guide_legend(reverse=T)) +
      theme_minimal() +
      theme(axis.line = element_line(),
            panel.spacing = unit(2, "lines")) +
      scale_x_continuous(limits=c(0, 200)) + 
      scale_y_continuous(limits=c(0.5, 1)) +
      labs(y = "Coverage",
           x = "Average site size",
           color = color_caption) +
      guides(size="none")
  }
  
  length_plot <- function(d, color_var, color_caption) {
    baseline_length <- full_study_J %>% 
      filter(method == "Single") %>% 
      group_by(nbar) %>% 
      summarize(interval_length = mean(q95 - q5)/2)
    
    d %>% 
      filter(method == "MLM") %>%
      mutate(across({{color_var}}, ~as.factor(.))) %>% 
      group_by(nbar, J, ICC, tau, tx_sd, method) %>% 
      summarize(interval_length = mean(q95 - q5)/2) %>% 
      ggplot(aes(x=nbar, y=interval_length, color={{color_var}})) +
      geom_point(alpha=0.5) +
      geom_line(aes(group={{color_var}}),
                size=1) +
      geom_line(data  = baseline_length,
                color = "black",
                lty   = "dashed") +
      # facet_wrap(~method) +
      scale_color_manual(values = c(pal[c(1,3,5)])) +
      coord_cartesian(ylim = c(0,1)) +
      theme_minimal() +
      theme(axis.line = element_line()) +
      labs(y = "Average margin of error",
           x = "Average site size",
           color = color_caption)
  }
}

### varying J
full_study_J <- read_csv("final_sims/full_study_J.csv") %>%
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= 0 | q95 <= 0,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05),
         method = case_when(
           method == "bayesnorm" ~ "MLM",
           method == "single" ~ "Single"))
caption <- TeX("$J$ value")

# power conditioned on tau_j
power_plot(full_study_J, J, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# power conditioned on tau_hat_j
power_plot(full_study_J, J, ATEhat_bin, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# coverage conditioned on tau_j
coverage_plot(full_study_J, J, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# coverage conditioned on tau_hat_j
coverage_plot(full_study_J, J, ATEhat_bin, caption, c(-0.4, -0.1, 0.2, 0.5, 0.8)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# average interval length: slightly decreases with J
length_plot(full_study_J, J, caption)
ggsave("writeup/images/simstudy_J_length.png")



### varying tau
full_study_tau <- read_csv("final_sims/full_study_tau.csv") %>%
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= 0 | q95 <= 0,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05),
         method = case_when(
           method == "bayesnorm" ~ "MLM",
           method == "single" ~ "Single"))
caption <- TeX("$\\tau$ value")

# power conditioned on tau_j: obviously differs with tau
power_plot(full_study_tau, tau, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# power conditioned on tau_hat_j
power_plot(full_study_tau, tau, ATEhat_bin, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# coverage conditioned on tau_j
coverage_plot(full_study_tau, tau, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# coverage conditioned on tau_hat_j
coverage_plot(full_study_tau,tau, ATEhat_bin, caption, c(-0.4, -0.1, 0.2, 0.5, 0.8)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# average interval length
length_plot(full_study_tau, tau, caption)



### varying tx_sd

full_study_txsd <- read_csv("final_sims/full_study_txsd.csv") %>%
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= 0 | q95 <= 0,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05),
         method = case_when(
           method == "bayesnorm" ~ "MLM",
           method == "single" ~ "Single"))
caption <- TeX("$\\sigma_\\tau$")

# power conditioned on tau_j
power_plot(full_study_txsd, tx_sd, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# power conditioned on tau_hat_j
power_plot(full_study_txsd, tx_sd, ATEhat_bin, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# coverage conditioned on tau_j
coverage_plot(full_study_txsd, tx_sd, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# coverage conditioned on tau_hat_j
coverage_plot(full_study_txsd, tx_sd, ATEhat_bin, caption, c(-0.4, -0.1, 0.2, 0.5, 0.8)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# average interval length
length_plot(full_study_txsd, tx_sd, caption)
ggsave("writeup/images/simstudy_txsd_length.png")



### varying ICC

full_study_icc <- read_csv("final_sims/full_study_ICC.csv") %>%
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= 0 | q95 <= 0,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05),
         method = case_when(
           method == "bayesnorm" ~ "MLM",
           method == "single" ~ "Single"))
caption <- TeX("$ICC$ value")

# power conditioned on tau_j
power_plot(full_study_icc, ICC, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# power conditioned on tau_hat_j
power_plot(full_study_icc, ICC, ATEhat_bin, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# coverage conditioned on tau_j
coverage_plot(full_study_icc, ICC, ATE, caption, c(0, 0.2, 0.4)) +
  facet_grid(ATE ~ method, scales="free")
# coverage conditioned on tau_hat_j
coverage_plot(full_study_icc, ICC, ATEhat_bin, caption, c(-0.4, -0.1, 0.2, 0.5, 0.8)) +
  facet_grid(ATEhat_bin ~ method, scales="free")

# average interval length
length_plot(full_study_icc, ICC, caption)






