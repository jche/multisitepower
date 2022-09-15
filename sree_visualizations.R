
require(tidyverse)
require(ggthemes)
require(wesanderson)
require(latex2exp)

pal <- wes_palette("Zissou1", 5, type="continuous")

tidy_results <- read_csv("results_sree/sree_sims.csv") %>% 
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= 0 | q95 <= 0,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05),
         method = case_when(
           method == "bayesnorm" ~ "MLM",
           method == "single" ~ "Single"))

MIN_SIZE <- 50




# power vs. sample size ---------------------------------------------------

# conditioned on tau_j
tidy_results %>% 
  filter(ATE %in% c(0, 0.2, 0.4)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
  summarize(n = n(),
            power = mean(reject_one)) %>% 
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
ggsave("results_sree/pp1.png")


# # visualize difference
# tidy_results %>% 
#   group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
#   summarize(n = n(),
#             power = mean(reject_one)) %>% 
#   filter(n >= MIN_SIZE) %>% 
#   pivot_wider(names_from=method, values_from=power) %>% 
#   mutate(diff = Bayes - Single) %>% 
#   ggplot(aes(x=nbar, y=diff, color = ATE)) +
#   geom_point(aes(size=n), alpha=0.1) +
#   # geom_line(aes(group=ATE)) +
#   geom_smooth(aes(group=ATE)) +
#   geom_hline(aes(yintercept=0)) +
#   scale_color_continuous(low="blue", high="orange") +
#   theme_minimal()


# conditioned on tau_hat_j
#  - note: hard to compare between Bayes and Single, since
#    range of tau_hat_j values is very different between them
tidy_results %>% 
  filter(ATEhat_bin %in% c(0, 0.2, 0.4)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATEhat_bin, method) %>% 
  summarize(n = n(),
            power = mean(reject_one)) %>% 
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
ggsave("results_sree/pp2.png")

# # visualize difference
# tidy_results %>% 
#   group_by(nbar, J, ICC, tau, tx_sd, ATEhat_bin, method) %>% 
#   summarize(power = mean(reject_one)) %>% 
#   # filter(ATEhat_bin >= -0.2, ATEhat_bin <= 0.6) %>%
#   pivot_wider(names_from=method, values_from=power) %>% 
#   mutate(diff = Bayes - Single) %>% 
#   
#   # ensure that only tau-hats that occur for many nbars
#   #  get geom_smoothed (to avoid errors)
#   filter(!is.na(diff)) %>% 
#   group_by(ATEhat_bin, J, ICC, tau, tx_sd) %>% 
#   filter(n() >= 25) %>% 
#   
#   ggplot(aes(x=nbar, y=diff, color = ATEhat_bin)) +
#   # geom_point(alpha=0.1) +
#   # geom_line(aes(group=ATEhat_bin)) +
#   geom_smooth(aes(group=ATEhat_bin)) +
#   geom_hline(aes(yintercept=0)) +
#   scale_color_continuous(low="blue", high="orange") +
#   theme_minimal()



# effect size vs sample size ----------------------------------------------

# sort of a weird plot: shows 
tidy_results %>% 
  group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
  summarize(n = n(),
            power = mean(reject_one)) %>% 
  filter(n >= MIN_SIZE) %>% 
  ggplot(aes(x=nbar, y=ATE)) +
  # geom_tile(aes(fill = power)) +
  # geom_contour(aes(z=power)) +
  geom_contour_filled(aes(z=power),
                      binwidth=0.1) +
  facet_wrap(~method) +
  # scale_fill_continuous(low="blue", high="orange") +
  theme_minimal()



# coverage vs sample size -------------------------------------------------

# conditioned on tau_j
tidy_results %>% 
  filter(ATE %in% c(0, 0.2, 0.4)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
  summarize(n = n(),
            coverage = mean(covered_one)) %>% 
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
  scale_y_continuous(limits=c(0, 1)) +
  labs(y = "Coverage",
       x = "Average site size",
       color = TeX("$\\tau_j$ value")) +
  guides(size="none")
ggsave("results_sree/cp1.png")

# conditioned on tau_hat_j
tidy_results %>% 
  filter(ATEhat_bin %in% c(-0.4, -0.1, 0.2, 0.5, 0.8)) %>%
  group_by(nbar, J, ICC, tau, tx_sd, ATEhat_bin, method) %>%
  summarize(n = n(),
            coverage = mean(covered_one)) %>% 
  filter(n >= MIN_SIZE) %>% 
  mutate(ATEhat_bin = factor(ATEhat_bin)) %>%
  ggplot(aes(x=nbar, y=coverage, color = ATEhat_bin)) +
  geom_point(aes(size=n), alpha=0.1) +
  # geom_line(aes(group=ATE)) +
  geom_smooth(aes(group=ATEhat_bin), se=F) +
  geom_hline(aes(yintercept=0.9), lty="dashed", color="black") +
  facet_wrap(~ method, scales="free") +
  # scale_color_continuous(low="blue", high="orange") +
  scale_color_manual(# values=pal[c(1,3,5)],
                     values = pal,
                     guide = guide_legend(reverse=T)) +
  theme_minimal() +
  theme(axis.line = element_line()) +
  labs(y = "Coverage",
       x = "Average site size",
       color = TeX("$\\hat{\\tau}_j$ value")) +
  guides(size="none")
ggsave("results_sree/cp2.png")






# conditioned on tau_j
tidy_results %>% 
  group_by(nbar, J, ICC, tau, tx_sd, ATE, method) %>% 
  summarize(n = n(),
            coverage = mean(covered_one)) %>% 
  filter(n >= MIN_SIZE) %>% 
  ggplot(aes(x=nbar, y=ATE)) +
  # geom_tile(aes(fill = coverage)) +
  # geom_contour(aes(z=coverage)) +
  geom_contour_filled(aes(z=coverage),
                      binwidth=0.1) +  
  facet_wrap(~method) +
  # scale_fill_continuous(low="blue", high="orange") +
  theme_minimal()

# conditioned on tau_hat_j
tidy_results %>% 
  group_by(nbar, J, ICC, tau, tx_sd, ATEhat_bin, method) %>% 
  summarize(n = n(),
            coverage = mean(covered_one)) %>% 
  filter(n >= MIN_SIZE) %>% 
  ggplot(aes(x=nbar, y=ATEhat_bin)) +
  # geom_tile(aes(fill = coverage)) +
  # geom_contour(aes(z=coverage)) +
  geom_contour_filled(aes(z=coverage),
                      binwidth=0.1) +  
  facet_wrap(~method) +
  # scale_fill_continuous(low="blue", high="orange") +
  theme_minimal()



# shrinkage example -------------------------------------------------------

source("simulation_functions.R")
require(rstan)

set.seed(123)   # for first example
set.seed(01720)   # for second example
set.seed(1)
if (T) {
  sdat = blkvar::generate_multilevel_data( n.bar = 50, 
                                           J = 25,
                                           variable.n = T,
                                           tau.11.star = 0.2^2,   # cross-site tx var
                                           gamma.10 = 0.2,      # cross-site ATE
                                           ICC = 0.2,           # ICC [really var(intercepts)]
                                           rho2.0W = 0,         # covariate has no explanatory power 
                                           rho2.1W = 0,
                                           zero.corr = T,       # don't correlate site intercepts and treatment effects (so treatment group is higher variance)
                                           return.sites = TRUE,
                                           verbose = FALSE) %>%
    mutate(sid = as.character(1:n()),
           beta.1 = round(beta.1/0.05) * 0.05)
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
    J = 25,
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
    sid = as.character(1:25),
    ATEhat_bayesnorm = apply(site_effects_norm, 2, mean),
    SE_bayesnorm = apply(site_effects_norm, 2, sd)) %>%
    cbind(t(apply(site_effects_norm, 2, 
                  function(x) quantile(x, probs=QUANTILES)))) %>%
    as_tibble(.name_repair = "minimal")
  
  QUANTILES <- c(0.05, 0.1, 0.95)
  NAME_VEC <- c("ATEhat", "SE", paste0("q", QUANTILES))
  names(res_bayesnorm) <- c("sid", paste0(NAME_VEC, "_bayesnorm"))
}


### plots

make.densodot.plot = function( X, group = NULL, binwidth, bw = binwidth, normal.density = FALSE ) {
  df = data.frame( X = X )
  if ( !is.null( group ) ) {
    df$group=group
  }
  
  # Hand-bin our dots
  scl = 1 / binwidth
  mn = round( scl * (min( X ) - binwidth/2 ) ) / scl
  breaks = seq( mn - binwidth/2, max( df$X + binwidth), by=binwidth )
  df = df %>% 
    # group_by(name) %>% 
    mutate(bin = cut( X, breaks=breaks ) )
  
  # mx <- df %>% 
  #   group_by(group, bin) %>% 
  #   count() %>% 
  #   pull(n) %>% 
  #   max()
  mx = max( table(df$bin ) )
  
  # Get density curve to plot
  if ( normal.density ) {
    dd = make.normal.density( df$X )
  } else {
    dens = density( df$X, bw=bw )
    dd = data.frame( x=dens$x, y=dens$y )
  }
  dmax = max( dd$y )
  
  # What fraction of density is in tallest histogram bar?
  frac = mx / nrow( df )
  
  # How high should density line be through the peak (to get relatively same area
  # under density curve (integrate curve over binwidth) vs. histogram bin (# dots in
  # the bin over total number of dots)
  ratio = (binwidth * dmax) / frac
  
  # Each unit of height is what in terms of dots? (The dots will stack up
  # without regard of y-axis, so we want to fix aspect ratio so the dots
  # correspond to the density line.)
  scaling = binwidth / ( (dmax / ratio) / (mx) )
  
  y.max = max( dmax, mx * binwidth/scaling )
  y.max = mx * binwidth/scaling
  
  
  if ( is.null( group ) ) {
    plt = ggplot( df )+
      geom_dotplot( aes(x=X), method="histodot",
                    binwidth = binwidth, stackgroups = TRUE)
  } else {
    plt = ggplot( df )+
      geom_dotplot( aes(x=X, fill=group, col=group), method="histodot",
                    binwidth = binwidth, stackgroups = TRUE)
  }
  
  print(y.max)
  plt = plt +
    # geom_line( data=dd, aes( x = x, y = y ) ) +
    # coord_fixed(ratio = scaling, ylim=c(0, y.max ) ) +
    coord_fixed(ratio = scaling, ylim=c(0, 1.6 ) ) +
    scale_y_continuous(name="", 
                       breaks=seq(0,
                                  by=binwidth/scaling, 
                                  length.out=(mx+1)), 
                       labels=c(0:mx) )
  
  plt
}
foo <- sdat %>% 
  select(sid, True = beta.1) %>% 
  bind_cols(Single = res_single$ATEhat_single,
            MLM = res_bayesnorm$ATEhat_bayesnorm) %>% 
  pivot_longer(True:MLM) %>% 
  mutate(name = factor(name, levels=c("Single", "MLM", "True")))

make.densodot.plot(foo$value, group=foo$name, binwidth=0.05) +
  facet_wrap(~foo$name, ncol=1) +
  guides(color="none", fill="none") +
  scale_color_manual(values=pal[c(1,3,5)]) +
  scale_fill_manual(values=pal[c(1,3,5)]) +
  labs(x="") +
  theme_minimal() +
  theme(axis.text.y=element_blank())
ggsave("results_sree/dotplot.png")



### intervals

# with 01720 seed
res_single %>% 
  mutate(method = "Single") %>% 
  select(sid, method, ATEhat = ATEhat_single, q0.1 = q0.1_single) %>% 
  bind_rows(
    res_bayesnorm %>%
      mutate(method = "MLM") %>% 
      select(sid, method, ATEhat = ATEhat_bayesnorm, q0.1 = q0.1_bayesnorm)) %>% 
  left_join(sdat %>% 
              select(sid, True = beta.1),
            by="sid") %>% 
  mutate(q0.1 = ifelse(q0.1 < -1, -1, q0.1),
         sid = fct_reorder(sid, -True)) %>% 
  ggplot(aes(x=as.numeric(sid))) +
  geom_point(aes(y=True)) +
  geom_errorbar(aes(ymin = q0.1, ymax=Inf)) +
  facet_wrap(~method, scales="free_y") +
  coord_flip() +
  labs(x = TeX("Site ID (sorted by $\\tau_j$)"),
       y = TeX("One-sided intervals for $\\tau_j$")) +
  theme_minimal() +
  theme(axis.line = element_line(),
        panel.spacing = unit(2, "lines"))
ggsave("results_sree/coverageplot.png")






# crazy shrinkage plot thing ----------------------------------------------

source("simulation_functions.R")
require(rstan)

set.seed(1)
set.seed(90210)
set.seed(01720)   # not bad


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


### plot!

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
ggsave("results_sree/shrinkageplot.png")

spdf %>% 
  mutate(in_range = tau >= YMIN1 & tau <= YMAX1) %>%
  ggplot(aes(x=tau, y=tau_hat)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=YMIN1, xmax=YMAX1),
            fill = "gray",
            alpha=0.01) +
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
ggsave("results_sree/shrinkageplot_slice1.png")

YMIN2 <- 1
YMAX2 <- 1.5
spdf %>% 
  mutate(in_range = tau_hat >= YMIN2 & tau_hat <= YMAX2) %>%
  ggplot(aes(x=tau, y=tau_hat)) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=YMIN2, ymax=YMAX2),
            fill = "gray",
            alpha=0.01) +
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
ggsave("results_sree/shrinkageplot_slice2.png")



YMIN2 <- 0.5
YMAX2 <- 1
spdf %>% 
  mutate(in_range = tau_hat >= YMIN2 & tau_hat <= YMAX2) %>%
  ggplot(aes(x=tau, y=tau_hat)) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=YMIN2, ymax=YMAX2),
            fill = "gray",
            alpha=0.01) +
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
ggsave("results_sree/shrinkageplot_slice3.png")


