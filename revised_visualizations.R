
#####
# Visualize simulation results
#####

require(tidyverse)
require(glue)
require(latex2exp)


#####
# Load simulation results
#####

if (F) {
  
  source("preprocess_cluster_results.R")
  
  # res <- read_csv("revised_results/revised_sim2-bf0803ddc413.csv")
  
  # hits <- res %>%
  #   mutate(nbar = as.factor(nbar),
  #          J = as.factor(J),
  #          ICC = as.factor(ICC),
  #          tau = as.factor(tau))
  
  hits <- res

  # Get ATEs per method (long data)
  ATEhats <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("ATEhat"), names_to = "method", values_to = "ATEhat") %>%
    dplyr::select(ATE, sid, method, ATEhat) %>%
    mutate(method = str_sub(method, 8)) %>%
    ungroup()
  
  # Get SEs per method
  SEs <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("SE_"), names_to = "method", values_to = "SE") %>%
    pull(SE)
  
  # get quantile estimates per method
  q5s <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("q0.05"), names_to = "method", values_to = "q5") %>%
    pull(q5)
  q10s <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("q0.1"), names_to = "method", values_to = "q10") %>%
    pull(q10)
  q95s <- hits %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID) %>%
    pivot_longer(contains("q0.95"), names_to = "method", values_to = "q95") %>%
    pull(q95)
  
  # join data
  tidy_results <- ATEhats %>%
    mutate(SE = SEs,
           q5 = q5s,
           q10 = q10s,
           q95 = q95s)
  
  # write_csv(tidy_results, "tidy_results.csv")
  
  # prep "overall" estimator results
  tidy_results_overall <- res_overall %>%
    mutate(method = case_when(method == "FIRC" ~ "firc_overall",
                              method == "RIRC" ~ "rirc_overall",
                              method == "Bayes" ~ "bayesnorm_overall",
                              method == "Single" ~ "single_overall"))
  # write_csv(tidy_results_overall, "tidy_results_overall.csv")
  
  # act as if we use overall estimate for each site
  ATE_lists <- tidy_results %>%
    group_by(nbar, J, ICC, tau, tx_sd, runID, method) %>%
    summarize(ATEs = list(ATE),
              sids = list(sid)) %>%
    summarize(ATEs = list(first(ATEs)),
              sids = list(first(sids)))
  tidy_results_overall <- tidy_results_overall %>%
    left_join(ATE_lists, by=c("nbar", "J", "ICC", "tau", "tx_sd", "runID")) %>%
    unnest(cols = c(ATEs, sids)) %>%
    rename(ATE = ATEs,
           sid = sids) %>%
    relocate(sid, .before = method) %>%
    relocate(ATE, .before = sid) %>%
    mutate(q5 = qnorm(0.05, mean=ATEhat, sd=SE),
           q10 = qnorm(0.10, mean=ATEhat, sd=SE),
           q95 = qnorm(0.95, mean=ATEhat, sd=SE))
  
  # get intervals and such
  # tidy_results <- tidy_results %>%
  #   mutate(reject_one = q10 >= 0,
  #          covered_one = q10 <= ATE,
  #          covered_two = (q5 <= ATE) & (ATE <= q95))
  # tidy_results_overall <- tidy_results_overall %>%
  #   mutate(reject_one = q10 >= 0,
  #          covered_one = q10 <= ATE,
    #        covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
    # mutate(sid = NA) %>%
    # relocate(sid, .before = method) %>%
    # relocate(ATE, .before = sid)
  
  final_df <- rbind(tidy_results, tidy_results_overall)
  write_csv(final_df, "tidy_results.csv")
  
} else {
  tidy_results <- read_csv("tidy_results.csv") %>%
    mutate(reject_one = q10 >= 0,
           covered_one = q10 <= ATE,
           covered_two = (q5 <= ATE) & (ATE <= q95))
}


#####
# set fixed simulation factors
#####

# 5 potential covariate factors to visualize:
#  * ATE (true)
#  * method
#  * nbar
#  * tx_sd
#  * tau
#  * J

INCLUDE_OVERALL <- F   # include overall estimators?
TAU_FIXED <- 0.2
J_FIXED <- 100
SIGMA_FIXED <- 0.2

# 5 potential outcomes to visualize:
#  * ATEhat
#  * SE
#  * MDES (one-sided, two-sided)
#  * power (one-sided, two-sided)
#  * coverage (one-sided, two-sided)

# plot a single trial
tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(nbar == 50, tx_sd == 0.3, tau == TAU_FIXED, J == 25, runID == 7) %>%
  ggplot(aes(x = ATE, color = method)) +
  geom_point(aes(y = ATEhat)) +
  geom_errorbar(aes(ymin = q5, 
                    ymax = q95)) +
  facet_wrap(~method) +
  geom_abline(slope = 1)


### prep summary dfs

MIN_SIZE <- 50

# prep power df
power_df <- tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            power = mean(reject_one)) %>%
  filter(n >= MIN_SIZE)

# prep coverage dfs
coverage2_df <- tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            coverage_two = mean(covered_two)) %>%
  filter(n >= MIN_SIZE)
coverage1_df <- tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            coverage_one = mean(covered_one)) %>%
  filter(n >= MIN_SIZE)


#####
# checking coverage
#####

# plot two-sided interval coverage vs. ATE size
# I UPDATED THIS
coverage2_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  filter(method %in% c("firc2", "rirc2", "bayesnorm", "single")) %>%
  mutate(method = case_when(
    method == "firc2" ~ "FIRC",
    method == "rirc2" ~ "RIRC",
    method == "bayesnorm" ~ "Bayes",
    method == "single" ~ "Single")) %>%
  ggplot() +
  geom_line(aes(x=ATE, y=coverage_two, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Conditional coverage rates (two-sided CI)",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Coverage",
       x = "True site ATE",
       color = "Method")
ggsave(glue("writeup/images/coverage_plot.png"), width=200, height=125, units="mm")

# plot one-sided interval coverage vs. ATE size
coverage1_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  ggplot() +
  geom_line(aes(x=ATE, y=coverage_one, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Conditional coverage rates (one-sided CI)",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Coverage",
       x = "True site ATE",
       color = "Method")
# ggsave(glue("writeup/images/coverage_plot.png"), width=200, height=125, units="mm")

# plot two-sided EB coverage
tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  group_by(nbar, J, ICC, tau, tx_sd, method) %>%
  summarize(coverage_two = mean(covered_two)) %>%
  ggplot(aes(x=nbar, y=coverage_two, color=method)) +
  geom_point() +
  geom_line(aes(group=method)) +
  facet_grid( ~ tx_sd, labeller=label_both, scales="free") +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Empirical Bayes coverage rates (two-sided CI)",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "EB Coverage",
       x = "Average site size",
       color = "Method")
# ggsave(glue("writeup/images/coverage_plot_eb.png"), width=200, height=75, units="mm")

# plot one-sided EB coverage
tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  group_by(nbar, J, ICC, tau, tx_sd, method) %>%
  summarize(coverage_one = mean(covered_one)) %>%
  ggplot(aes(x=nbar, y=coverage_one, color=method)) +
  geom_point() +
  geom_line(aes(group=method)) +
  facet_grid( ~ tx_sd, labeller=label_both, scales="free") +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Empirical Bayes coverage rates (one-sided CI)",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "EB Coverage",
       x = "Average site size",
       color = "Method")
# ggsave(glue("writeup/images/coverage_plot_eb.png"), width=200, height=75, units="mm")


#####
# visualizing MDES
#####

# just get simulation settings
settings_df <- tidy_results %>%
  group_by(nbar, J, tau, tx_sd, method) %>%
  summarize() %>%
  ungroup()

mdes_fun <- function(pow) {
  # power_df is grouped already
  mdes_df <- power_df %>%
    filter(n >= 10,
           power >= pow) %>%
    summarize(MDES = min(ATE))   # per simulation setting, grab min ATE that satisfies power
    
  
  # NA if the sim setting never gets to power=pow
  settings_df %>%
    left_join(mdes_df, by=names(settings_df)) %>%
    mutate(power_level = pow) %>%
    mutate(MDES = replace_na(MDES, 1))   # if we can't find any size, say MDES = 1 for clearer(?) visuals
}

# for a fixed power level, plot MDES
mdes_fun(0.8) %>%
  filter(J == J_FIXED) %>%
  mutate(nbar = as.factor(nbar)) %>%
  ggplot(aes(x=nbar, y=MDES, color=method)) +
  geom_point(position = position_dodge(width=0.5),
             size = 3) +
  geom_hline(aes(yintercept = 0.2), lty="dashed") +
  facet_grid(tau ~ tx_sd, labeller=label_both ) +
  labs(title = "Minimum detectable effect size (power = 0.8)",
       subtitle = glue("J = {J_FIXED}"),
       x = "Average site size",
       color = "Method")

# plot MDES against power threshold
map_df(seq(0.5, 0.9, by=0.05), mdes_fun) %>%
  bind_rows() %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  ggplot(aes(y=MDES, x=power_level, color=method)) +
  geom_point() +
  geom_line(aes(group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  labs(title = "Minimum detectable effect size",
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       x = "Power level",
       color = "Method")



#####
# checking power
#####

# Plot power vs. ATE size
# I UPDATED THIS
power_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  filter(method %in% c("firc2", "rirc2", "bayesnorm", "single")) %>%
  mutate(method = case_when(
    method == "firc2" ~ "FIRC",
    method == "rirc2" ~ "RIRC",
    method == "bayesnorm" ~ "Bayes",
    method == "single" ~ "Single")) %>%
  
  ggplot(aes(x = ATE, y = power, col = method, group=method ) ) +
  facet_grid(tx_sd ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = 0.1, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 0.6)) +
  labs(title = glue("Power (\u03B1 = 0.1) vs. true site effect"),
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/power_plot_J100.png"), width=200, height=125, units="mm")

# Plot power gaps
power_df_single <- power_df %>%
  filter(method == "single") %>%
  rename(power_single = power) %>%
  ungroup() %>%
  dplyr::select(-method)

power_df %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  left_join(power_df_single, by=c("nbar", "J", "tau", "tx_sd", "ATE")) %>%
  mutate(power_gap = power - power_single) %>%
  ggplot(aes(x = ATE, y = power_gap, col = method, group=method ) ) +
  facet_grid(tx_sd ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  geom_hline( yintercept = 0) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  labs(title = glue("Power (\u03B1 = 0.1) vs. true site effect"),
       subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method")
# ggsave(glue("writeup/images/power_plot_J300_diff.png"), width=200, height=125, units="mm")


### FOR POWER ANALYSIS EXAMPLE

# ALPHA <- 0.1
# tidy_results %>%
#   mutate(nbar = as.numeric(as.character(nbar))) %>%
#   filter(round(ATE, 2) %in% c(0.2, 0.3, 0.4)) %>%
#   mutate(reject = pvalue_one < ALPHA,
#          ATE = as.factor(ATE)) %>%
#   group_by(nbar, tx_sd, tau, J, method, ATE) %>%
#   summarize(power = mean(reject)) %>%
#   ggplot(aes(x = nbar, y = power, col = method, group=method ) ) +
#   facet_grid(~ ATE, labeller=label_both ) +
#   geom_line(alpha = 0.7) +
#   geom_hline( yintercept = 0.8, lty = "dashed" ) +
#   labs(title = glue("Power (\u03B1 = {ALPHA}) vs. site size"),
#        subtitle = "Overall average effect: 0.2",
#        y = "Power",
#        x = "Site size",
#        color = "Method")
# ggsave(glue("writeup/images/power_plot_ex.png"), width=200, height=75, units="mm")


#####
# visualizing power for sites with ATE=0.2
#####

ATE_of_interest <- 0.2

# for sites with ATE=0.2, plot power (across sim settings)
circ <- power_df %>% 
  filter( abs(ATE-ATE_of_interest) < 0.001 )

circ %>%
  filter(tau == TAU_FIXED, J == J_FIXED) %>%
  ggplot( aes( x=nbar, y=power, col=method, group=method ) ) +
  geom_point() +
  geom_line(aes(group=method)) +
  geom_hline(yintercept=0.9, lty="dashed") +
  facet_grid(~ tx_sd, labeller=label_both) +
  labs( title = glue("Power (\u03B1 = 0.1) to detect a site with ATE={ATE_of_interest}"), 
        x = "nbar" )
# ggsave("writeup/images/power_plot_ATE03.png", width=200, height=75, units="mm")

# for sites with ATE=0.2, plot estimated ATEs (across sim settings)
#  - note: # of sites with ATE=0.2 depends on J and tau
tidy_results %>% {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(method != "single") %>%
  filter(abs(ATE-ATE_of_interest) < 0.01) %>%
  filter(tau == TAU_FIXED, J == J_FIXED) %>%
  ggplot( aes( x=ATEhat ) ) +
  geom_density(aes(color=method, group=method), 
               position="identity", fill="transparent") +
  geom_vline( xintercept = ATE_of_interest, col="red" ) +
  facet_grid( tx_sd ~ nbar, labeller = label_both ) +
  labs(title = glue("Estimated ATEs for sites with ATE={ATE_of_interest}"))
# ggsave("writeup/images/power_plot_ATE03_dens.png", width=200, height=J_FIXED, units="mm")


#####
# RMSE plots
#####

# boxplots of RMSE values
tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  group_by(nbar, tx_sd, tau, J, method, runID) %>%
  summarize(rmse = sqrt(mean((ATE - ATEhat)^2))) %>%
  ggplot(aes(x=method, y=rmse, color=method)) +
  geom_boxplot() +
  facet_grid(tx_sd ~ nbar, labeller = label_both, scales="free_y") +
  labs(subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
       y = "RMSE",
       x = "Method") +
  guides(color = "none")
# ggsave("writeup/images/rmse_plot.png", width=200, height=150, units="mm")



#####
# check out standard error (SE) reductions
#####

# always better than "single", as expected.
tidy_results %>% 
  {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
  filter(tau == TAU_FIXED, J == J_FIXED) %>%
  group_by(nbar, J, ICC, tau, tx_sd, method) %>%
  summarize(SE = mean(SE)) %>%
  ggplot(aes(x = nbar, y = SE, color = method)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_grid( ~ tx_sd, labeller=label_both) +
  labs(subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"))




