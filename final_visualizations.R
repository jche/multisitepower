
require(tidyverse)
require(glue)
require(latex2exp)


tidy_results <- read_csv("tidy_results.csv")
tidy_results <- tidy_results %>%
  mutate(reject_one = q10 >= 0,
         covered_one = q10 <= ATE,
         reject_two = q5 >= ATE | q95 <= ATE,
         covered_two = (q5 <= ATE) & (ATE <= q95)) %>%
  mutate(ATEhat_bin = plyr::round_any(ATEhat, 0.05))




TAU_FIXED <- 0.2
J_FIXED <- 100

# firc1 is using se.fixef + se.ranef, firc2 uses arm() samples
tidy_results_sm <- tidy_results %>%
  filter(J == J_FIXED, tau == TAU_FIXED) %>%
  filter(method %in% c("firc2", "rirc2", "bayesnorm", "single")) %>%
  mutate(method = case_when(
    method == "firc2" ~ "FIRC",
    method == "rirc2" ~ "RIRC",
    method == "bayesnorm" ~ "Bayes",
    method == "single" ~ "Single"))



# coverage/power, one-sided ------------------------------------------------

MIN_SIZE <- 50

# prep power df
power_df <- tidy_results_sm %>% 
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            power = mean(reject_one)) %>%
  filter(n >= MIN_SIZE)

# prep coverage dfs
coverage1_df <- tidy_results_sm %>% 
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            coverage_one = mean(covered_one)) %>%
  filter(n >= MIN_SIZE)

# plot one-sided interval coverage vs. ATE size
coverage1_df %>%
  ggplot() +
  geom_line(aes(x=ATE, y=coverage_one, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Conditional coverage rates (one-sided CI)",
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Coverage",
       x = "True site ATE",
       color = "Method") +
  theme_classic()
ggsave(glue("writeup/images/coverage_plot.png"), width=200, height=125, units="mm")

# Plot power vs. ATE size
# I UPDATED THIS
power_df %>%
  ggplot(aes(x = ATE, y = power, col = method, group=method ) ) +
  facet_grid(tx_sd ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  # geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = 0.1, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 0.6)) +
  labs(title = glue("Power (\u03B1 = 0.1) vs. true site effect (one-sided CI)"),
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method") +
  theme_classic()
ggsave(glue("writeup/images/power_plot_J100.png"), width=200, height=125, units="mm")



# coverage/power, two-sided -----------------------------------------------

MIN_SIZE <- 50

# prep power df
power_df2 <- tidy_results_sm %>% 
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            power = mean(reject_two)) %>%
  filter(n >= MIN_SIZE)

# prep coverage dfs
coverage_df2 <- tidy_results_sm %>% 
  group_by(nbar, tx_sd, tau, J, method, ATE) %>%
  summarize(n = n(),
            coverage_two = mean(covered_two)) %>%
  filter(n >= MIN_SIZE)

# plot one-sided interval coverage vs. ATE size
coverage_df2 %>%
  ggplot() +
  geom_line(aes(x=ATE, y=coverage_two, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Conditional coverage rates (two-sided CI)",
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Coverage",
       x = "True site ATE",
       color = "Method") +
  theme_classic()
# ggsave(glue("writeup/images/coverage_plot2.png"), width=200, height=125, units="mm")

# Plot power vs. ATE size
power_df2 %>%
  ggplot(aes(x = ATE, y = power, col = method, group=method ) ) +
  facet_grid(tx_sd ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  # geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = 0.1, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 0.6)) +
  labs(title = glue("Power (\u03B1 = 0.1) vs. true site effect (two-sided CI)"),
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method") +
  theme_classic()
# ggsave(glue("writeup/images/power_plot2.png"), width=200, height=125, units="mm")



# what if we use ATEhat_bin instead? --------------------------------------

MIN_SIZE <- 50

# prep power df
power_df2hat <- tidy_results_sm %>% 
  group_by(nbar, tx_sd, tau, J, method, ATEhat_bin) %>%
  summarize(n = n(),
            power = mean(reject_two)) %>%
  filter(n >= MIN_SIZE)

# prep coverage dfs
coverage_df2hat <- tidy_results_sm %>% 
  group_by(nbar, tx_sd, tau, J, method, ATEhat_bin) %>%
  summarize(n = n(),
            coverage_two = mean(covered_two)) %>%
  filter(n >= MIN_SIZE)

# plot two-sided interval coverage vs. ATE size
coverage_df2hat %>%
  ggplot() +
  geom_line(aes(x=ATEhat_bin, y=coverage_two, color=method, group=method)) +
  facet_grid(tx_sd ~ nbar, labeller=label_both) +
  geom_hline(yintercept=0.9, lty="dashed") +
  labs(title = "Conditional coverage rates (two-sided CI)",
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Coverage",
       x = "Estimated site ATE",
       color = "Method") +
  theme_classic()
# ggsave(glue("writeup/images/coverage_plot2.png"), width=200, height=125, units="mm")

# Plot power vs. ATE size
power_df2hat %>%
  ggplot(aes(x = ATEhat_bin, y = power, col = method, group=method ) ) +
  facet_grid(tx_sd ~ nbar, labeller=label_both ) +
  geom_line(alpha = 0.7) +
  # geom_hline( yintercept = 0.8, lty = "dashed" ) +
  geom_hline( yintercept = 0.1, lty = "dashed" ) +
  geom_vline( xintercept = 0 ) +
  coord_cartesian(xlim = c(-0.5, 0.6)) +
  labs(title = glue("Power (\u03B1 = 0.1) vs. true site effect (two-sided CI)"),
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Power",
       x = "Site-level ATE",
       color = "Method") +
  theme_classic()


# interval length, two-sided -----------------------------------------------

tidy_results_sm %>%
  mutate(length = q95 - q5) %>%
  group_by(nbar, J, ICC, tau, tx_sd, method) %>%
  summarize(mn_length = mean(length)) %>%
  ggplot(aes(x=nbar, y=mn_length, color=method, group=method)) +
  geom_point() +
  geom_line() +
  facet_grid(~tx_sd, labeller="label_both") +
  labs(title = "Average interval lengths (two-sided CI)",
       subtitle = glue("J = {J_FIXED}, \u03C4 = {TAU_FIXED}"),
       y = "Average interval length",
       x = "Average site size",
       color = "Method") +
  theme_classic()
ggsave(glue("writeup/images/length_plot.png"), width=200, height=75, units="mm")








# EB results...
# 
# 
# tidy_results_sm %>% 
#   {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
#   filter(J == J_FIXED, tau == TAU_FIXED) %>%
#   group_by(nbar, J, ICC, tau, tx_sd, method) %>%
#   summarize(coverage_one = mean(covered_one)) %>%
#   ggplot(aes(x=nbar, y=coverage_one, color=method)) +
#   geom_point() +
#   geom_line(aes(group=method)) +
#   facet_grid( ~ tx_sd, labeller=label_both, scales="free") +
#   geom_hline(yintercept=0.9, lty="dashed") +
#   labs(title = "Empirical Bayes coverage rates (one-sided CI)",
#        subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
#        y = "EB Coverage",
#        x = "Average site size",
#        color = "Method")
# 
# # EB power
# #  - NO: this is currently proportion of all samples that get rejected
# #  - we're only really interested in proportion of true rejections...
# tidy_results_sm %>% 
#   {if (!INCLUDE_OVERALL) filter(., !str_detect(method, "_overall")) else .} %>%
#   filter(J == J_FIXED, tau == TAU_FIXED) %>%
#   filter(ATE > 0) %>%
#   group_by(nbar, J, ICC, tau, tx_sd, method) %>%
#   summarize(# reject_one = mean(reject_one),
#             reject_one = mean(reject_one)) %>%
#   ggplot(aes(x=nbar, y=reject_one, color=method)) +
#   geom_point() +
#   geom_line(aes(group=method)) +
#   facet_grid( ~ tx_sd, labeller=label_both, scales="free") +
#   geom_hline(yintercept=0.9, lty="dashed") +
#   geom_hline(yintercept=0.1, lty="dashed") +
#   labs(title = "Empirical Bayes power (one-sided CI)",
#        subtitle = glue("J = {J_FIXED}, tau = {TAU_FIXED}"),
#        y = "EB Power",
#        x = "Average site size",
#        color = "Method")



