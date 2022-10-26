
# analyzing case study results

# note: lots of divergences and ESS warnings...


require(tidyverse)
require(wesanderson)
require(latex2exp)

res <- read_csv("case_study/case_study_results.csv")
site_size_key <- tibble(
  sid = 1:10,
  n = c(551, 412, 343, 173, 464, 544, 499, 396, 197, 116)
)

res <- read_csv("case_study/case_study_states_results.csv")
site_size_key <- tibble(
  sid = 1:5,
  n = c(551, 928, 895, 1008, 309)
)



my_theme <- theme_minimal() +
  theme(axis.line = element_line())
pal <- wes_palette("Zissou1", 5, type="continuous")




# plotting gamma distributions --------------------------------------------

# mean: shape/rate, calibrate to 0.03
# skew: 2/sqrt(shape), allow this to vary

mn <- 0.03
beta <- 1
alpha <- mn*beta

tibble(
  mn = 0.03,
  beta = c(50,150,250)
) %>% 
  mutate(alpha = mn*beta) %>% 
  rowwise() %>% 
  mutate(samps = list(rgamma(500000, shape=alpha, rate=beta))) %>% 
  unnest(samps) %>% 
  mutate(beta = as.factor(beta)) %>% 
  ggplot(aes(x=samps, color=beta)) +
  geom_density(aes(group=beta), size=1) +
  scale_color_manual(values=pal[c(1,3,5)]) +
  # facet_grid(~beta, scales="free_x") +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8, 0.6),
        legend.background = element_rect(color = "black",
                                         fill = "white")) +
  labs(y = "",
       x = TeX("$\\tau_j$"),
       color = "     b")
ggsave("writeup/images/gamma.png")


# plotting raw intervals --------------------------------------------------

res %>% 
  filter(sig_tau == 0.01, rho == 0.6) %>% 
  sample_n(1000) %>% 
  ggplot(aes(x=tau_j, y=tau_j_hat)) +
  geom_pointrange(aes(ymin = tau_j_hat-1.96*se_j,
                      ymax = tau_j_hat+1.96*se_j)) +
  geom_abline() +
  facet_wrap(~method)



# plotting summaries ------------------------------------------------------

# summarize results
res_sum <- res %>% 
  mutate(reject = (q5 >= 0) | (q95 <= 0),
         covered = (q5 <= tau_j) & (q95 >= tau_j),
         moe = (q95-q5)/2) %>% 
  group_by(tau, sig_tau, alpha, sig_alpha, rho, sid, method) %>% 
  summarize(coverage = mean(covered),
            reject_rate = mean(reject),
            avg_moe = mean(moe)) %>% 
  left_join(site_size_key, by="sid")

### basic plots

# average moe
res_sum %>% 
  filter(rho == 0, sig_tau == 0.01) %>% 
  ggplot(aes(x=n, y=avg_moe, color=method)) +
  geom_point() +
  geom_line()

# average coverage
res_sum %>% 
  filter(rho == 0, sig_tau == 0.01) %>% 
  ggplot(aes(x=n, y=coverage, color=method)) +
  geom_point() +
  geom_hline(aes(yintercept=0.9), lty="dashed")

# power
res_sum %>% 
  filter(rho == 0, sig_tau == 0.01) %>% 
  ggplot(aes(x=n, y=reject_rate, color=method)) +
  geom_point() +
  geom_line()


### varying rho

# average moe
res_sum %>% 
  filter(sig_tau == 0.01) %>% 
  ggplot(aes(x=n, y=avg_moe, color=rho)) +
  geom_point() +
  geom_line(aes(group=rho)) +
  facet_wrap(~method)

# average coverage
res_sum %>% 
  filter(sig_tau == 0.01) %>% 
  ggplot(aes(x=n, y=coverage, color=rho)) +
  geom_point() +
  geom_hline(aes(yintercept=0.9), lty="dashed") +
  facet_wrap(~method)

# power
res_sum %>% 
  filter(sig_tau == 0.01) %>% 
  ggplot(aes(x=n, y=reject_rate, color=rho)) +
  geom_point() +
  geom_line(aes(group=rho)) +
  facet_wrap(~method)


### varying sig_tau

temp <- res_sum %>% 
  filter(rho == 0, method=="single", sig_tau==0.01)
# average moe
res_sum %>% 
  filter(rho == 0, method=="MLM") %>% 
  mutate(sig_tau = as.factor(sig_tau)) %>% 
  ggplot(aes(x=n, y=avg_moe)) +
  geom_point(aes(color=sig_tau)) +
  geom_line(aes(color=sig_tau, group=sig_tau)) +
  
  # add in baseline single-site stuff
  geom_point(data=temp) +
  geom_line(data=temp, lty="dashed") +
  
  scale_color_manual(values=pal) +
  labs(y = "Average margin of error",
       x = "Site size",
       color = TeX("$\\sigma_\\tau$")) +
  my_theme
ggsave("writeup/images/case_study_moe.png")

# average coverage
res_sum %>% 
  filter(rho == 0) %>% 
  ggplot(aes(x=n, y=coverage, color=sig_tau)) +
  geom_point() +
  # geom_line(aes(group=sig_tau)) +
  geom_hline(aes(yintercept=0.9), lty="dashed") +
  facet_wrap(~method)

# power
res_sum %>% 
  filter(rho == 0) %>% 
  ggplot(aes(x=n, y=reject_rate, color=sig_tau)) +
  geom_point() +
  geom_line(aes(group=sig_tau)) +
  facet_wrap(~method)



