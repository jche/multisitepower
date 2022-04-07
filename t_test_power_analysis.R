
# power analysis for t-test

# conditional on tau

# conditional on tau-hat

# TAU <- c(0.1, 0.2, 0.5)
# NUMSAMP <- seq(10, 100, by=10)

sim_dat <- function(tau, n) {
  tibble(X = rnorm(n),
         Y = rnorm(n) + tau)
}

# power_fun <- function(tau, n, NUMSAMP=1000) {
#   map_dbl(1:NUMSAMP, function(x) {
#     dat <- sim_dat(tau, n)
#     t.test(dat)$p.value < 0.05
#   }) %>%
#     mean()
# }
# power_fun(0.1, 10)
# 
# res <- expand_grid(tau = TAU, n = NUMSAMP) %>%
#   rowwise() %>%
#   mutate(res = power_fun(tau, n))
# 
# # power curve!
# res %>%
#   mutate(tau = as.factor(tau)) %>%
#   ggplot(aes(x=n, y=res, color=tau, group=tau)) +
#   geom_line()


# TODO: condition on tau-hat?

# run t-test
# record tau-hat


# full sim ----------------------------------------------------------------


# returns:
#  - estimate
#  - whether true tau is covered
#  - whether null is rejected
power_fun_full <- function(tau, n, NUMSAMP=1000) {
  map_dfr(1:NUMSAMP, function(x) {
    dat <- sim_dat(tau, n)
    tt <- t.test(dat$X, dat$Y)
    
    est <- (tt$estimate[2] - tt$estimate[1]) %>% as.numeric()
    cover <- (tt$conf.int[1] < tau) & (tt$conf.int[2] > tau)
    reject <- tt$p.value < 0.05
    
    return(c(estimate = est, cover=cover, reject=reject))
  })
}
foo <- power_fun_full(0.1, 100)

# for each true value of tau, we have a bunch of tau-hats!

TAU <- seq(0.1, 0.5, by=0.05)
N <- 100

res <- expand_grid(tau = TAU, n = N) %>%
  rowwise() %>%
  mutate(res = list(power_fun_full(tau, n))) %>%
  unnest(res)

res_bin <- res %>%
  mutate(estimate = round(estimate, 2))

# vertical averages
vert_avg <- res_bin %>% 
  group_by(tau) %>%
  summarize(estimate = mean(estimate))

# horizontal averages(per ATEhat, average ATE)
hor_avg <- res_bin %>% 
  group_by(estimate) %>%
  summarize(tau = mean(tau))

# final plot!
res_bin %>%
  ggplot(aes(x=tau, y=estimate)) +
  geom_point() +
  geom_point(data=vert_avg, 
             color="green") +
  geom_point(data=hor_avg, 
             color="red") +
  geom_abline(aes(intercept=0, slope=1))


# coverage
res_bin %>%
  group_by(tau) %>%
  summarize(prop = mean(cover)) %>%
  ggplot(aes(x=tau, y=prop)) +
  geom_line()
res_bin %>%
  group_by(estimate) %>%
  summarize(prop = mean(cover)) %>%
  ggplot(aes(x=estimate, y=prop)) +
  geom_line()

# power
res_bin %>%
  group_by(tau) %>%
  summarize(prop = mean(reject)) %>%
  ggplot(aes(x=tau, y=prop)) +
  geom_line()
res_bin %>%
  group_by(estimate) %>%
  summarize(prop = mean(reject)) %>%
  ggplot(aes(x=estimate, y=prop)) +
  geom_line()



# averaging vertically: for each true tau value, what is my estimate?

# averaging horizontally: for each estimate, what is the mean of the distribution
#  of the tau values?
#  - in this case, we're doing a uniform distribution of tau values


# we're not attaching a distributional assumption to the tau values!
#  - we could if we wanted to!

# CONCLUSION: depending on the true distribution of tau values...
#  we can say that conditional on my point estimate being 0.2, 
#   I know that my average

