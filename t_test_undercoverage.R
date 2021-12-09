
# test coverage of a t-test

# CONCLUSION: t-test "undercovers" when I approximate t dist with normal dist
#  - since "stderr" output from a t test is just the denominator of the t statistic

# run a t-test on two normal samples of size n
# returns whether the t-test rejects H0: m1 = m2
sim_t_test <- function(n, ALPHA=0.1) {
  vec1 <- rnorm(n, mean=1)
  vec2 <- rnorm(n)
  
  t <- t.test(vec1, vec2, conf.level = 1-ALPHA)
  CI <- t$conf.int
  
  # CI
  # mean(CI) - qt(1-ALPHA/2, df=t$parameter) * t$stderr
  # mean(CI) + qt(1-ALPHA/2, df=t$parameter) * t$stderr
  
  lb <- CI[1]
  ub <- CI[2]
  
  # lb <- mean(CI) - qnorm(1-ALPHA/2) * t$stderr
  # ub <- mean(CI) + qnorm(1-ALPHA/2) * t$stderr
  
  # (CI[1] <= 1) & (CI[2] >= 1)
  (lb <= 1) & (ub >= 1)
  
  # SE <- t$stderr
  # pt_est <- t$estimate["mean of x"] - t$estimate["mean of y"]
  # (pt_est - qnorm(1-ALPHA/2)*SE <= 1) & (pt_est + qnorm(1-ALPHA/2)*SE >= 1)
}

# TODO: need to record df... or just directly record CI values...?

(1:10000) %>%
  map_dbl(~sim_t_test(n=10)) %>%
  mean()

# For a sequence of sample sizes, average rejection rates over 10000 samples
#  - note: takes about a minute to run
set.seed(90210)
samp_sizes <- seq(10, 50, by=10)
res <- samp_sizes %>%
  map_dbl(function(x) {
    (1:10000) %>%
      map_dbl(~sim_t_test(n=x)) %>%
      mean()
  })

# plot results
tibble(
  n = samp_sizes,
  rr = res
) %>%
  ggplot(aes(x=n, y=rr)) +
  geom_point() +
  geom_line()

