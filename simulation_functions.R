
#####
# simulation helper functions
#####

# given a df for a single site, run a one-sided t-test
run_t_test <- function(df) {
  
  df %>%
    summarize(t = list(t.test(Y1[Z==1], Y0[Z==0], alternative = "greater"))) %>%
    mutate(ATEhat_single = t[[1]]$estimate["mean of x"] - t[[1]]$estimate["mean of y"],
           SE_single = t[[1]]$stderr
           # t_single = t[[1]]$statistic,
           # pvalue_single = t[[1]]$p.value
    ) %>%
    select(-t)
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

