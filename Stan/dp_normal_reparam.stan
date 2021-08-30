data{
    int<lower=0> J;                     // number of sites
    vector[J] site_mn_obs;              // observed site-specific mean
    vector<lower=0>[J] site_sd_obs;     // observed site-specific sd
}
    
parameters{
    real pop_mn;            // grand mean
    real<lower=0> pop_sd;   // grand sd
    
    vector[J] eta;
}

transformed parameters{
    vector[J] site_mn;
    site_mn = pop_mn + pop_sd * eta;    // so site_mn ~ N(pop_mn, pop_sd^2)
}

model{
    // some weakly informative priors
    pop_mn ~ normal(0, 5);
    pop_sd ~ normal(0, 1);
    
    eta ~ normal(0,1);
    site_mn_obs ~ normal(site_mn, site_sd_obs);
}

generated quantities{
    // posterior predictive draw of true site impact
    real eta_new = normal_rng(0, 1);
    real y_site_pred = pop_mn + pop_sd * eta_new;
}
