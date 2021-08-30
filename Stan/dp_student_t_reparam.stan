data{
    int<lower=0> J;                     // number of sites
    vector[J] y_site_obs;               // observed site-specific mean
    vector<lower=0>[J] sigma_site_obs;  // observed site-specific sd
    real<lower=3> nu;                   // df for random effects
}
    
parameters{
    real mu_true;                 // grand mean
    real<lower=0> sigma_true;     // grand sd

    vector[J] eta;
}

transformed parameters{
    vector[J] y_site_true;
    y_site_true = mu_true + sqrt((nu-2)/nu) * sigma_true * eta;
}


model{
    // some weakly informative priors
    mu_true ~ normal(0, 5);
    sigma_true ~ normal(0, 1);
    
    eta ~ student_t(nu, 0, 1);   // so y_site_true ~ t(nu, mu_true, nu/(nu-2)*sigma_true
    y_site_obs ~ normal(y_site_true, sigma_site_obs);
}

generated quantities{
    // posterior predictive draw of true site impact
    real eta_new = student_t_rng(nu, 0, 1);
    real y_site_new = mu_true + sqrt((nu-2)/nu) * sigma_true * eta_new;
    
    // posterior distribution of standard deviation statistic
    // (we could include all statistics here; just sd for now)
    real<lower=0> sdev = sigma_true;
}
