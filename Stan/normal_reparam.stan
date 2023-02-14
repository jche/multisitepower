data{
    int<lower=0> J;             // number of sites
    vector[J] tau_j_hat;        // observed site-specific mean
    vector<lower=0>[J] se_j;    // observed site-specific sd
}
    
parameters{
    real tau;                   // grand mean
    real<lower=0> sig_tau;      // grand sd
    
    vector[J] eta;
}

transformed parameters{
    vector[J] tau_j;
    tau_j = tau + sig_tau * eta;    // so tau_j ~ N(tau, sig_tau^2)
}

model{
    tau ~ normal(0, 0.1);           // prior: tau between +/- 0.3
    sig_tau ~ normal(0, 0.1);       // prior: sig_tau less than 0.3
    
    eta ~ normal(0,1);
    tau_j_hat ~ normal(tau_j, se_j);
}

generated quantities{
    // posterior predictive draw of true site impact
    real eta_new = normal_rng(0, 1);
    real y_site_pred = tau + sig_tau * eta_new;
}
