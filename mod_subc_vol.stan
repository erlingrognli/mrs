data{
  int<lower=1> N, n_obs, n_str;
  array [n_obs] int ind_id, ind_str;   // indices of subjects and brain structures
  vector[n_obs] mri;
  vector[n_obs] age;
  vector[n_obs] female;
  vector[N] icv;
  vector[n_obs] ocd;
}
transformed data{
 vector[n_obs] log_mri = log(mri);
}
parameters{
// hyperparameters
  real beta_icv;
  real<lower=0> sigma_alpha_id;
  real mu_beta_age;
  real<lower=0> sigma_beta_age;
  real mu_beta_female;
  real<lower=0> sigma_beta_female;
  real mu_beta_ocd;
  real<lower=0> sigma_beta_ocd;
// parameters
  real alpha;
  vector<lower=0>[n_str] sigma;
  vector[n_str] beta_age;
  vector[n_str] beta_female;
  vector[n_str] beta_ocd;
  vector[n_str-1] alpha_str_raw;
  vector[N] alpha_id;
}
transformed parameters{
  // identifying the intercepts per structure by fixing one to log(1)
  vector[n_str] alpha_str = append_row(zeros_vector(1), alpha_str_raw);
}
model{
  vector[n_obs] mri_pred;
  
// hyperpriors
  beta_icv ~ normal(0, .3);
  sigma_alpha_id ~ normal(0, .3);

  mu_beta_age ~ normal(0, .5);
  sigma_beta_age ~ normal(0, .3);

  mu_beta_female ~ normal(0, .5);
  sigma_beta_female ~ normal(0, .3);

  mu_beta_ocd ~ normal(0, .5);
  sigma_beta_ocd ~ normal(0, .3);

//priors
  alpha_str_raw ~ normal(0, .75); // structure intercepts defined as deviations from structure 1 for identification
  beta_age ~ normal(mu_beta_age, sigma_beta_age);
  beta_female ~ normal(mu_beta_female, sigma_beta_female);
  beta_ocd ~ normal(mu_beta_ocd, sigma_beta_ocd);
  
// individual intercept terms informed by data on intracranial volume
  alpha_id ~ normal(icv * beta_icv, sigma_alpha_id);
  alpha ~ normal(8, .1);
  sigma ~ normal(0, 1);
  
  mri_pred = alpha + alpha_str[ind_str] + alpha_id[ind_id] +
  
             age .* beta_age[ind_str] + female .* beta_female[ind_str] +
             
             ocd .* beta_ocd[ind_str];
  
  log_mri ~ normal(mri_pred, sigma[ind_str]);
  
}
generated quantities{
  array [n_obs] real ppc;
  vector [n_obs] log_lik;
  {vector[n_obs] mri_pred = alpha + alpha_str[ind_str] + alpha_id[ind_id] +
  
                            age .* beta_age[ind_str] + female .* beta_female[ind_str] +
             
                            ocd .* beta_ocd[ind_str];
  
  ppc = normal_rng(mri_pred, sigma[ind_str]);
  
  for(n in 1:n_obs){

    log_lik[n] = normal_lpdf(log_mri[n] | mri_pred[n], sigma[ind_str][n]);}
  }
}
