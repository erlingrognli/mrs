data{
  int<lower=1> N, n_obs, n_beta, n_region;
  array [n_obs] int ind_id, ind_region;   // indices of subjects and brain regions
  array [n_region] int ind_measure; // index of which measure used for each region
  matrix [n_obs, n_beta] pred; // predictor variable matrix
  vector[n_obs] mri;
}
transformed data{
 vector[n_obs] log_mri = log(mri);
}
parameters{
  real<lower=0> icpt;
  vector[n_beta] beta;
  vector<lower=0>[3] sigma;
  vector[2] region_icpt_mu_raw;
  vector<lower=0>[3] region_icpt_sd; // measurewise hyperparameter on variance of region intercepts
  vector[n_region-3] region_icpt_raw;
  vector[N-1] id_icpt_raw;
  real<lower=0> id_icpt_sd;
}
transformed parameters{
  // fixing icpts at each level to log(1) for identification
  // regions 1-3 are selected thickness, area and volume regions
  // icpt is hence the thickness of the lateraloccipital cortex of participant 1
  // but is not an estimate of interest
  vector[3] measure_icpt = append_row(zeros_vector(1), measure_icpt_raw);
  vector[n_region] region_icpt = append_row(zeros_vector(3), region_icpt_raw);
  vector[N] id_icpt = append_row(zeros_vector(1), id_icpt_raw * id_icpt_sd);
}
model{
  vector[n_obs] mri_mod;
 
  // priors
  id_icpt_raw ~ std_normal(); // multiplier of .5 gives expected individual average deviations between 44 and 230 %
  id_icpt_sd ~ gamma(.1, .1);
  measure_icpt_raw[1] ~ normal(6.9, .5); // average area 1000 times thickness
  measure_icpt_raw[2] ~ normal(7.6, .5);  // average volume 2000 times thickness
  region_icpt_sd ~ normal(0, 1.6); // average multiplicative error for regions within 10
  region_icpt ~ normal(measure_icpt[ind_measure], region_icpt_sd[ind_measure]); // region intercepts centered on measure intercept, with hyperprior by measure
  
  icpt ~ normal(.7, .3); // weakly informative centered on a thickness of 2 mm
  beta ~ normal(0, 0.45); // assuming multiplicative effect of predictors from 0.48 - 2.1
  sigma ~ normal(0, 1.5); // multiplicative error assumed to be within 11 times the predictor

  // compute model predictions combining varying effects of 
  // brain area and individual, and regression coefficients
  
  mri_mod =  icpt + region_icpt[ind_region] + 
  
             id_icpt[ind_id] + 
             
             pred * beta;
             
  log_mri ~ normal(mri_mod, sigma[ind_measure[ind_region]]); // residuals allowed to vary over measure  
  
}
generated quantities{
  array [n_obs] real ppc;
  vector [n_obs] log_lik;
  vector [n_beta] beta_exp = exp(beta);
  
  {vector[n_obs] mri_mod;
  
  mri_mod =  icpt + measure_icpt[ind_measure[ind_region]] + region_icpt[ind_region] + 
  
             id_icpt[ind_id] + 
             
             pred * beta;
  
  ppc = normal_rng(mri_mod, sigma[ind_measure[ind_region]]);
  
  for(n in 1:n_obs){
    
    log_lik[n] = normal_lpdf(log_mri[n] | mri_mod[n], sigma[ind_measure[ind_region]][n]);}
  }
}
