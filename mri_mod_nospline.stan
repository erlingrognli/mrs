data{
  int<lower=1> N, n_obs, n_beta, n_region;
  array [n_obs] int ind_id, ind_region, ind_spl;   // indices of subjects, brain areas and splines
  matrix [n_obs, n_beta] ind_pred; // predictor variable matrix
  vector[n_obs] mri; // mri measurements and subject age
}
transformed data{
	// log transform the mri data
	vector[n_obs] log_mri = log(mri);
}
parameters{
  real<lower=1> nu;
  real<lower=0> alpha;
  vector[n_beta] beta;
  vector<lower=0>[3] sigma;
  vector<multiplier=2.5>[n_region - 1] region_icpt_raw;
  vector<multiplier=0.5>[N] id_icpt;
}
transformed parameters{
	// set the icpt of the last area to log(1)
  vector[n_region] region_icpt = append_row(region_icpt_raw, zeros_vector(1));
}
model{
  vector[n_obs] mri_mod;
 
  // priors
  nu ~ gamma(3,.3);
  id_icpt ~ std_normal();
  region_icpt_raw ~ std_normal();
  alpha ~ normal(0, 10);
  beta ~ normal(0, 0.3);
  sigma ~ student_t(4,0,2);

  // compute model predictions combingn varying effects of 
  // brain area and individual, and regression coefficients
  
  mri_mod = alpha +
  
            region_icpt[ind_region] + id_icpt[ind_id] + 
             
            ind_pred * beta;
  
  log_mri ~ student_t(nu, mri_mod, sigma[ind_spl]); // likelihood
}
generated quantities{
  array [n_obs] real ppc;
  vector[n_obs] log_lik;
  vector[n_beta] beta_exp = exp(beta);
  
  {vector[n_obs] mri_mod;
  
  mri_mod = alpha + 
   
             region_icpt[ind_region] + id_icpt[ind_id] + 
             
             ind_pred * beta;
  
  ppc = normal_rng(mri_mod, sigma[ind_spl]);
  
  for(n in 1:n_obs){
    
  log_lik[n] = normal_lpdf(log_mri[n] | mri_mod[n], sigma[ind_spl]);
  
  }}
}
