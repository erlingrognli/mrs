data{
  int<lower=1> N, n_obs, n_str;
  array [n_obs] int ind_id, ind_str; // indices of subjects and brain structures
  vector[n_obs] mri;
  vector[n_obs] age;
  vector[n_obs] gender; // gender coded +/- .5
  vector[N] icv;
  vector[n_obs] ocd; // one-hot encoded ocd diagnostic status
  vector[n_obs] eos; // one-hot encoded eop diagnostic status
  vector[2] alpha_params;
}
transformed data{
 vector[n_obs] log_mri = log(mri);
 real sigma_alpha_str_multiplier = sqrt(n_str * inv(n_str - 1));
 // multiplier for the variance of a normally distributed zero-sum vector, 
 // as shown by Sean Pinkney, with reference to "Fraser, D. A. S. (1951). 
 // Normal Samples With Linear Constraints and Given Variances. 
 // Canadian Journal of Mathematics, 3, 363–366. doi:10.4153/CJM-1951-041-9".
}
parameters{
// hyperparameters
  real beta_icv;
  real<lower=0> sigma_alpha_id;
  real<lower=0> sigma_alpha_str;
  real<lower=0> sigma_beta_age;
  real<lower=0> sigma_beta_gender;
  real mu_beta_age;
  real mu_beta_gender;
// parameters
  real alpha;
  vector<lower=0>[n_str] sigma;
  vector[n_str] beta_age;
  vector[n_str] beta_gender;
  vector[N] alpha_id;
  sum_to_zero_vector[n_str] alpha_str; 
// intercepts for structures get a sum to zero constraint, for model identification
}                                 
model{
  vector[n_obs] mri_pred;
  
// hyperpriors
  sigma_alpha_id ~ normal(0, .25);
  // assuming that the average multiplicative deviation from the expectation for 
  // subject random effect, conditional on intracranial volume, is within 1.6
  sigma_alpha_str ~ normal(0, .8);
  // assuming that the average multiplicative deviation from the mean for the
  // intercepts for each structure is within 5
  sigma_beta_age ~ lognormal(-1, 1);
  sigma_beta_gender ~ lognormal(-1, 1);
  beta_icv ~ normal(0, .4);
  mu_beta_age ~ normal(0, .4);
  mu_beta_gender ~ normal(0, .4);
  // hyperpriors encode an assumption that the average multiplicative effect of 
  // 0-1 scaled age, gender difference, ocd and z-score icv across structures 
  // are between .44 and 2.22 with 95% certainty
 
// priors
  alpha_str ~ normal(0, sigma_alpha_str * sigma_alpha_str_multiplier);
  beta_age ~ normal(mu_beta_age, sigma_beta_age);
  beta_gender ~ normal(mu_beta_gender, sigma_beta_gender);
  // variance of .4 encodes a general assumption that multiplicative variability  
  // in age, gender and ocd effects across structures is no larger than 2.22
  
  alpha_id ~ normal(icv * beta_icv, sigma_alpha_id);
  // individual intercept terms informed by data on intracranial volume
  alpha ~ normal(alpha_params[1], alpha_params[2]);
  // overall intercept scales the rest of the model parameters, as they are multiplicative
  // and is supplied as data, to allow for modeling of thickness, volume and area
  // with same code
  sigma ~ normal(0, 1);
  
  mri_pred = alpha + alpha_str[ind_str] + alpha_id[ind_id] +
  
             age .* beta_age[ind_str] + gender .* beta_gender[ind_str];
  
  log_mri ~ normal(mri_pred, sigma[ind_str]);
}
generated quantities{
  vector [n_obs] log_lik;

  {vector[n_obs] mri_pred = alpha + alpha_str[ind_str] + alpha_id[ind_id] +
  
                            age .* beta_age[ind_str] + gender .* beta_gender[ind_str];
    
   for(n in 1:n_obs){
      
      log_lik[n] = normal_lpdf(log_mri[n] | mri_pred[n], sigma[ind_str][n]);}}
}
