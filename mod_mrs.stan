data{
  int<lower=1> N, n_obs, n_str;
  array [n_obs] int ind_id, ind_str; // indices of subjects and brain structures
  vector[n_obs] mri;
  vector[n_obs] age;
  vector[n_obs] female;
  vector[N] icv;
  vector[n_obs] ocd; // one-hot encoded ocd diagnostic status
  vector[n_obs] eop; // one-hot encoded eop diagnostic status
  vector[2] alpha_params;
  int<lower=0, upper=1> predict;
}
transformed data{
 vector[n_obs] log_mri = log(mri);
}
parameters{
// hyperparameters
  real beta_icv;
  real<lower=0> sigma_alpha_id;
  real mu_beta_age;
  real mu_beta_female;
  real mu_beta_ocd;
  real mu_beta_eop;
// parameters
  real alpha;
  vector<lower=0>[n_str] sigma;
  vector[n_str] beta_age;
  vector[n_str] beta_female;
  vector[n_str] beta_ocd;
  vector[n_str-1] alpha_str_raw;
  vector[n_str] beta_eop;
  vector[N] alpha_id;
}
transformed parameters{
  // identifying the intercepts per structure by fixing one to log(1)
  vector[n_str] alpha_str = append_row(zeros_vector(1), alpha_str_raw);
}
model{
  vector[n_obs] mri_pred;
  
// hyperpriors
  sigma_alpha_id ~ normal(0, .25);
  // assuming that the average multiplicative deviation from the expectation for 
  // subject random effect, conditional on intracranial volume, is within 1.65
  
  beta_icv ~ normal(0, .4);
  mu_beta_age ~ normal(0, .4);
  mu_beta_female ~ normal(0, .4);
  mu_beta_ocd ~ normal(0, .4);
  mu_beta_eop ~ normal(0, .4);
  // hyperpriors encode an assumption that the average multiplicative effect of 
  // 0-1 scaled age, female gender, ocd and z-score icv across structures 
  // are between .44 and 2.22 with 95% certainty
 
// priors
  alpha_str_raw ~ normal(0, .8);
  // structure intercepts defined as multiplicative deviations from structure #1 
  // for model identifiability, and prior assumes structure intercepts to have
  // multiplicative deviations from #1 between .2 and 4.95
  beta_age ~ normal(mu_beta_age, .4);
  beta_female ~ normal(mu_beta_female, .4);
  beta_ocd ~ normal(mu_beta_ocd, .4);
  beta_eop ~ normal(mu_beta_eop, .4);
  // variance of .4 encodes a general assumption that multiplicative variability  
  // in age, gender and ocd effects across structures is no larger than 2.22
  
  alpha_id ~ normal(icv * beta_icv, sigma_alpha_id);
  // individual intercept terms informed by data on intracranial volume
  alpha ~ normal(alpha_params[1], alpha_params[2]);
  // overall intercept scales the rest of the model parameters, as they are multiplicative
  // and is supplied as data, to allow for modeling of thickness, volume and area
  sigma ~ normal(0, 1);
  
  mri_pred = alpha + alpha_str[ind_str] + alpha_id[ind_id] +
  
             age .* beta_age[ind_str] + female .* beta_female[ind_str] +
             
             eop .* beta_eop[ind_str] + 
             
             ocd .* beta_ocd[ind_str];
  
  log_mri ~ normal(mri_pred, sigma[ind_str]);
}
generated quantities{
  
  if(predict == 0){
  
    array [n_obs] real ppc;
    vector [n_obs] log_lik;
    vector [n_str] exp_beta_ocd = exp(beta_ocd);
    vector [n_str] exp_beta_eop = exp(beta_eop);
 
    {vector[n_obs] mri_pred = alpha + alpha_str[ind_str] + alpha_id[ind_id] +
  
                            age .* beta_age[ind_str] + female .* beta_female[ind_str] +
                            
                            eop .* beta_eop[ind_str] + 
             
                            ocd .* beta_ocd[ind_str];
  
    ppc = normal_rng(mri_pred, sigma[ind_str]);
  
    for(n in 1:n_obs){

      log_lik[n] = normal_lpdf(log_mri[n] | mri_pred[n], sigma[ind_str][n]);}
    }}
    
    else{
      
      array[3, n_str] real predictions;
      real mean_icv = mean(icv);
      real mean_age = mean(age);
      
      for(k in 1:n_str){ // gender data is encoded +/- .5
                         // random intercepts are based on icv, which is centered
                         // both are omitted from predictions
          
          predictions[1,k] = normal_rng(alpha + alpha_str[k] + mean_age * beta_age[k], sigma[k]);
          predictions[2,k] = normal_rng(alpha + alpha_str[k] + mean_age * beta_age[k] + beta_ocd[k], sigma[k]);
          predictions[3,k] = normal_rng(alpha + alpha_str[k] + mean_age * beta_age[k] + beta_eop[k], sigma[k]);
      }}
    
}
