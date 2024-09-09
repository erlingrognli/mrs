functions{
  // code for splines by Sergey Koposov, doi: 10.5281/zenodo.7193909
  
  // get the vector of spacings between nodes
  vector spline_geths(vector nodes)
  {
    int n = size(nodes)-1;
    vector[n] hs;
    for (i in 1:n)
      {
        hs[i] = nodes[i+1]-nodes[i];
      }
    return hs;
  }
  
  // obtain the vector of spline coefficients given the location
  // of the nodes and values there
  // We are using natural spline definition           
  vector spline_getcoeffs(vector nodes, vector vals)
  {
    int n_nodes = size(nodes);
    int n=n_nodes-1;
    vector[n] hi;
    vector[n] bi;
    vector[n-1] vi;
    vector[n-1] ui;
    vector[n_nodes] ret;
    vector[n-1] zs;
    matrix[n-1,n-1] M = rep_matrix(0, n-1, n-1);

    n = n_nodes-1;

    for (i in 1:n)
      {
        hi[i] = nodes[i+1]-nodes[i];
        bi[i] =  1/hi[i]*(vals[i+1]-vals[i]);
      }
    for (i in 2:n)
      {
        vi[i-1] = 2*(hi[i-1]+hi[i]);
        ui[i-1] = 6*(bi[i] - bi[i-1]);
      }
    for (i in 1:n-1)
      {
        M[i,i] = vi[i];
      }
    for (i in 1:n-2)
      {
        M[i+1,i] = hi[i+1];
        M[i,i+1] = hi[i+1];
      }
    zs = M \ ui ;
    
    ret[1]=0;
    ret[n_nodes] =0;
    ret[2:n_nodes-1]=zs;

    return ret;

  }
  // Evaluate the spline, given nodes, values at the nodes
  // spline coefficients, locations of evaluation points
  // and integer bin ids of each point            
  vector spline_eval(vector nodes,
                     vector vals, vector zs,
                     vector x, array[] int i)
  {
    int n_nodes = size(nodes);
    int n_dat = size(x);
    vector[n_nodes-1] h;
    vector[n_dat] ret;
    array[n_dat] int i1;
    for (ii in 1:n_dat)
      {
        i1[ii] = i[ii] + 1;
      }
    h = spline_geths(nodes);

    ret = (
           zs[i1] ./ 6 ./ h[i] .* square(x-nodes[i]) .* (x-nodes[i])+
           zs[i]  ./ 6 ./ h[i] .* square(nodes[i1]-x) .* (nodes[i1]-x)+
           (vals[i1] ./ h[i] - h[i] .* zs[i1] ./ 6) .* (x-nodes[i])+
           (vals[i] ./ h[i] - h[i] .* zs[i] ./ 6) .* (nodes[i1]-x)
           );
    return ret;
  }
  // find in which node interval we should place each point of the vector                   
array[] int spline_findpos(vector nodes, vector x)
  {
    int n_nodes = size(nodes);
    int n_dat = size(x);
    array[n_dat] int ret;
    for (i in 1:n_dat)
      {
	int success = 0;
        for (j in 1:n_nodes-1)
          {
            if ((x[i]>=nodes[j]) && (x[i]<nodes[j+1]))
              {
                ret[i] = j;
		success = 1;
		break;
              }
          }
        if (success==0)
	  {
	    reject("Point outside knot");
	  }
      }
    return ret;
  }
}
data{
  int<lower=1> N;
  int<lower=1> n_obs;
  int<lower=1> n_beta;
  int<lower=1> n_area;
  array [n_obs] int ind_id;   // index of subjects
  array [n_obs] int ind_area; // index of brain areas
  matrix [n_obs, n_beta] ind_pred; // predictor variable matrix
  vector[n_obs] mri; // mri measurements, sorted by subject
  vector[n_obs] age;
  // input data for controlling spline function
  int n_knots;
  vector [n_knots] age_knots;
}
transformed data{
	// log transform the mri data
	vector[n_obs] log_mri = log(mri);
	
	// determine which knots the point belong to
	array[n_obs] int age_pos_knots= spline_findpos(age_knots, age);
}
parameters{
  real alpha;
  vector[n_knots] knot_values;
  vector[N] id_icpt;
  real<lower=0> var_area_icpt;
  ordered[n_area] area_icpt;
  vector[n_beta] beta;
  real<lower=0> sigma;
}
transformed parameters{
	// these are the spline coefficients corresponding to the current model
	vector[n_knots] spl_coeffs = spline_getcoeffs(age_knots, knot_values);
}
model{
  vector[n_obs] mri_mod;
  
  // priors
  alpha ~ normal(0,10);
  knot_values ~ normal(0,10);
  id_icpt ~ normal(0,10);
  area_icpt ~ normal(0, var_area_icpt);
  var_area_icpt ~ student_t(3,0,10);
  beta ~ normal(0,10);
  sigma ~ student_t(3,0,1);
  
  // compute model predictions combining age spline, varying effects of 
  // brain area and individual, and regression coefficients

  mri_mod = spline_eval(age_knots, knot_values, spl_coeffs, age, age_pos_knots) .*
  
             area_icpt[ind_area] .* id_icpt[ind_id] + 
             
             ind_pred * beta;
  
  log_mri ~ normal(mri_mod, sigma); // likelihood
}
generated quantities{
  array [n_obs] real ppc;
  vector[N] log_lik;
  vector[n_obs] mri_mod;
  
   mri_mod = rep_vector(alpha, n_obs) + 
   
             spline_eval(age_knots, knot_values, spl_coeffs, age, age_pos_knots) +
             
             area_icpt[ind_area] + id_icpt[ind_id] + 
             
             ind_pred * beta;
  
  ppc = normal_rng(mri_mod, rep_array(sigma, n_obs)); 
  
  {int ind = 0; // intialise count variable for segmenting mri data 
                // and model predictions, to calculate log-likelihood per subject
  
    for(n in 1:N){
      
      log_lik[n] = normal_lpdf(segment(log_mri, ind+1, n_area)| segment(mri_mod, ind+1, n_area), sigma);
      
      ind += n_area;}
  }
}
