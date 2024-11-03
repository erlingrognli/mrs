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
  int<lower=1> N, n_obs, n_beta, n_region;
  array [3] int n_spl;
  array [n_obs] int ind_id, ind_region, ind_spl;   // indices of subjects, brain areas and splines
  matrix [n_obs, n_beta] ind_pred; // predictor variable matrix
  vector[n_obs] mri, age; // mri measurements and subject age
  // input data for controlling spline function
  int n_knots;
  vector [n_knots] age_knots;
}
transformed data{
	// log transform the mri data
	vector[n_obs] log_mri = log(mri);
	
	// define indices for separate splines for cortical thickness, area and subcortical volume
	array[n_spl[1]] int ind_thick;
	array[n_spl[2]] int ind_area;
	array[n_spl[3]] int ind_volume;
	
	{array [3] int pos = zeros_int_array(3);
	
	  for(n in 1:n_obs){
	    
	    if(ind_spl[n] == 1){
	      
	      ind_thick[pos[1]+1] = n;
	      
	      pos[1] += 1;}
	    
	    else if(ind_spl[n] == 2){
	      
	      ind_area[pos[2]+1] = n;
	      
	      pos[2] += 1;}
	      
	   else if(ind_spl[n] == 3){
	     
	     ind_volume[pos[3]+1] = n;
	     
	     pos[3] += 1;}
	  }}
	
	// determine which knots points belong to
	array[n_spl[1]] int thick_pos_knots = spline_findpos(age_knots, age[ind_thick]);
	array[n_spl[2]] int area_pos_knots = spline_findpos(age_knots, age[ind_area]);
	array[n_spl[3]] int volume_pos_knots = spline_findpos(age_knots, age[ind_volume]);
}
parameters{
  vector[n_knots] thick_knot_values;
  vector[n_knots] area_knot_values;
  vector[n_knots] volume_knot_values;
  real<lower=0> alpha;
  vector[n_beta] beta;
  real<lower=0> sigma;
  vector<multiplier=2.5>[n_region - 1] region_icpt_raw;
  vector<multiplier=0.5>[N] id_icpt;
}
transformed parameters{
	vector[n_obs] age_spline;
	array [3] vector[n_knots] spl_coeffs;
	
	spl_coeffs[1] = spline_getcoeffs(age_knots, thick_knot_values);
	spl_coeffs[2] = spline_getcoeffs(age_knots, area_knot_values);
	spl_coeffs[3] = spline_getcoeffs(age_knots, volume_knot_values);
	
	// compute and combine spline functions of age for thickness, area and volume
  
  age_spline[ind_thick] = spline_eval(age_knots, thick_knot_values, spl_coeffs[1], age[ind_thick], thick_pos_knots);
  age_spline[ind_area] = spline_eval(age_knots, area_knot_values, spl_coeffs[2], age[ind_area], area_pos_knots);
  age_spline[ind_volume] = spline_eval(age_knots, volume_knot_values, spl_coeffs[3], age[ind_volume], volume_pos_knots);

	// set the icpt of the last area to log(1)
  vector[n_region] region_icpt = append_row(region_icpt_raw, zeros_vector(1));
}
model{
  vector[n_obs] mri_mod;
 
  // priors
  thick_knot_values ~ normal(0, 1);
  area_knot_values ~ normal(0, 1);
  volume_knot_values ~ normal(0, 1);
  id_icpt ~ std_normal();
  region_icpt_raw ~ std_normal();
  alpha ~ normal(0, 10);
  beta ~ normal(0, 0.3);
  sigma ~ student_t(4,0,2);
  
  // compute model predictions combining age spline, varying effects of 
  // brain area and individual, and regression coefficients
  
  mri_mod = alpha + age_spline + 
  
            region_icpt[ind_region] + id_icpt[ind_id] + 
             
            ind_pred * beta;
  
  log_mri ~ normal(mri_mod, sigma); // likelihood
}
generated quantities{
  array [n_obs] real ppc;
  vector[n_obs] log_lik;
  vector[n_obs] mri_mod;
  vector[n_beta] beta_exp = exp(beta);
  real sigma_exp = exp(sigma);
  
  mri_mod = alpha + 
   
             age_spline +
             
             region_icpt[ind_region] + id_icpt[ind_id] + 
             
             ind_pred * beta;
  
  ppc = normal_rng(mri_mod, rep_array(sigma, n_obs));
  
  for(n in 1:n_obs){
    
  log_lik[n] = normal_lpdf(log_mri[n] | mri_mod[n], sigma);
  
  }
}
