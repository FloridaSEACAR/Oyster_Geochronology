//CZ325 Task 3: Evaluation of sample size to assess managed-area-level trends Deliverable 3a:

// random effects linear model 
// where means have a geostatistical model
 

data {
  int<lower=1> N;
  vector[N] y; // shell heights
  int<lower=1> G; // number of groups/sites
  int group[N];  // map data to group
  vector[2] positions[G];   // distance matrix (in km)
  //int M; // number of new positions
 //vector[2] new_positions[M];
}

transformed data {
  // Small offset to ensure the covariance matrix is positive definite
  real delta = 1e-9;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[G] global_mean;

  
  real<lower=0> sigma;

  vector[G] site_mean;

  // these are for the gaussian process  
  real<lower=0> length_scale;
  real<lower=0> gp_sigma;
  vector[G] eta;

  //vector[N] w; 
  real<lower=0> sig;
  
}

transformed parameters {
  // Calculate the latent Gaussian process function (phi)
  matrix[G,G] L_chol = gp_exp_quad_cov( positions,gp_sigma,length_scale) + diag_matrix(rep_vector(delta, G));
  vector[G] gp_function = cholesky_decompose(L_chol) * eta;
}



model {
  
  length_scale ~ normal(0, 2.5);
  // Prior for the GP covariance magnitude
  gp_sigma ~ std_normal();
  // Multiplier for non-centred GP parameterisation
  eta ~ std_normal();
  sig ~ cauchy(0,1);

  global_mean ~ normal(0,10);
  
  
  site_mean ~ normal(global_mean + gp_function, sig);
  for(nn in 1:N){
    y[nn] ~ normal(site_mean[group[nn]], sigma);
    // y[nn] ~ normal(global_mean[group[nn]] + gp_function[group[nn]], sigma);
  }
}






