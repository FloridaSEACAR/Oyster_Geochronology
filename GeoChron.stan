// Multilevel GeoChron model
// Age and uncertainties at the sample level
// The input data are vectors 'x' (age), of length N.

data {
  int<lower=1> N; // the number of specimens
  int<lower=1> Locality[N]; // maps holes to locality
  int<lower=1> Reef[N];  // maps holes to reefs
  real x[N];
  int<lower=1> nL; // number of localities
  int<lower=1> nR; // number of unique reefs 
  int<lower=1> Reef2Locality[nR]; // maps reefs to localities
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu_region;
  real mu_locality[nL];
  real mu_reef[nR];
  real<lower=0> sigma_hole;
  real<lower=0> sigma_reef;
  real<lower=0> sigma_locality;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
   mu_region ~ normal(0,100); 
  for( ii in 1:nL) 
     mu_locality[ii] ~ normal(0,100); 
     
   for( ii in 1:nR) 
     mu_reef[ii] ~ normal(0,100);
  
  sigma_hole ~ cauchy(0,1); 
  sigma_reef ~ cauchy(0,1); 
  sigma_locality ~ cauchy(0,1);
  

  
  for( ll in 1:nL){
    mu_locality[ll] ~ normal(mu_region, sigma_locality);
     }
  
  for( rr in 1:nR){
    mu_reef[rr] ~ normal(mu_locality[Reef2Locality[rr]], sigma_reef);
     }
     
  for ( ii in 1:N){
    x[ii] ~ normal(mu_reef[Reef[ii]], sigma_hole);

  }
}

