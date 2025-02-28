data {
  int<lower = 1> n_sites; // number of random effect levels (sites) 

  int<lower = 1> N; // Sample size for micro data 
  array[N] int<lower = 0> micro_sites; // id of random effect (sites)
  
  vector[N] yeari;
  vector[N] yMicroi; // Observed micro diversity
}

parameters{
  // Traits
  real mu_grand; // grand mean 
  array[n_sites] real b_muSite;
  
  array[n_sites] real alphaYear;
  real muYear;
  real<lower = 0> sigmaYear;

  real<lower = 0> sigma_Microy; // sd general
  real<lower = 0> sigma_site; // sd site
  
 }
//
transformed parameters{
  vector[N] y_hat;

  for (i in 1:N){
    y_hat[i] =   mu_grand  + b_muSite[micro_sites[i]] + alphaYear[micro_sites[i]] * yeari[i]
    ;
  }

}

model{
  alphaYear ~ normal(muYear, sigmaYear);

  b_muSite ~ normal(0, sigma_site);
  mu_grand ~ normal(10,10);
  
  muYear ~ normal(0,15);
  sigmaYear ~ normal(0,15);
  sigma_site ~ normal(4,5);

  sigma_Microy ~ normal(3, 5);
  
  yMicroi ~ normal(y_hat, sigma_Microy);

}


generated quantities {
} 
