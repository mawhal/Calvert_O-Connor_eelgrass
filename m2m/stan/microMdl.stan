data {
  int<lower = 1> n_site; // number of random effect levels (species) 
 
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  array[N] int<lower = 0> micro_site; // id of random effect (species)
  
  vector[N] yeari;
  vector[N] yMicroi; // Observed trait
}

parameters{
  // Traits
  real mu_grand; // grand mean for trait value 
  array[n_site] real b_muSite;
  
  real b_year;
  
  //real<lower = 0> sigma_tranlat; // sd pop
  real<lower = 0> sigma_Microy; // sd general
  real<lower = 0> sigma_site; // sd species
  
 }
//
transformed parameters{
  // Traits
  vector[N] y_hat;
  

  // 
  // // Trait

  for (i in 1:N){
    y_hat[i] =   mu_grand + b_year * yeari[i]  + b_muSite[micro_site[i]]
    ;
  }


}

model{
  // Traits
  //// likelihood
  b_muSite ~ normal(0, sigma_site);
  mu_grand ~ normal(10,10);
  b_year ~ normal(0,5);
 // sigma_tranE ~ normal(0,10);
  
  //// priors
 // mu_grand ~ normal(10,10);
  sigma_site ~ normal(4,5);
  sigma_Microy ~ normal(3, 5);
  
  yMicroi ~ normal(y_hat, sigma_Microy);


}


generated quantities {
} 
