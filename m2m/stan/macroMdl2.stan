data {
  int<lower = 1> n_site; // number of random effect levels (sites) 

  int<lower = 1> N; // Sample size for micro data 
  array[N] int<lower = 0> micro_site; // id of random effect (sites)
  
  vector[N] yeari;
  vector[N] yMicroi; // Observed micro diversity
  
    vector[N] tempi;
    vector[N] salini;
    vector[N] biomassi;
    vector[N] macroAi;
    vector[N] bedAreai;
    vector[N] depthi;
    vector[N] distEi;
}

parameters{
  // Traits
  real mu_grand; // grand mean 
  array[n_site] real b_muSite;
  
  real b_year;
  
  real<lower = 0> sigma_Microy; // sd general
  real<lower = 0> sigma_site; // sd site
  
    real mu_force;
    real mu_salin;
    real mu_biomass;
    real mu_macroA;
    real mu_bedArea;
    real mu_depth;
    real mu_distE;

 }
//
transformed parameters{
  vector[N] y_hat;

  for (i in 1:N){
    y_hat[i] =   mu_grand + b_year * yeari[i]  + b_muSite[micro_site[i]]
    + mu_force * tempi[i] + mu_salin * salini[i] 
    + mu_biomass * biomassi[i] + mu_macroA * macroAi[i] 
    + mu_bedArea * bedAreai[i] + mu_depth * depthi[i] 
    + mu_distE * distEi[i] 
    ;
  }

}

model{
  b_muSite ~ normal(0, sigma_site);
  mu_grand ~ normal(10,10);
  b_year ~ normal(0,5);
  sigma_site ~ normal(0,10);
  sigma_Microy ~ normal(0,10);
  
  yMicroi ~ normal(y_hat, sigma_Microy);

  mu_force ~  normal(0,10);
  mu_salin ~ normal(0,10);
  mu_biomass ~  normal(0,10);
  mu_macroA ~ normal(0,10);
  mu_bedArea ~  normal(0,10);
  mu_depth ~ normal(0,10);
  mu_distE ~  normal(0,10);

}


generated quantities {
} 
