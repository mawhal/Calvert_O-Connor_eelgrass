data {
  int<lower = 1> n_sites; // number of random effect levels---sites
 // Micro
 int<lower = 1> N; // Sample size for micro data 
  array[N] int<lower = 1, upper = n_sites> micro_sites; // id of random effect (sites)
  
  vector[N] yeari;
  vector[N] yMicroi; // Observed micro diversity
  
 // Macro
  int<lower = 1> Nmacro; 
  array[Nmacro] int<lower = 1, upper = n_sites> macro_sites; 
  vector[Nmacro] yMacroi; // Outcome phenology
  // vector[Nmacro] Var1i; 
  // vector[Nmacro] Var2i;
  // vector[Nmacro] Var3i;
  vector[Nmacro] Var4i;
  vector[Nmacro] Var5i;
  vector[Nmacro] Var6i;
  vector[Nmacro] Var7i;
  vector[Nmacro] Var8i;
}

parameters{
  // Micro
  real mu_grand; // grand mean 
  array[n_sites] real b_muSite;
  
  array[n_sites] real alphaYear;
  real muYear;
  real<lower = 0> sigmaYear;

  real<lower = 0> sigma_Microy; // sd general
  real<lower = 0> sigma_site; // sd site
  
  // Macro
  // array[n_sites] real alphaVar1;
  // real muVar1;
  // real<lower = 0> sigmaVar1;
  // 
  // array[n_sites] real alphaVar2;
  // real muVar2;
  // real<lower = 0> sigmaVar2;
  // 
  //  array[n_sites] real alphaVar3;
  // real muVar3;
  // real<lower = 0> sigmaVar3;

   array[n_sites] real alphaVar4;
  real muVar4;
  real<lower = 0> sigmaVar4;

   array[n_sites] real alphaVar5;
  real muVar5;
  real<lower = 0> sigmaVar5;

   array[n_sites] real alphaVar6;
  real muVar6;
  real<lower = 0> sigmaVar6;

   array[n_sites] real alphaVar7;
  real muVar7;
  real<lower = 0> sigmaVar7;

   array[n_sites] real alphaVar8;
  real muVar8;
  real<lower = 0> sigmaVar8;

  array[n_sites] real alphaMacroSites;
  real muMacroSites;
  real<lower = 0> sigmaMacroSites;

  real<lower = 0> sigmaMacro_y;
  
  // real betaMicroxV1;
  // real betaMicroxV2;
  // real betaMicroxV3;
  // 
 }
 
transformed parameters{
  vector[N] y_hat;
  vector[n_sites] mu_grand_site;
  
  // array[n_sites] real betaV1;    //site level temp
  // array[n_sites] real betaV2;    
  // array[n_sites] real betaV3;  
  
  
    for(i in 1:n_sites){
    mu_grand_site[i] = mu_grand + b_muSite[i];
  }
  

  for (i in 1:N){
    y_hat[i] =   mu_grand  + b_muSite[micro_sites[i]] + alphaYear[micro_sites[i]] * yeari[i]
    ;
  }

// Macrofauna
  // for (ist in 1:n_sites){
  //   betaV1[ist] = alphaVar1[ist] + betaMicroxV1 * (mu_grand_site[ist]);
  // }
  // for (ist in 1:n_sites){
  //   betaV2[ist] = alphaVar2[ist] + betaMicroxV2 * (mu_grand_site[ist]);
  // }
  // for (ist in 1:n_sites){
  //   betaV3[ist] = alphaVar3[ist] + betaMicroxV3 * (mu_grand_site[ist]);
  // }
  // 
}

model{
  
 // Phenology likelihood
 for (i in 1:Nmacro){
   yMacroi[i] ~ normal(alphaMacroSites[macro_sites[i]] 
  // + betaV1[macro_sites[i]] * Var1i[i] 
  // + betaV2[macro_sites[i]] * Var2i[i]
  // + betaV3[macro_sites[i]] * Var3i[i]
  + alphaVar4[macro_sites[i]] * Var4i[i]
  + alphaVar5[macro_sites[i]] * Var5i[i]
 + alphaVar6[macro_sites[i]] * Var6i[i]
 + alphaVar7[macro_sites[i]] * Var7i[i]
 + alphaVar8[macro_sites[i]] * Var8i[i]
   , sigmaMacro_y);
 }
 
 //Micro 
 alphaYear ~ normal(muYear, sigmaYear);
 b_muSite ~ normal(0, sigma_site);
 mu_grand ~ normal(10,10);
 muYear ~ normal(0,15);
 sigmaYear ~ normal(0,15);
 sigma_site ~ normal(4,5);
 sigma_Microy ~ normal(3, 5);
 yMicroi ~ normal(y_hat, sigma_Microy);
  
  //Macro
 alphaMacroSites ~ normal(muMacroSites, sigmaMacroSites);
 // alphaVar1 ~ normal(muVar1, sigmaVar1);
 // alphaVar2 ~ normal(muVar2, sigmaVar2);
 // alphaVar3 ~ normal(muVar3, sigmaVar3);
 alphaVar4 ~ normal(muVar4, sigmaVar4);
 alphaVar5 ~ normal(muVar5, sigmaVar5);
 alphaVar6 ~ normal(muVar6, sigmaVar6);
 alphaVar7 ~ normal(muVar7, sigmaVar7);
 alphaVar8 ~ normal(muVar8, sigmaVar8);


 //// priors
 muMacroSites ~ normal(50,35);
 sigmaMacroSites ~ normal(0,15);

 sigmaMacro_y ~ normal(0,15);

 // muVar1 ~ normal(0,15);
 // sigmaVar1 ~ normal(0,15);
 // 
 // muVar2 ~ normal(0,15);
 // sigmaVar2 ~ normal(0,15);
 // 
 //  muVar3 ~ normal(0,15);
 // sigmaVar3 ~ normal(0,15);

  muVar4 ~ normal(0,15);
 sigmaVar4 ~ normal(0,15);

 muVar5 ~ normal(0,15);
sigmaVar5 ~ normal(0,15);

 muVar6 ~ normal(0,15);
sigmaVar6 ~ normal(0,15);

 muVar7 ~ normal(0,15);
sigmaVar7 ~ normal(0,15);

 muVar8 ~ normal(0,15);
sigmaVar8 ~ normal(0,15);

  // betaMicroxV1 ~ normal(0,5);
  // betaMicroxV2 ~ normal(0,5);
  // betaMicroxV3 ~ normal(0,5);

}


generated quantities {
} 
