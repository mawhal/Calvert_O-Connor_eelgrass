data {
  int<lower = 1> n_sites; // number of random effect levels (species) 

  
 // Phenology
  int<lower = 1> Nmacro; // Sample size for phenology data
  array[Nmacro] int<lower = 1, upper = n_sites> macro_sites; // id of random effect (species)
  vector[Nmacro] yMacroi; // Outcome phenology
  vector[Nmacro] Var1i; 
  vector[Nmacro] Var2i;
  vector[Nmacro] Var3i;
  vector[Nmacro] Var4i;
  vector[Nmacro] Var5i;
  vector[Nmacro] Var6i;
  vector[Nmacro] Var7i;
  vector[Nmacro] Var8i;
}

parameters{

  array[n_sites] real alphaVar1Sp;
  real muVar1Sp;
  real<lower = 0> sigmaVar1Sp;

  array[n_sites] real alphaVar2Sp;
  real muVar2Sp;
  real<lower = 0> sigmaVar2Sp;

   array[n_sites] real alphaVar3Sp;
  real muVar3Sp;
  real<lower = 0> sigmaVar3Sp;

   array[n_sites] real alphaVar4Sp;
  real muVar4Sp;
  real<lower = 0> sigmaVar4Sp;

   array[n_sites] real alphaVar5Sp;
  real muVar5Sp;
  real<lower = 0> sigmaVar5Sp;

   array[n_sites] real alphaVar6Sp;
  real muVar6Sp;
  real<lower = 0> sigmaVar6Sp;

   array[n_sites] real alphaVar7Sp;
  real muVar7Sp;
  real<lower = 0> sigmaVar7Sp;

   array[n_sites] real alphaVar8Sp;
  real muVar8Sp;
  real<lower = 0> sigmaVar8Sp;

  array[n_sites] real alphaMacroSites;
  real muMacroSites;
  real<lower = 0> sigmaMacroSites;

  real<lower = 0> sigmaMacro_y;
 }
//

model{
  
 // Phenology likelihood
 for (i in 1:Nmacro){
   yMacroi[i] ~ normal(alphaMacroSites[macro_sites[i]] 
  + alphaVar1Sp[macro_sites[i]] * Var1i[i] 
  + alphaVar2Sp[macro_sites[i]] * Var2i[i]
   + alphaVar3Sp[macro_sites[i]] * Var3i[i]
  + alphaVar4Sp[macro_sites[i]] * Var4i[i]
  + alphaVar5Sp[macro_sites[i]] * Var5i[i]
 + alphaVar6Sp[macro_sites[i]] * Var6i[i]
 + alphaVar7Sp[macro_sites[i]] * Var7i[i]
 + alphaVar8Sp[macro_sites[i]] * Var8i[i]
   , sigmaMacro_y);
 }
 alphaMacroSites ~ normal(muMacroSites, sigmaMacroSites);
alphaVar1Sp ~ normal(muVar1Sp, sigmaVar1Sp);
alphaVar2Sp ~ normal(muVar2Sp, sigmaVar2Sp);
 alphaVar3Sp ~ normal(muVar3Sp, sigmaVar3Sp);
 alphaVar4Sp ~ normal(muVar4Sp, sigmaVar4Sp);
 alphaVar5Sp ~ normal(muVar5Sp, sigmaVar5Sp);
 alphaVar6Sp ~ normal(muVar6Sp, sigmaVar6Sp);
 alphaVar7Sp ~ normal(muVar7Sp, sigmaVar7Sp);
 alphaVar8Sp ~ normal(muVar8Sp, sigmaVar8Sp);


 //// priors
 muMacroSites ~ normal(50,35);
 sigmaMacroSites ~ normal(0,15);

 sigmaMacro_y ~ normal(0,15);

 muVar1Sp ~ normal(0,15);
 sigmaVar1Sp ~ normal(0,15);
 
 muVar2Sp ~ normal(0,15);
 sigmaVar2Sp ~ normal(0,15);

  muVar3Sp ~ normal(0,15);
 sigmaVar3Sp ~ normal(0,15);

  muVar4Sp ~ normal(0,15);
 sigmaVar4Sp ~ normal(0,15);

 muVar5Sp ~ normal(0,15);
sigmaVar5Sp ~ normal(0,15);

 muVar6Sp ~ normal(0,15);
sigmaVar6Sp ~ normal(0,15);

 muVar7Sp ~ normal(0,15);
sigmaVar7Sp ~ normal(0,15);

 muVar8Sp ~ normal(0,15);
sigmaVar8Sp ~ normal(0,15);


}


generated quantities {
} 
