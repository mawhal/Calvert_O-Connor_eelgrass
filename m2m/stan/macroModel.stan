// started Dec 6 2024 by Deirdre

// joint model linking microbial diversity to macro diversity and the environmental factors that shape it

data {
   // Macrofauna
  int<lower = 1> n_site2; 
  int<lower = 1> Nmacro; // Sample size for macrofauna data
  array[Nmacro] int<lower = 1, upper = n_site2> macro_sites; // id of random effect (sites)
  vector[Nmacro] yMacroi; // Outcome macro div
  vector[Nmacro] temp; // predictor linking micro to macro
  vector[Nmacro] salinity; // predictor linking micro to macro
  vector[Nmacro] biomass; // predictor linking micro to macro
  vector[Nmacro] depth; // pred macro only
  vector[Nmacro] year2; // pred macro only
  vector[Nmacro] mAlgae; // pred macro only
  vector[Nmacro] bedArea; // pred macro only
  vector[Nmacro] distEdge; // pred macro only
}

parameters{
  //Macro
  real betaYear2;
  
  array[n_site2] real alphaTempSite;
  real muTempSite;
  real<lower = 0> sigmaTempSite;
  array[n_site2] real alphaSalinSite;
  real muSalinSite;
  real<lower = 0> sigmaSalinSite;
  array[n_site2] real alphaBiomassSite;
  real muBiomassSite;
  real<lower = 0> sigmaBiomassSite;
  
  array[n_site2] real alphaMAlgaeSite;
  real muMAlgaeSite;
  real<lower = 0> sigmaMAlgaeSite;
  array[n_site2] real alphaBedAreaSite;
  real muBedAreaSite;
  real<lower = 0> sigmaBedAreaSite;
  array[n_site2] real alphaDepthSite;
  real muDepthSite;
  real<lower = 0> sigmaDepthSite;
  array[n_site2] real alphaDistEdgeSite;
  real muDistEdgeSite;
  real<lower = 0> sigmaDistEdgeSite;
  
  array[n_site2] real alphaMacroSite;
  real muSite2;
  real<lower = 0> sigma_site2;
  
  real betaMicroxTemp;
  real betaMicroxSalin;
  real betaMicroxBiomass;
  real<lower = 0> sigmaMacro_y;
 }
//
transformed parameters{
  // Macro
  array[n_site2] real betaTemp;    //site level temp
  array[n_site2] real betaSalin;    
  array[n_site2] real betaBiomass;   
  array[n_site2] real betaMAlgae;    //site level temp
  array[n_site2] real betaBedArea;    
  array[n_site2] real betaDepth;   
  array[n_site2] real betaDistEdge;   
}

model{
// Macro likelihood
  for (i in 1:Nmacro){
    yMacroi[i] ~ normal(alphaMacroSite[macro_sites[i]] + betaYear2 * year2[i] + betaTemp[macro_sites[i]] * temp[i] + 
                 betaSalin[macro_sites[i]] * salinity[i] + betaBiomass[macro_sites[i]] * biomass[i] +
                 betaMAlgae[macro_sites[i]] * mAlgae[i] + betaBedArea[macro_sites[i]] * bedArea[i] + 
                 betaDepth[macro_sites[i]] * depth[i] + betaDistEdge[macro_sites[i]] * distEdge[i], 
                 sigmaMacro_y);
  }
  
  alphaMacroSite ~ normal(muSite2, sigma_site2);
  alphaTempSite ~ normal(muTempSite, sigmaTempSite);
  alphaSalinSite ~ normal(muSalinSite, sigmaSalinSite);
  alphaBiomassSite ~ normal(muBiomassSite, sigmaBiomassSite);
  
  alphaMAlgaeSite ~ normal(muMAlgaeSite, sigmaMAlgaeSite);
  alphaBedAreaSite ~ normal(muBedAreaSite, sigmaBedAreaSite);
  alphaDepthSite ~ normal(muDepthSite, sigmaDepthSite);
  alphaDistEdgeSite ~ normal(muDistEdgeSite, sigmaDistEdgeSite);
  
  //// priors
  muSite2 ~ normal(40,10);
  sigma_site2 ~ normal(5,5);

  sigmaMacro_y ~ normal(10,5);

  muTempSite ~ normal(-15,10);
  sigmaTempSite ~ normal(5,5);

  muSalinSite ~ normal(-15,10);
  sigmaSalinSite ~ normal(5,5);

  muBiomassSite ~ normal(-15,10);
  sigmaBiomassSite ~ normal(5,5);
  
  muMAlgaeSite ~ normal(40,10);
  sigmaMAlgaeSite ~ normal(5,5);

  muBedAreaSite ~ normal(-15,10);
  sigmaBedAreaSite ~ normal(5,5);

  muDepthSite ~ normal(-15,10);
  sigmaDepthSite ~ normal(5,5);

  muDistEdgeSite ~ normal(-15,10);
  sigmaDistEdgeSite ~ normal(5,5);

  betaMicroxTemp ~ normal(0,1);
  betaMicroxSalin ~ normal(0,1);
  betaMicroxBiomass ~ normal(0,1);

}

generated quantities {
} 
