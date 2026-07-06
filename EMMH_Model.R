rm(list=ls())
setwd("")
source("OrganiseDataframes.R")

library(rstan)
library(MCMCglmm)
library(dplyr)
library(forcats)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


#######################################################################################
### EMMH Model: caterpillar mean timing, height and width as predictors of blue tit ###
###   fitness function optimum timing, maximum fitness and fitness function width   ###
###      (fitness as hurdle model - Binomial and truncated-Generalised Poisson)     ###
#######################################################################################

#Model
write("functions {
      //truncated generalised poisson (Wang and Famoye 1997 parameterisation)
      real tgenpoiss_lpmf(int y, real lambda1, real lambda2) {
      return y*(log(lambda1)-log1p(lambda1*lambda2))+(y-1)*log1p(lambda2*y)-lgamma(y+1)-lambda1*(1+lambda2*y)/(1+lambda1*lambda2)-log1m(exp(-lambda1/(1+lambda1*lambda2)));
      }
      }
      
      data{
      // CATERPILLAR DATA (c)
      int<lower=0> N_c;     // length of cater data
      int<lower=0> y_c[N_c];     // cater abundance data
      vector[N_c] date_c;     // hatch date data
      int<lower=0> N_site;     // number of sites 
      int<lower=0> N_year;     // number of years 
      int<lower=0> N_siteyear;     // number of site-year combinations 
      int<lower=0> N_siteday;     // number of day-site-year combinations   
      int<lower=0> N_treeID;     // number of unique tree ids   
      int<lower=0> N_rec;     // number of recorders
      int<lower=0,upper=N_site> site_id_c[N_c];     // cater site id data 
      int<lower=0,upper=N_year> year_id_c[N_c];     // cater year id data
      int<lower=0,upper=N_siteyear> siteyear_id_c[N_c];     // cater site-year id data  
      int<lower=0,upper=N_siteday> siteday_id[N_c];     // day-site-year data
      int<lower=0,upper=N_treeID> treeID_id[N_c];     // tree id data
      int<lower=0,upper=N_rec> rec_id[N_c];     // recorder data
      
      // FOR CALCULATING SITEYEAR PARAMETER ESTIMATES
      int<lower=0,upper=N_siteyear> unq_sy[N_siteyear];     // ordered site-year combinations
      int<lower=0,upper=N_site> unq_s[N_siteyear];     // corresponding factor level for sites associated with each site-year
      int<lower=0,upper=N_year> unq_y[N_siteyear];     // corresponding factor level for years associated with each site-year
      
      // BIRD BERNOULLI DATA (bb)
      int<lower=0> N_bb;     // length of bird data, all nests 
      int<lower=0,upper=1> y_bb[N_bb];     // fledge 0/1 bb data
      vector[N_bb] date_bb;     // hatch date bb data
      int<lower=0> N_fem_bb;     // number of females in bb
      int<lower=0,upper=N_site> site_id_bb[N_bb];     // bb site id data
      int<lower=0,upper=N_year> year_id_bb[N_bb];     // bb year id data 
      int<lower=0,upper=N_siteyear> siteyear_id_bb[N_bb];     // bb site-year id data
      int<lower=0,upper=N_fem_bb> fem_id_bb[N_bb];     // bb female id data 
      
      //BIRD TRUNCATED POISSON DATA (bp)
      int<lower=0> N_bp;     // length of truncated bird data, 0 fledged nests removed
      int<lower=1> y_bp[N_bp];     // fledge abundance bp data
      vector[N_bp] date_bp;     // hatch date bp data
      int<lower=0> N_site_bp;     // number of sites in bp data
      int<lower=0> N_year_bp;     // number of years in bp data
      int<lower=0> N_siteyear_bp;     // number of site-year combinations in bp data
      int<lower=0> N_fem_bp;     // number of females in bp data 
      int<lower=0,upper=N_site_bp> site_id_bp[N_bp];     // bp site id data 
      int<lower=0,upper=N_year_bp> year_id_bp[N_bp];     // bp year id data 
      int<lower=0,upper=N_siteyear_bp> siteyear_id_bp[N_bp];     // bp site-year data 
      int<lower=0,upper=N_fem_bp> fem_id_bp[N_bp];     // bp female data 
      
      //FOR LINKING TRUNCATED POIS TO CATER AND BERNOULLI RANDOM EFFECTS
      int<lower=0,upper=N_siteyear> unq_sy_bp[N_siteyear_bp];     //bb and c site-year factor level number that corresponds to each bp site-year
      int<lower=0,upper=N_site> unq_s_bp[N_siteyear_bp];     // bb and c site factor level that corresponds to each bp site-year
      int<lower=0,upper=N_year> unq_y_bp[N_siteyear_bp];     // bb and c year factor level that corresponds to each bp site-year
      }
      
      parameters{
      // CATERPILLAR PARAMETERS
      real mu;     //  cater mean timing intercept
      real logsigma;     // cater log width intercept
      real logmax;     // cater log height intercept
      cholesky_factor_corr[3] L_c;     // cholesky decomposition L matrix for cater gaussian parameters across random effects (used for s, y and sy)
      matrix[3, N_site] site_scaled_c;     // cater scaled site random effects, row for each gaussian parameter 
      matrix[3, N_year] year_scaled_c;     // cater scaled year random effects 
      matrix[3, N_siteyear] siteyear_scaled_c;     // cater scaled site-year random effects
      vector[N_siteday] siteday_scaled;     // scaled day-site-year random effects
      vector[N_treeID] treeID_scaled;     //  scaled tree id random effects
      vector[N_rec] rec_scaled;     // scaled recorder random effects 
      vector<lower=0>[3] sd_site_c;     // standard deviations of cater site random effects for each gaussian parameter
      vector<lower=0>[3] sd_year_c;     // standard deviations of cater year random effects for each gaussian parameter 
      vector<lower=0>[3] sd_siteyear_c;     // standard deviations of cater site-year random effects for each gaussian parameter
      real<lower=0> sd_siteday;     // standard deviation of day-site-year random effects   
      real<lower=0> sd_treeID;     // standard deviation of tree id random effects
      real<lower=0> sd_rec;     // standard deviation of recorder random effects
      real<lower=0> sd_resid_c;     // standard deviation of cater observation random effects
      vector<lower=0>[N_c] l_y_1c;     // data scale cater estimates from lognormal dist with resid sd
      
      
      // BIRD BERNOULLI PARAMETERS (bb)
      real theta_bb;     // bb optimum timing intercept (at mu=0)
      real logomega_bb;     // bb log width intercept (at logsigma=0 and logmax=0)
      real logWmax_bb;     // bb logit height intercept (at logmax=0) - code throughout says logWmax for bb but estimated as logit
      real mu_theta_bb;     // bb slope for change in theta with mu
      real ls_logomega_bb;     // bb slope for change in logomega with logsigma   
      real lm_logomega_bb;     // bb slope for change in logomega with logmax
      real lm_logWmax_bb;     // bb slope for change in logWmax with logmax        
      cholesky_factor_corr[3] L_bb;     // cholesky decomp L matrix for bb gaussian parameters across random effects (s,y and sy)
      matrix[3, N_site] site_scaled_bb;     // bb scaled site random effects, row for each parameter 
      matrix[3, N_year] year_scaled_bb;     // bb scaled year random effects 
      matrix[3, N_siteyear] siteyear_scaled_bb;     // bb scaled site-year random effects
      vector[N_fem_bb] fem_scaled_bb;     // bb scaled female random effects
      vector<lower=0>[3] sd_site_bb;     // standard deviations of bb site random effects for each gaussian parameter
      vector<lower=0>[3] sd_year_bb;     // standard deviations of bb year random effects for each gaussian parameter
      vector<lower=0>[3] sd_siteyear_bb;     // standard deviations of bb site-year random effects for each gaussian parameter
      real<lower=0> sd_fem_bb;     // standard deviation of bb female random effects
      
      // REGRESSION PARAMETERS BETWEEN BIRD BERNOULLI AND TRUNCATED POISSON (bbp)
      vector[2] bbp_reg;     // one coefficient for each gaussian function parameter - theta and Wmax (not included for logomega)
      
      // BIRD TRUNCATED POISSON PARAMETERS (bp)
      real theta_bp;     // bp optimum timing intercept (at mu=0)
      real logomega_bp;     // bp log width intercept (at logsigma=0 and logmax=0)
      real logWmax_bp;     // bp log height intercept (at logmax=0)
      real mu_theta_bp;     // bp slope for change in theta with mu
      real ls_logomega_bp;     // bp slope for change in logomega with logsigma   
      real lm_logomega_bp;     // bp slope for change in logomega with logmax
      real lm_logWmax_bp;     // bp slope for change in logWmax with logmax        
      cholesky_factor_corr[3] L_bp;     // cholesky decomp L matrix for bp gaussian parameters across random effects (s,y and sy)
      matrix[3, N_site_bp] site_scaled_bp;     // bp scaled site residuals from regression against bb effects, row for each parameter  
      matrix[3, N_year_bp] year_scaled_bp;     // bp scaled year residuals from regression against bb effects, row for each parameter   
      matrix[3, N_siteyear_bp] siteyear_scaled_bp;     // bp scaled site-year residuals from regression against bb effects, row for each parameter  
      vector[N_fem_bp] fem_scaled_bp;     // bp scaled female random effects
      vector<lower=0>[3] sd_site_bp;     // standard deviation of bp site residuals from regression against bb effects for each gaussian parameter
      vector<lower=0>[3] sd_year_bp;     // standard deviation of bp year residuals from regression against bb effects for each gaussian parameter 
      vector<lower=0>[3] sd_siteyear_bp;     // standard deviation of bp site-year residuals from regression against bb effects for each gaussian parameter
      real<lower=0> sd_fem_bp;     // standard deviation of bp female random effects
      real lambda;     // over- or under-dispersion parameter (number fledged is under-dispersed for a Poisson distribution)
      
      }
      
      transformed parameters{
      matrix[N_site,3] site_effs_c;     // cater unscaled site random effects, column for each parameter
      matrix[N_year,3] year_effs_c;     // cater unscaled year random effects
      matrix[N_siteyear,3] siteyear_effs_c;     // cater unscaled site-year random effects
      
      vector[N_siteyear] c_mu_bb;     // cater timing estimate for each site-year in bb
      vector[N_siteyear] c_ls_bb;     // cater log width estimate for each site-year in bb
      vector[N_siteyear] c_lm_bb;     // cater log height estimate for each site-year in bb
      
      vector[N_siteyear_bp] c_mu_bp;     // cater timing estimate for each site-year in bp
      vector[N_siteyear_bp] c_ls_bp;     // cater log width estimate for each site-year in bp
      vector[N_siteyear_bp] c_lm_bp;     // cater log height estimate for each site-year in bp
      
      matrix[N_site,3] site_effs_bb;     // bb unscaled site random effects, column for each parameter
      matrix[N_year,3] year_effs_bb;     // bb unscaled year random effects
      matrix[N_siteyear,3] siteyear_effs_bb;     // bb unscaled site-year random effects  
      
      matrix[N_site_bp,3] site_bp_r;     // bp unscaled site residuals from regression against bb effects, column for each parameter 
      matrix[N_year_bp,3] year_bp_r;     // bp unscaled year residuals from regression against bb effects
      matrix[N_siteyear_bp,3] siteyear_bp_r;     // bp unscaled site-year residuals from regression against bb effects
      
      matrix[N_site_bp,2] site_bp_p;     // bp predicted site random effects from regression against bb effects, column for each parameter
      matrix[N_year_bp,2] year_bp_p;     // bp predicted year random effects from regression against bb effects
      matrix[N_siteyear_bp,2] siteyear_bp_p;     // bp predicted site-year random effects from regression against bb effects
      
      matrix[N_site_bp,3] site_effs_bp;     // bp site random effects (predicted + residual), column for each parameter
      matrix[N_year_bp,3] year_effs_bp;     // bp year random effects (predicted + residual), column for each parameter
      matrix[N_siteyear_bp,3] siteyear_effs_bp;     // bp site-year random effects (predicted + residual), column for each parameter   
      
      site_effs_c = (diag_pre_multiply(sd_site_c, L_c) * site_scaled_c)';     // sds as a diag matrix multipled by L matrix then by scaled effects and transposed  
      year_effs_c = (diag_pre_multiply(sd_year_c, L_c) * year_scaled_c)';
      siteyear_effs_c = (diag_pre_multiply(sd_siteyear_c, L_c) * siteyear_scaled_c)';
      
      c_mu_bb = mu+site_effs_c[unq_s,1]+year_effs_c[unq_y,1]+siteyear_effs_c[unq_sy,1]-8;     // -8 so 'centred' roughly, sum the appropiate intercept, site, year and siteyear effects for each site-year in bb 
      c_ls_bb = logsigma+site_effs_c[unq_s,2]+year_effs_c[unq_y,2]+siteyear_effs_c[unq_sy,2]-2.5;
      c_lm_bb = logmax+site_effs_c[unq_s,3]+year_effs_c[unq_y,3]+siteyear_effs_c[unq_sy,3]+2.5;
      
      c_mu_bp = mu+site_effs_c[unq_s_bp,1]+year_effs_c[unq_y_bp,1]+siteyear_effs_c[unq_sy_bp,1]-8;     // -8 so 'centred' roughly, sum the appropiate intercept, site, year and siteyear effects for each site-year in bp 
      c_ls_bp = logsigma+site_effs_c[unq_s_bp,2]+year_effs_c[unq_y_bp,2]+siteyear_effs_c[unq_sy_bp,2]-2.5; // -2.5 so 'centred' roughly
      c_lm_bp = logmax+site_effs_c[unq_s_bp,3]+year_effs_c[unq_y_bp,3]+siteyear_effs_c[unq_sy_bp,3]+2.5; // +2.5 so 'centred' roughly
      
      site_effs_bb = (diag_pre_multiply(sd_site_bb, L_bb) * site_scaled_bb)';      // sds as a diag matrix multipled by L matrix then by scaled effects and transposed
      year_effs_bb = (diag_pre_multiply(sd_year_bb, L_bb) * year_scaled_bb)';
      siteyear_effs_bb = (diag_pre_multiply(sd_siteyear_bb, L_bb) * siteyear_scaled_bb)';
      
      site_bp_r = (diag_pre_multiply(sd_site_bp, L_bp) * site_scaled_bp)';     // sds as a diag matrix multipled by L matrix then by scaled effects and transposed 
      year_bp_r = (diag_pre_multiply(sd_year_bp, L_bp) * year_scaled_bp)';
      siteyear_bp_r = (diag_pre_multiply(sd_siteyear_bp, L_bp) * siteyear_scaled_bp)';
      
      site_bp_p = diag_post_multiply(site_effs_bb[,{1,3}],bbp_reg);     // bb effects multipled by diagonal matrix of the regression coefficients 
      year_bp_p = diag_post_multiply(year_effs_bb[,{1,3}],bbp_reg);
      siteyear_bp_p = diag_post_multiply(siteyear_effs_bb[unq_sy_bp,{1,3}],bbp_reg);
      
      site_effs_bp[,{1,3}] = site_bp_r[,{1,3}] + site_bp_p;     // sum predicted and resid
      year_effs_bp[,{1,3}] = year_bp_r[,{1,3}] + year_bp_p;
      siteyear_effs_bp[,{1,3}] = siteyear_bp_r[,{1,3}] + siteyear_bp_p;
      
      site_effs_bp[,2] = site_bp_r[,2];     // no covariance with bb for logomega effects so residual is full effect
      year_effs_bp[,2] = year_bp_r[,2];
      siteyear_effs_bp[,2] = siteyear_bp_r[,2];
      
      }
      
      model{
      //Caterpillar
      vector[N_c] y_1c;     // logscale cater estimates for mean of normal dist (with estimated residual sd)
      vector[N_siteday] siteday_effs = sd_siteday * siteday_scaled;     // unscale random effects that have no covariance
      vector[N_treeID] treeID_effs = sd_treeID * treeID_scaled; 
      vector[N_rec] rec_effs = sd_rec * rec_scaled; 
      
      // Caterpillar Gaussian parameter calculations
      vector[N_c] mu_obs = mu+site_effs_c[site_id_c,1]+year_effs_c[year_id_c,1]+siteyear_effs_c[siteyear_id_c,1];
      vector[N_c] sigma_obs = exp(logsigma+site_effs_c[site_id_c,2]+year_effs_c[year_id_c,2]+siteyear_effs_c[siteyear_id_c,2]);
      vector[N_c] logmax_obs = logmax + site_effs_c[site_id_c,3] + year_effs_c[year_id_c,3] + siteyear_effs_c[siteyear_id_c,3] 
      + siteday_effs[siteday_id] + treeID_effs[treeID_id] + rec_effs[rec_id];
      
      //Bird Bernoulli
      vector[N_bb] y_1bb;     // logit scale bb estimates
      vector[N_fem_bb] fem_effs_bb = sd_fem_bb * fem_scaled_bb;     // unscale random effects that have no covariance
      
      // Bird Bernoulli Gaussian parameter calculations
      vector[N_bb] theta_bb_obs = theta_bb+mu_theta_bb*c_mu_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,1]+year_effs_bb[year_id_bb,1]+siteyear_effs_bb[siteyear_id_bb,1];
      vector[N_bb] omega_bb_obs = exp(logomega_bb+ls_logomega_bb*c_ls_bb[siteyear_id_bb]+lm_logomega_bb*c_lm_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,2]+year_effs_bb[year_id_bb,2]+siteyear_effs_bb[siteyear_id_bb,2]);
      vector[N_bb] logWmax_bb_obs = logWmax_bb + lm_logWmax_bb*c_lm_bb[siteyear_id_bb] 
      + site_effs_bb[site_id_bb,3] + year_effs_bb[year_id_bb,3] + siteyear_effs_bb[siteyear_id_bb,3] 
      + fem_effs_bb[fem_id_bb];
      
      // Bird truncated-Poisson
      vector[N_bp] y_1bp;     // logscale bird trunc-Pois estimates for mean of normal dist
      vector[N_fem_bp] fem_effs_bp = sd_fem_bp * fem_scaled_bp;     // unscale random effects that have no covariance
      
      // Bird truncated-Poisson Gaussian parameter calculations
      vector[N_bp] theta_bp_obs = theta_bp+mu_theta_bp*c_mu_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,1]+year_effs_bp[year_id_bp,1]+siteyear_effs_bp[siteyear_id_bp,1];
      vector[N_bp] omega_bp_obs = exp(logomega_bp+ls_logomega_bp*c_ls_bp[siteyear_id_bp]+lm_logomega_bp*c_lm_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,2]+year_effs_bp[year_id_bp,2]+siteyear_effs_bp[siteyear_id_bp,2]);
      vector[N_bp] logWmax_bp_obs = logWmax_bp + lm_logWmax_bp*c_lm_bp[siteyear_id_bp] 
      + site_effs_bp[site_id_bp,3] + year_effs_bp[year_id_bp,3] + siteyear_effs_bp[siteyear_id_bp,3] 
      + fem_effs_bp[fem_id_bp];
      
      // CATERPILLAR MODEL
      y_1c = logmax_obs - square(date_c-mu_obs) ./ (2 * square(sigma_obs));
      l_y_1c ~ lognormal(y_1c, sd_resid_c);
      y_c ~ poisson(l_y_1c);
      
      // BIRD BERNOULLI MODEL
      y_1bb = logWmax_bb_obs - square(date_bb-theta_bb_obs) ./ (2 * square(omega_bb_obs));
      y_bb ~ bernoulli_logit(y_1bb);
      
      // BIRD TRUNCATED-POISSON MODEL
      y_1bp = exp(logWmax_bp_obs - square(date_bp-theta_bp_obs) ./ (2 * square(omega_bp_obs)));
      
      for(i in 1:N_bp)     
      target += tgenpoiss_lpmf(y_bp[i] | y_1bp[i], lambda);
      //print(omega_bp_obs);
      
      // CATERPILLAR PRIORS      
      to_vector(site_scaled_c) ~ normal(0,1); 
      to_vector(year_scaled_c) ~ normal(0,1); 
      to_vector(siteyear_scaled_c) ~ normal(0,1); 
      siteday_scaled ~ normal(0,1); 
      treeID_scaled ~ normal(0,1); 
      rec_scaled ~ normal(0,1); 
      mu ~ normal(0,20);
      logsigma ~ normal(0,5);
      logmax ~ normal(0,10);
      L_c ~ lkj_corr_cholesky(2.0);  
      to_vector(sd_site_c) ~ cauchy(0,10);
      to_vector(sd_year_c) ~ cauchy(0,10);
      to_vector(sd_siteyear_c) ~ cauchy(0,10);
      sd_siteday ~ cauchy(0,10);
      sd_treeID ~ cauchy(0,10);
      sd_rec ~ cauchy(0,10);
      sd_resid_c ~ cauchy(0,10);
      
      // BIRD BERNOULLI PRIORS
      to_vector(site_scaled_bb) ~ normal(0,1); 
      to_vector(year_scaled_bb) ~ normal(0,1); 
      to_vector(siteyear_scaled_bb) ~ normal(0,1); 
      fem_scaled_bb ~ normal(0,1);
      theta_bb ~ normal(0,20);
      logomega_bb ~ normal(0,5);
      logWmax_bb ~ normal(0,sqrt(square(pi())/3));
      mu_theta_bb ~ normal(0,5);
      ls_logomega_bb ~ normal(0,5);
      lm_logomega_bb ~ normal(0,5);
      lm_logWmax_bb ~ normal(0,5);       
      L_bb ~ lkj_corr_cholesky(2.0);  
      sd_site_bb[1] ~ cauchy(0,10);
      sd_site_bb[2] ~ cauchy(0,3);
      sd_site_bb[3] ~ cauchy(0,3);
      sd_year_bb[1] ~ cauchy(0,10);
      sd_year_bb[2] ~ cauchy(0,3);
      sd_year_bb[3] ~ cauchy(0,3);
      sd_siteyear_bb[1] ~ cauchy(0,10);
      sd_siteyear_bb[2] ~ cauchy(0,3);
      sd_siteyear_bb[3] ~ cauchy(0,3);
      sd_fem_bb ~ cauchy(0,10);
      
      to_vector(bbp_reg) ~ normal(0,10);
      
      // BIRD TRUNCATED POISSON PRIORS
      to_vector(site_scaled_bp) ~ normal(0,1); 
      to_vector(year_scaled_bp) ~ normal(0,1); 
      to_vector(siteyear_scaled_bp) ~ normal(0,1); 
      fem_scaled_bp ~ normal(0,1);
      theta_bp ~ normal(0,20);
      logomega_bp ~ normal(0,5);
      logWmax_bp ~ normal(0,10);
      mu_theta_bp ~ normal(0,5);
      ls_logomega_bp ~ normal(0,5);
      lm_logomega_bp ~ normal(0,5);
      lm_logWmax_bp ~ normal(0,5);       
      L_bp ~ lkj_corr_cholesky(2.0);  
      sd_site_bp[1] ~ cauchy(0,10);
      sd_site_bp[2] ~ cauchy(0,3);
      sd_site_bp[3] ~ cauchy(0,3);
      sd_year_bp[1] ~ cauchy(0,10);
      sd_year_bp[2] ~ cauchy(0,3);
      sd_year_bp[3] ~ cauchy(0,3);
      sd_siteyear_bp[1] ~ cauchy(0,10);
      sd_siteyear_bp[2] ~ cauchy(0,3);
      sd_siteyear_bp[3] ~ cauchy(0,3);
      sd_fem_bp ~ cauchy(0,10);
      lambda ~ normal(0,10);
      }    
      
      generated quantities {
      
      matrix[3, 3] omega_c;  // correlation matrix = L times L transposed
      matrix[3, 3] omega_bb;
      matrix[3, 3] omega_bp;
      
      omega_c = L_c * L_c'; 
      omega_bb = L_bb * L_bb'; 
      omega_bp = L_bp * L_bp'; 
      
      }
      ",
      
      "FullMMH.stan"
)
stanc("FullMMH.stan")
FullMMH_stan <- stan_model("FullMMH.stan")

data_FullMMH <- list(
  N_c=nrow(cater),
  y_c=cater$caterpillars,
  date_c=cater$date_cent,
  N_site=length(unique(cater$site)),
  N_siteyear=length(unique(cater$siteyear)),
  N_year=length(unique(cater$year)),
  N_siteday=length(unique(cater$siteday)),
  N_treeID=length(unique(cater$treeID)),
  N_rec=length(unique(cater$recorder)),
  site_id_c=as.numeric(as.factor(cater$site)),
  siteyear_id_c=as.numeric(as.factor(cater$siteyear)),
  year_id_c=as.numeric(as.factor(cater$year)),
  siteday_id=as.numeric(as.factor(cater$siteday)),
  treeID_id=as.numeric(as.factor(cater$treeID)),
  rec_id=as.numeric(as.factor(cater$recorder)),
  unq_sy=as.numeric(as.factor(SY_id$siteyear_id)),
  unq_s=as.numeric(as.factor(SY_id$site_id)),
  unq_y=as.numeric(as.factor(SY_id$year_id)),
  unq_sy_bp=as.numeric(SY_id_p$siteyear_id_c),
  unq_s_bp=as.numeric(SY_id_p$site_id_c),
  unq_y_bp=as.numeric(SY_id_p$year_id_c),
  N_bb=nrow(SYBhd),
  y_bb=SYBhd$fledged.binom,
  date_bb=SYBhd$hd_cent,
  N_fem_bb=length(unique(SYBhd$fem)),
  site_id_bb=as.numeric(as.factor(SYBhd$site)),
  siteyear_id_bb=as.numeric(as.factor(SYBhd$siteyear)),
  year_id_bb=as.numeric(as.factor(SYBhd$year)),
  fem_id_bb=as.numeric(as.factor(SYBhd$fem)),
  N_bp=nrow(SYBhd_p),
  y_bp=SYBhd_p$fledged,
  date_bp=SYBhd_p$hd_cent,
  N_site_bp=length(unique(SYBhd_p$site)),
  N_siteyear_bp=length(unique(SYBhd_p$siteyear)),
  N_year_bp=length(unique(SYBhd_p$year)),
  N_fem_bp=length(unique(SYBhd_p$fem)),
  site_id_bp=as.numeric(as.factor(SYBhd_p$site)),
  siteyear_id_bp=as.numeric(as.factor(SYBhd_p$siteyear)),
  year_id_bp=as.numeric(as.factor(SYBhd_p$year)),
  fem_id_bp=as.numeric(as.factor(SYBhd_p$fem))
)

#Running model

init=list(list(logomega_bp=log(runif(1,5,15)), ls_logomega_bp=runif(1,-0.1,0.1), lm_logomega_bp=runif(1,-0.1,0.1), site_scaled_bp=matrix(runif(3*44,-0.1,0.1),3,44), year_scaled_bp=matrix(runif(3*12,-0.1,0.1),3,12), siteyear_scaled_bp=matrix(runif(3*381,-0.1,0.1),3,381), bbp_reg=runif(2,-0.1,0.1)),
          list(logomega_bp=log(runif(1,5,15)), ls_logomega_bp=runif(1,-0.1,0.1), lm_logomega_bp=runif(1,-0.1,0.1), site_scaled_bp=matrix(runif(3*44,-0.1,0.1),3,44), year_scaled_bp=matrix(runif(3*12,-0.1,0.1),3,12), siteyear_scaled_bp=matrix(runif(3*381,-0.1,0.1),3,381), bbp_reg=runif(2,-0.1,0.1)),
          list(logomega_bp=log(runif(1,5,15)), ls_logomega_bp=runif(1,-0.1,0.1), lm_logomega_bp=runif(1,-0.1,0.1), site_scaled_bp=matrix(runif(3*44,-0.1,0.1),3,44), year_scaled_bp=matrix(runif(3*12,-0.1,0.1),3,12), siteyear_scaled_bp=matrix(runif(3*381,-0.1,0.1),3,381), bbp_reg=runif(2,-0.1,0.1)),
          list(logomega_bp=log(runif(1,5,15)), ls_logomega_bp=runif(1,-0.1,0.1), lm_logomega_bp=runif(1,-0.1,0.1), site_scaled_bp=matrix(runif(3*44,-0.1,0.1),3,44), year_scaled_bp=matrix(runif(3*12,-0.1,0.1),3,12), siteyear_scaled_bp=matrix(runif(3*381,-0.1,0.1),3,381), bbp_reg=runif(2,-0.1,0.1)))

FullMMH_mod <- sampling(object=FullMMH_stan, data=data_FullMMH,
                                chains=4, cores=4, warmup= 2000, iter=24000, thin=20,
                                init=init, control = list(adapt_delta = 0.94,stepsize = 0.1, max_treedepth = 12))  
saveRDS(FullMMH_mod,"FullMMH_mod.rds")

mod <- FullMMH_mod

rstan::traceplot(mod, pars=c("mu", "logsigma", "logmax",
                             "theta_bb", "logWmax_bb", "logomega_bb",
                             "theta_bp", "logWmax_bp", "logomega_bp","mu_theta_bb","ls_logomega_bb","lm_logomega_bb","lm_logWmax_bb",
                             "mu_theta_bp","ls_logomega_bp","lm_logomega_bp","lm_logWmax_bp"))

rstan::traceplot(mod, pars=c("sd_site_bb","sd_year_bb","sd_siteyear_bb","sd_site_bp","sd_year_bp","sd_siteyear_bp","bbp_reg"))

