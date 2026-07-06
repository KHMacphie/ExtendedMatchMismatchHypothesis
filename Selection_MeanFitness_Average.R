rm(list=ls())
setwd("")
source("OrganisingDataframes.R")
source("ResultsFunctions.R")
library(rstan)
library(rstanarm)
library(MCMCglmm)
library(dplyr)
library(forcats)
library(mvtnorm)

mod <- readRDS("FullMMH_mod.rds")


## Posteriors
cater_mu_bb <- stanpost(model=mod, parameters=c("c_mu_bb"))
cater_ls_bb <- stanpost(model=mod, parameters=c("c_ls_bb"))
cater_lm_bb <- stanpost(model=mod, parameters=c("c_lm_bb"))

th_bb <- stanpost(model=mod, parameters=c("theta_bb"))
lo_bb <- stanpost(model=mod, parameters=c("logomega_bb"))
lw_bb <- stanpost(model=mod, parameters=c("logWmax_bb"))

mu_th_bb <- stanpost(model=mod, parameters=c("mu_theta_bb"))
ls_lo_bb <- stanpost(model=mod, parameters=c("ls_logomega_bb"))
lm_lo_bb<- stanpost(model=mod, parameters=c("lm_logomega_bb"))
lm_lw_bb<- stanpost(model=mod, parameters=c("lm_logWmax_bb"))

bb_s_sd <- stanpost(model=mod, parameters=c("sd_site_bb"))
bb_y_sd <- stanpost(model=mod, parameters=c("sd_year_bb"))
bb_sy_sd <- stanpost(model=mod, parameters=c("sd_siteyear_bb"))
bb_fem_sd <- stanpost(model=mod, parameters=c("sd_fem_bb"))

th_bp <- stanpost(model=mod, parameters=c("theta_bp"))
lo_bp <- stanpost(model=mod, parameters=c("logomega_bp"))
lw_bp <- stanpost(model=mod, parameters=c("logWmax_bp"))

mu_th_bp <- stanpost(model=mod, parameters=c("mu_theta_bp"))
ls_lo_bp <- stanpost(model=mod, parameters=c("ls_logomega_bp"))
lm_lo_bp <- stanpost(model=mod, parameters=c("lm_logomega_bp"))
lm_lw_bp <- stanpost(model=mod, parameters=c("lm_logWmax_bp"))

bp_s_sd <- stanpost(model=mod, parameters=c("sd_site_bp"))
bp_y_sd <- stanpost(model=mod, parameters=c("sd_year_bp"))
bp_sy_sd <- stanpost(model=mod, parameters=c("sd_siteyear_bp"))
bp_fem_sd <- stanpost(model=mod, parameters=c("sd_fem_bp"))

omega_bb <- stanpost(model=mod, parameters=c("omega_bb"))
omega_bp <- stanpost(model=mod, parameters=c("omega_bp"))

## Hatch date model - run in Lag_Prediction
load(file="hd_mod")

# store trait mean and sd
hd_mean_post <- posterior_samples(hd_mod,pars=c("b_Intercept"))
hd_sd_post <- posterior_samples(hd_mod,pars=c("sigma"))

N.h <- 10000 

meanfit_dirsel_av_pred <- function(N,link, hd_mean, hd_sd, mu, logsigma, logmax, omega,s_sd,y_sd,sy_sd,fem_sd,theta,logomega,lWmax,mu_th,ls_lo,lm_lo,lm_lw){
  
  pb <- txtProgressBar(min = 1, max = nrow(theta), style = 3)
  
  mf_ds <- data.frame(matrix(NA,nrow=2,ncol=4))
  colnames(mf_ds) <- c("metric","mean","lci","uci")
  
  mf <- c()
  ds <- c()
  
  # Simulate predictions from the model
  for (i in 1:nrow(theta)){ # for each iteration
    
    Esim <- c()
    
    dat.cnt <- rnorm(N, mean=hd_mean[i,1], sd=hd_sd[i,1])
    
    mu_mn=mean(as.numeric(mu[i,]))
    logsigma_mn=mean(as.numeric(logsigma[i,]))
    logmax_mn=mean(as.numeric(logmax[i,]))
    
    cor.matrix <- matrix(c(omega[i,1], omega[i,2], omega[i,3],
                           omega[i,4], omega[i,5], omega[i,6],
                           omega[i,7], omega[i,8], omega[i,9]), nrow=3, ncol=3) #correlation matrix- called omega in model code
    
    site.sd <- diag(c(s_sd[i,1], s_sd[i,2], s_sd[i,3])) # for theta, logomega and logWmax respectively
    year.sd <- diag(c(y_sd[i,1], y_sd[i,2], y_sd[i,3]))
    styr.sd <- diag(c(sy_sd[i,1], sy_sd[i,2], sy_sd[i,3]))
    
    site.covar <- site.sd%*%cor.matrix%*%site.sd  
    year.covar <- year.sd%*%cor.matrix%*%year.sd
    styr.covar <- styr.sd%*%cor.matrix%*%styr.sd
    
    site.ef <- rmvnorm(n=N, mean=c(0,0,0), sigma=site.covar) #rmvnorm random number generator for multivariate with VCV
    year.ef <- rmvnorm(n=N, mean=c(0,0,0), sigma=year.covar)
    styr.ef <- rmvnorm(n=N, mean=c(0,0,0), sigma=styr.covar)
    
    fem.ef <- rnorm(N.h, mean=0, sd=fem_sd[i,1])
    
    f1 <- lWmax[i,1] + lm_lw[i,1]*logmax_mn + site.ef[,3] + year.ef[,3] + styr.ef[,3]
    f2 <- dat.cnt - theta[i,1] - mu_th[i,1]*mu_mn - site.ef[,1] - year.ef[,1] - styr.ef[,1]
    f3 <- sqrt(2) * exp(logomega[i,1] + ls_lo[i,1]*logsigma_mn + lm_lo[i,1]*logsigma_mn + site.ef[,2] + year.ef[,2] + styr.ef[,2])
    
    f4 <- fem.ef
    
    if(link=="logit"){
      Esim <- plogis(f1+f4-(f2/f3)^2)
    } else if(link=="log"){
      Esim <- exp(f1+f4-(f2/f3)^2)
    } 
    
    dat.cnt.cnt <- dat.cnt-mean(dat.cnt)
    
    EW.rel <- Esim/mean(Esim)
    
    mod <- lm(EW.rel~dat.cnt.cnt+I(dat.cnt.cnt^2))
    
    
    mf[i] <- mean(Esim)
    ds[i] <- as.numeric(coefficients(mod)[2])
    
    f1 <- NULL
    f2 <- NULL
    f3 <- NULL
    f4 <- NULL
    
    setTxtProgressBar(pb, i)
    
  } 
  
  mf_ds[1,] <- c("Mean fitness", mean(mf), posterior_interval(matrix(mf), prob=0.95))
  mf_ds[2,] <- c("Directional selection", mean(ds), posterior_interval(matrix(ds), prob=0.95))
  
  close(pb)  
  return(mf_ds)
}



bb_pred <- meanfit_dirsel_av_pred(N=N.h,
                                  link="logit",
                                  hd_mean=hd_mean_post, 
                                  hd_sd=hd_sd_post,
                                  omega=omega_bb,
                                  s_sd=bb_s_sd,
                                  y_sd=bb_y_sd,
                                  sy_sd=bb_sy_sd,
                                  fem_sd=bb_fem_sd,
                                  mu=cater_mu_bb,
                                  logsigma=cater_ls_bb,
                                  logmax=cater_lm_bb,
                                  theta=th_bb,
                                  logomega=lo_bb,
                                  lWmax=lw_bb,
                                  mu_th=mu_th_bb,
                                  ls_lo=ls_lo_bb,
                                  lm_lo=lm_lo_bb,
                                  lm_lw=lm_lw_bb)

bp_pred <- meanfit_dirsel_av_pred(N=N.h,
                                  link="log",
                                  hd_mean=hd_mean_post, 
                                  hd_sd=hd_sd_post,
                                  omega=omega_bp,
                                  s_sd=bp_s_sd,
                                  y_sd=bp_y_sd,
                                  sy_sd=bp_sy_sd,
                                  fem_sd=bp_fem_sd,
                                  mu=cater_mu_bb,
                                  logsigma=cater_ls_bb,
                                  logmax=cater_lm_bb,
                                  theta=th_bp,
                                  logomega=lo_bp,
                                  lWmax=lw_bp,
                                  mu_th=mu_th_bp,
                                  ls_lo=ls_lo_bp,
                                  lm_lo=lm_lo_bp,
                                  lm_lw=lm_lw_bp)

for (i in 2:4){
  bb_pred[,i] <- round(as.numeric(bb_pred[,i]),4)
  bp_pred[,i] <- round(as.numeric(bp_pred[,i]),4)
}


View(bb_pred)
View(bp_pred)
