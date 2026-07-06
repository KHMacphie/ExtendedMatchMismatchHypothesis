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
rstan::get_num_divergent(mod)


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


## Predicted fitness by hatch date ##
#random draw from random terms inc. covariance for s, y and sy
#into equations for daily estimates
#simulate multiple times from same variance estimates- take mean 
#repeat for every iteration

# f1 = logWmax + mn_lm*lm_lw + random(s+y+sy)
# f2 = date - theta - mn_mu*mu_th - random(s+y+sy)
# f3 = sqrt(2) * exp(logsig + temp*t + random(s+y+sy))
# f4 = sum other ran effs

# Values to calculate from within sim
dat.cnt <- seq(min(SYBhd$hd_cent),max(SYBhd$hd_cent),0.5)
N.h <- 1000 

# Set up dfs for simulation
mu=as.numeric(quantile(colMeans(cater_mu_bb),probs=c(0.025,0.25,0.5,0.75,0.975)))
logsigma=as.numeric(quantile(colMeans(cater_ls_bb),probs=c(0.025,0.25,0.5,0.75,0.975)))
logmax=as.numeric(quantile(colMeans(cater_lm_bb),probs=c(0.025,0.25,0.5,0.75,0.975)))

pred_range <- data.frame(para=rep(c("mu","logsigma","logmax"),each=length(dat.cnt)*length(logsigma)),
                         dat.cnt=rep(dat.cnt,length(mu)+length(logsigma)+length(logmax)),
                         mu=c(rep(mu,each=length(dat.cnt),1),rep(mu[3],length(dat.cnt)*(length(logsigma)+length(logmax)))),
                         logsigma=c(rep(logsigma[3],length(dat.cnt)*(length(mu))),rep(logsigma,each=length(dat.cnt),1),rep(logsigma[3],length(dat.cnt)*(length(logmax)))),
                         logmax=c(rep(logmax[3],length(dat.cnt)*(length(mu)+length(logsigma))),rep(logmax,each=length(dat.cnt),1)))

pred_range$sigma <- exp(pred_range$logsigma)
pred_range$max <- exp(pred_range$logmax)

gaussianfunc_pred <- function(N,x_range,link,omega,s_sd,y_sd,sy_sd,fem_sd,theta,logomega,lWmax,mu_th,ls_lo,lm_lo,lm_lw){
  
  Esim <- data.frame(matrix(NA, nrow = nrow(x_range), ncol = N+ncol(x_range))) #for mean expectations
  Emean <- data.frame(matrix(NA, nrow = nrow(x_range), ncol = nrow(theta)+ncol(x_range))) #for mean expectations
  
  #### Need to bring in x_range columns
  colnames(Esim)[1:ncol(x_range)] <- colnames(x_range)
  colnames(Emean)[1:ncol(x_range)] <- colnames(x_range)
  
  Esim[,1:ncol(x_range)] <- x_range
  Emean[,1:ncol(x_range)] <- x_range
  
  pb <- txtProgressBar(min = 1, max = nrow(theta), style = 3)
  
  # Simulate predictions from the model
  for (i in 1:nrow(theta)){ # for each iteration
    #one correlation matrix between the 3 parameters, difference variance for each random term
    ##covar matrix: 1=theta, 2=logomega, 3=logit/logWmax
    cor.matrix <- matrix(c(omega[i,1], omega[i,2], omega[i,3],
                           omega[i,4], omega[i,5], omega[i,6],
                           omega[i,7], omega[i,8], omega[i,9]), nrow=3, ncol=3) #correlation matrix- called omega in model code
    
    site.sd <- diag(c(s_sd[i,1], s_sd[i,2], s_sd[i,3])) # for theta, logomega and logWmax respectively
    year.sd <- diag(c(y_sd[i,1], y_sd[i,2], y_sd[i,3]))
    styr.sd <- diag(c(sy_sd[i,1], sy_sd[i,2], sy_sd[i,3]))
    
    site.covar <- site.sd%*%cor.matrix%*%site.sd  
    year.covar <- year.sd%*%cor.matrix%*%year.sd
    styr.covar <- styr.sd%*%cor.matrix%*%styr.sd
    
    for(h in 1:N){ 
      site.ef <- rmvnorm(n=1, mean=c(0,0,0), sigma=site.covar) #rmvnorm random number generator for multivariate with VCV
      year.ef <- rmvnorm(n=1, mean=c(0,0,0), sigma=year.covar)
      styr.ef <- rmvnorm(n=1, mean=c(0,0,0), sigma=styr.covar)
      
      fem.ef <- rnorm(1, mean=0, sd=fem_sd[i,1])
      
      f1 <- lWmax[i,1] + lm_lw[i,1]*Esim$logmax + site.ef[3] + year.ef[3] + styr.ef[3]
      f2 <- Esim$dat.cnt - theta[i,1] - mu_th[i,1]*Esim$mu - site.ef[1] - year.ef[1] - styr.ef[1]
      f3 <- sqrt(2) * exp(logomega[i,1] + ls_lo[i,1]*Esim$logsigma + lm_lo[i,1]*Esim$logmax + site.ef[2] + year.ef[2] + styr.ef[2])
      
      f4 <- fem.ef
      
      if(link=="logit"){
        Esim[,h+ncol(x_range)] <- plogis(f1+f4-(f2/f3)^2)
      } else if(link=="log"){
        Esim[,h+ncol(x_range)] <- exp(f1+f4-(f2/f3)^2)
      } 
      
      f1 <- NULL
      f2 <- NULL
      f3 <- NULL
      f4 <- NULL
      
    } 
    
    Emean[,i+ncol(x_range)] <- apply(Esim[,(ncol(x_range)+1):ncol(Esim)], 1, mean)
    
    Esim[,c((ncol(x_range)+1):ncol(Esim))] <- NA 
    setTxtProgressBar(pb, i)
  } 
  close(pb)  
  return(Emean)
}


## Bernoulli part prediction

bb_pred <- gaussianfunc_pred(N=N.h,
                             x_range=pred_range,
                             link="logit",
                             omega=omega_bb,
                             s_sd=bb_s_sd,
                             y_sd=bb_y_sd,
                             sy_sd=bb_sy_sd,
                             fem_sd=bb_fem_sd,
                             theta=th_bb,
                             logomega=lo_bb,
                             lWmax=lw_bb,
                             mu_th=mu_th_bb,
                             ls_lo=ls_lo_bb,
                             lm_lo=lm_lo_bb,
                             lm_lw=lm_lw_bb)

# Mean and CIs for estimate at each date or temperature
bb_pred$mean <- apply(bb_pred[,(ncol(pred_range)+1):ncol(bb_pred)], 1, mean)
bb_pred$lwci <- posterior_interval(t(as.matrix(bb_pred[1:nrow(bb_pred),(ncol(pred_range)+1):(ncol(bb_pred)-1)])), prob=0.95)[,1]
bb_pred$upci <- posterior_interval(t(as.matrix(bb_pred[1:nrow(bb_pred),(ncol(pred_range)+1):(ncol(bb_pred)-2)])), prob=0.95)[,2]
bb_pred$dat <- bb_pred$dat.cnt+146

#write.csv(bb_pred,"bb_pred_function.csv",row.names = F)



## Poisson part prediction

bp_pred <- gaussianfunc_pred(N=N.h,
                             x_range=pred_range,
                             link="log",
                             omega=omega_bp,
                             s_sd=bp_s_sd,
                             y_sd=bp_y_sd,
                             sy_sd=bp_sy_sd,
                             fem_sd=bp_fem_sd,
                             theta=th_bp,
                             logomega=lo_bp,
                             lWmax=lw_bp,
                             mu_th=mu_th_bp,
                             ls_lo=ls_lo_bp,
                             lm_lo=lm_lo_bp,
                             lm_lw=lm_lw_bp)

# Mean and CIs for estimate at each date or temperature
bp_pred$mean <- apply(bp_pred[,(ncol(pred_range)+1):ncol(bp_pred)], 1, mean)
bp_pred$lwci <- posterior_interval(t(as.matrix(bp_pred[1:nrow(bp_pred),(ncol(pred_range)+1):(ncol(bp_pred)-1)])), prob=0.95)[,1]
bp_pred$upci <- posterior_interval(t(as.matrix(bp_pred[1:nrow(bp_pred),(ncol(pred_range)+1):(ncol(bp_pred)-2)])), prob=0.95)[,2]
bp_pred$dat <- bp_pred$dat.cnt+146

#write.csv(bp_pred,"bp_pred_function.csv",row.names = F)