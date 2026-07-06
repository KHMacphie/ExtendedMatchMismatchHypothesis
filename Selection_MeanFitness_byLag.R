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
library(lme4)
library(gridExtra)
library(grid)
library(brms)

mod <- readRDS("FullMMH_mod.rds")


## Posteriors
cater_mu_bb <- stanpost(model=mod, parameters=c("c_mu_bb"))
cater_ls_bb <- stanpost(model=mod, parameters=c("c_ls_bb"))
cater_lm_bb <- stanpost(model=mod, parameters=c("c_lm_bb"))

cater_mu_bp <- stanpost(model=mod, parameters=c("c_mu_bp"))
cater_ls_bp <- stanpost(model=mod, parameters=c("c_ls_bp"))
cater_lm_bp <- stanpost(model=mod, parameters=c("c_lm_bp"))

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

bb_s <- stanpost(model=mod, parameters=c("site_effs_bb"))
bb_y <- stanpost(model=mod, parameters=c("year_effs_bb"))
bb_sy <- stanpost(model=mod, parameters=c("siteyear_effs_bb"))

bp_s <- stanpost(model=mod, parameters=c("site_effs_bp"))
bp_y <- stanpost(model=mod, parameters=c("year_effs_bp"))
bp_sy <- stanpost(model=mod, parameters=c("siteyear_effs_bp"))

omega_bb <- stanpost(model=mod, parameters=c("omega_bb"))
omega_bp <- stanpost(model=mod, parameters=c("omega_bp"))

## Predicted fitness by hatch date - requires simulation under model ##
#random draw from random terms inc. covariance for s, y and sy
#into equations for daily estimates
#simulate multiple times from same variance estimates- take mean 
#repeat for every iteration

# f1 = logWmax + mn_lm*lm_lw + random(s+y+sy)
# f2 = date - theta - mn_mu*mu_th - random(s+y+sy)
# f3 = sqrt(2) * exp(logsig + temp*t + random(s+y+sy))
# f4 = sum other ran effs

### Selection gradients ###

# Want selection gradients for different caterpillar peaks
# So use fixed optimum, mean trait and trait variance
# Vary caterpillar peak width and peak height - 0.025 / 0.5 / 0.975 quantiles and other at 0.5

#### Load model with SY mean hatch dates - run in Lag_Prediction
load(file="hd_mod")

# store trait mean and sd
hd_mean_post <- posterior_samples(hd_mod,pars=c("b_Intercept"))
hd_sd_post <- posterior_samples(hd_mod,pars=c("sigma"))

# Organise caterpillar peak values for selection to be estimated at
mu=as.numeric(quantile(colMeans(cater_mu_bb),probs=c(0.1,0.5,0.9)))
logsigma=as.numeric(quantile(colMeans(cater_ls_bb),probs=c(0.1,0.5,0.9)))
logmax=as.numeric(quantile(colMeans(cater_lm_bb),probs=c(0.1,0.5,0.9)))

lag <- c(0,7,14)

selection_cater <- data.frame(para=rep(rep(c("logsigma","logmax"),each=length(logsigma)),length(lag)),
                              mu=mu[2],
                              logsigma=rep(c(logsigma,rep(logsigma[2],length(logmax))),length(lag)),
                              logmax=rep(c(rep(logmax[2],length(logsigma)),logmax),length(lag)),
                              lag=rep(lag,each=(length(c(logsigma,rep(logsigma[2],length(logmax))))),1),
                              bb_d_mn=NA,
                              bb_d_lci=NA,
                              bb_d_uci=NA,
                              bb_d2_mn=NA,
                              bb_d2_lci=NA,
                              bb_d2_uci=NA,
                              bb_pmf_mn=NA,
                              bb_pmf_lci=NA,
                              bb_pmf_uci=NA,
                              bp_d_mn=NA,
                              bp_d_lci=NA,
                              bp_d_uci=NA,
                              bp_d2_mn=NA,
                              bp_d2_lci=NA,
                              bp_d2_uci=NA,
                              bp_pmf_mn=NA,
                              bp_pmf_lci=NA,
                              bp_pmf_uci=NA
)

selection_bb_d <- data.frame(matrix(NA,nrow=nrow(selection_cater),ncol=nrow(th_bb)+5))
colnames(selection_bb_d)[1:5] <- colnames(selection_cater)[1:5]
selection_bb_d[,1:5] <- selection_cater[,1:5]

selection_bp_d <- selection_bb_d
popmeanfit_bb <- selection_bb_d
popmeanfit_bp <- selection_bb_d

N.h <- 10000 
lag_rows <- unname(table(selection_cater$lag)[1])

pb <- txtProgressBar(min = 1, max = nrow(th_bb), style = 3)

# Simulate predictions from the model
for (i in 1:nrow(th_bb)){ # for each iteration
  
  ## draw from random term variances
  
  # one correlation matrix between the 3 parameters, difference variance for each random term
  ## covar matrix: 1=theta, 2=logomega, 3=logit/logWmax
  
  # Bernoulli terms
  cor.matrix.bb <- matrix(c(omega_bb[i,1], omega_bb[i,2], omega_bb[i,3],
                            omega_bb[i,4], omega_bb[i,5], omega_bb[i,6],
                            omega_bb[i,7], omega_bb[i,8], omega_bb[i,9]), nrow=3, ncol=3) #correlation matrix- called omega in model code
  
  site.sd.bb <- diag(c(bb_s_sd[i,1], bb_s_sd[i,2], bb_s_sd[i,3])) # for theta, logomega and logWmax respectively
  year.sd.bb <- diag(c(bb_y_sd[i,1], bb_y_sd[i,2], bb_y_sd[i,3]))
  styr.sd.bb <- diag(c(bb_sy_sd[i,1], bb_sy_sd[i,2], bb_sy_sd[i,3]))
  
  site.covar.bb <- site.sd.bb%*%cor.matrix.bb%*%site.sd.bb  
  year.covar.bb <- year.sd.bb%*%cor.matrix.bb%*%year.sd.bb
  styr.covar.bb <- styr.sd.bb%*%cor.matrix.bb%*%styr.sd.bb
  
  site.ef.bb <- rmvnorm(n=N.h, mean=c(0,0,0), sigma=site.covar.bb) #rmvnorm random number generator for multivariate with VCV
  year.ef.bb <- rmvnorm(n=N.h, mean=c(0,0,0), sigma=year.covar.bb)
  styr.ef.bb <- rmvnorm(n=N.h, mean=c(0,0,0), sigma=styr.covar.bb)
  
  fem.ef.bb <- rnorm(N.h, mean=0, sd=bb_fem_sd[i,1])
  
  # TGPoisson terms
  
  cor.matrix.bp <- matrix(c(omega_bp[i,1], omega_bp[i,2], omega_bp[i,3],
                            omega_bp[i,4], omega_bp[i,5], omega_bp[i,6],
                            omega_bp[i,7], omega_bp[i,8], omega_bp[i,9]), nrow=3, ncol=3) #correlation matrix- called omega in model code
  
  site.sd.bp <- diag(c(bp_s_sd[i,1], bp_s_sd[i,2], bp_s_sd[i,3])) # for theta, logomega and logWmax respectively
  year.sd.bp <- diag(c(bp_y_sd[i,1], bp_y_sd[i,2], bp_y_sd[i,3]))
  styr.sd.bp <- diag(c(bp_sy_sd[i,1], bp_sy_sd[i,2], bp_sy_sd[i,3]))
  
  site.covar.bp <- site.sd.bp%*%cor.matrix.bp%*%site.sd.bp  
  year.covar.bp <- year.sd.bp%*%cor.matrix.bp%*%year.sd.bp
  styr.covar.bp <- styr.sd.bp%*%cor.matrix.bp%*%styr.sd.bp
  
  site.ef.bp <- rmvnorm(n=N.h, mean=c(0,0,0), sigma=site.covar.bp) #rmvnorm random number generator for multivariate with VCV
  year.ef.bp <- rmvnorm(n=N.h, mean=c(0,0,0), sigma=year.covar.bp)
  styr.ef.bp <- rmvnorm(n=N.h, mean=c(0,0,0), sigma=styr.covar.bp)
  
  fem.ef.bp <- rnorm(N.h, mean=0, sd=bp_fem_sd[i,1])
  
  for(j in 1:length(lag)){ # for each level of lag draw N.h hatch dates
    
    theta.bb <- th_bb[i,1] + mu_th_bb[i,1]*mu[2] + site.ef.bb[,1] + year.ef.bb[,1] + styr.ef.bb[,1]
    hd.mean.bb <- lag[j]+theta.bb
    hd.cent.bb <- rnorm(N.h, mean=hd.mean.bb, sd=hd_sd_post$sigma[i])
    hd.cent.cent.bb <- hd.cent.bb-mean(hd.cent.bb)
    
    theta.bp <- th_bp[i,1] + mu_th_bp[i,1]*mu[2] + site.ef.bp[,1] + year.ef.bp[,1] + styr.ef.bp[,1]
    hd.mean.bp <- lag[j]+theta.bp
    hd.cent.bp <- rnorm(N.h, mean=hd.mean.bp, sd=hd_sd_post$sigma[i])
    hd.cent.cent.bp <- hd.cent.bp-mean(hd.cent.bp)
    
    for(h in 1:lag_rows){ # for each combination of theta, logmax and logomega estimate fitness across hatch dates
      
      f1.bb <- lw_bb[i,1] + lm_lw_bb[i,1]*selection_cater$logmax[(h+lag_rows*(j-1))] + site.ef.bb[,3] + year.ef.bb[,3] + styr.ef.bb[,3]
      f2.bb <- hd.cent.bb - th_bb[i,1] - mu_th_bb[i,1]*selection_cater$mu[(h+lag_rows*(j-1))] - site.ef.bb[,1] - year.ef.bb[,1] - styr.ef.bb[,1]
      f3.bb <- sqrt(2) * exp(lo_bb[i,1] + ls_lo_bb[i,1]*selection_cater$logsigma[(h+lag_rows*(j-1))] + lm_lo_bb[i,1]*selection_cater$logmax[(h+lag_rows*(j-1))] + site.ef.bb[,2] + year.ef.bb[,2] + styr.ef.bb[,2])
      f4.bb <- fem.ef.bb
      
      f1.bp <- lw_bp[i,1] + lm_lw_bp[i,1]*selection_cater$logmax[(h+lag_rows*(j-1))] + site.ef.bp[,3] + year.ef.bp[,3] + styr.ef.bp[,3]
      f2.bp <- hd.cent.bp - th_bp[i,1] - mu_th_bp[i,1]*selection_cater$mu[(h+lag_rows*(j-1))] - site.ef.bp[,1] - year.ef.bp[,1] - styr.ef.bp[,1]
      f3.bp <- sqrt(2) * exp(lo_bp[i,1] + ls_lo_bp[i,1]*selection_cater$logsigma[(h+lag_rows*(j-1))] + lm_lo_bp[i,1]*selection_cater$logmax[(h+lag_rows*(j-1))] + site.ef.bp[,2] + year.ef.bp[,2] + styr.ef.bp[,2])
      f4.bp <- fem.ef.bp
      
      EW.bb <- plogis(f1.bb+f4.bb-(f2.bb/f3.bb)^2)
      EW.bp <- exp(f1.bp+f4.bp-(f2.bp/f3.bp)^2)
      
      EW.rel.bb <- EW.bb/mean(EW.bb)
      EW.rel.bp <- EW.bp/mean(EW.bp)
      
      mod.bb <- lm(EW.rel.bb~hd.cent.cent.bb+I(hd.cent.cent.bb^2)) # estimate selection gradient
      mod.bp <- lm(EW.rel.bp~hd.cent.cent.bp+I(hd.cent.cent.bp^2))
      
      selection_bb_d[(h+lag_rows*(j-1)),(i+5)] <- as.numeric(coefficients(mod.bb)[2])
      selection_bb_d2[(h+lag_rows*(j-1)),(i+5)] <- as.numeric(coefficients(mod.bb)[3])
      selection_bp_d[(h+lag_rows*(j-1)),(i+5)] <- as.numeric(coefficients(mod.bp)[2])
      selection_bp_d2[(h+lag_rows*(j-1)),(i+5)] <- as.numeric(coefficients(mod.bp)[3])
      popmeanfit_bb[(h+lag_rows*(j-1)),(i+5)] <- mean(EW.bb) # and mean fitness
      popmeanfit_bp[(h+lag_rows*(j-1)),(i+5)] <- mean(EW.bp)
      
      rm(f1.bb,f2.bb,f3.bb,f4.bb,f1.bp,f2.bp,f3.bp,f4.bp,EW.bb,EW.bp,EW.rel.bb,EW.rel.bp,mod.bb,mod.bp)
    }
    
    rm(hd.cent.bb,hd.cent.cent.bb,theta.bb,hd.mean.bb,hd.cent.bp,hd.cent.cent.bp,theta.bp,hd.mean.bp)
  } 
  
  rm(cor.matrix.bb,site.sd.bb,year.sd.bb,styr.sd.bb,site.covar.bb,year.covar.bb,styr.covar.bb,
     site.ef.bb,year.ef.bb,styr.ef.bb,fem.ef.bb,cor.matrix.bp,site.sd.bp,year.sd.bp,styr.sd.bp,
     site.covar.bp,year.covar.bp,styr.covar.bp,site.ef.bp,year.ef.bp,styr.ef.bp,fem.ef.bp
  )
  
  setTxtProgressBar(pb, i)
} 

close(pb)  

selection_cater$bb_d_mn <- rowMeans(selection_bb_d[,6:ncol(selection_bb_d)])
selection_cater$bb_d_lci <- posterior_interval(t(as.matrix(selection_bb_d[1:nrow(selection_bb_d),6:ncol(selection_bb_d)])), prob=0.95)[,1]
selection_cater$bb_d_uci <- posterior_interval(t(as.matrix(selection_bb_d[1:nrow(selection_bb_d),6:ncol(selection_bb_d)])), prob=0.95)[,2]

selection_cater$bb_d2_mn <- rowMeans(selection_bb_d2[,6:ncol(selection_bb_d2)])
selection_cater$bb_d2_lci <- posterior_interval(t(as.matrix(selection_bb_d2[1:nrow(selection_bb_d2),6:ncol(selection_bb_d2)])), prob=0.95)[,1]
selection_cater$bb_d2_uci <- posterior_interval(t(as.matrix(selection_bb_d2[1:nrow(selection_bb_d2),6:ncol(selection_bb_d2)])), prob=0.95)[,2]

selection_cater$bp_d_mn <- rowMeans(selection_bp_d[,6:ncol(selection_bp_d)])
selection_cater$bp_d_lci <- posterior_interval(t(as.matrix(selection_bp_d[1:nrow(selection_bp_d),6:ncol(selection_bp_d)])), prob=0.95)[,1]
selection_cater$bp_d_uci <- posterior_interval(t(as.matrix(selection_bp_d[1:nrow(selection_bp_d),6:ncol(selection_bp_d)])), prob=0.95)[,2]

selection_cater$bp_d2_mn <- rowMeans(selection_bp_d2[,6:ncol(selection_bp_d2)])
selection_cater$bp_d2_lci <- posterior_interval(t(as.matrix(selection_bp_d2[1:nrow(selection_bp_d2),6:ncol(selection_bp_d2)])), prob=0.95)[,1]
selection_cater$bp_d2_uci <- posterior_interval(t(as.matrix(selection_bp_d2[1:nrow(selection_bp_d2),6:ncol(selection_bp_d2)])), prob=0.95)[,2]

selection_cater$bb_pmf_mn <- rowMeans(popmeanfit_bb[,6:ncol(popmeanfit_bb)])
selection_cater$bb_pmf_lci <- posterior_interval(t(as.matrix(popmeanfit_bb[1:nrow(popmeanfit_bb),6:ncol(popmeanfit_bb)])), prob=0.95)[,1]
selection_cater$bb_pmf_uci <- posterior_interval(t(as.matrix(popmeanfit_bb[1:nrow(popmeanfit_bb),6:ncol(popmeanfit_bb)])), prob=0.95)[,2]

selection_cater$bp_pmf_mn <- rowMeans(popmeanfit_bp[,6:ncol(popmeanfit_bp)])
selection_cater$bp_pmf_lci <- posterior_interval(t(as.matrix(popmeanfit_bp[1:nrow(popmeanfit_bp),6:ncol(popmeanfit_bp)])), prob=0.95)[,1]
selection_cater$bp_pmf_uci <- posterior_interval(t(as.matrix(popmeanfit_bp[1:nrow(popmeanfit_bp),6:ncol(popmeanfit_bp)])), prob=0.95)[,2]

write.csv(selection_cater,"Selection&PopFitEstimates.csv",row.names = F)
#selection_cater <- read.csv("Selection&PopFitEstimates.csv")

# Difference in popfit and selection by cater peak
cater_fitsel_dif <-data.frame(matrix(NA,nrow=24,ncol=7))
colnames(cater_fitsel_dif) <-  c("mod","para","coef","lag","mean","lci","uci") 

d_bb_logsigma_0 <- subset(selection_bb_d,para=="logsigma"&lag==0)
d_bb_logsigma_7 <- subset(selection_bb_d,para=="logsigma"&lag==7)
d_bb_logsigma_14 <- subset(selection_bb_d,para=="logsigma"&lag==14)

pf_bb_logsigma_0 <- subset(popmeanfit_bb,para=="logsigma"&lag==0)
pf_bb_logsigma_7 <- subset(popmeanfit_bb,para=="logsigma"&lag==7)
pf_bb_logsigma_14 <- subset(popmeanfit_bb,para=="logsigma"&lag==14)

d_bp_logsigma_0 <- subset(selection_bp_d,para=="logsigma"&lag==0)
d_bp_logsigma_7 <- subset(selection_bp_d,para=="logsigma"&lag==7)
d_bp_logsigma_14 <- subset(selection_bp_d,para=="logsigma"&lag==14)

pf_bp_logsigma_0 <- subset(popmeanfit_bp,para=="logsigma"&lag==0)
pf_bp_logsigma_7 <- subset(popmeanfit_bp,para=="logsigma"&lag==7)
pf_bp_logsigma_14 <- subset(popmeanfit_bp,para=="logsigma"&lag==14)

d_bb_logmax_0 <- subset(selection_bb_d,para=="logmax"&lag==0)
d_bb_logmax_7 <- subset(selection_bb_d,para=="logmax"&lag==7)
d_bb_logmax_14 <- subset(selection_bb_d,para=="logmax"&lag==14)

pf_bb_logmax_0 <- subset(popmeanfit_bb,para=="logmax"&lag==0)
pf_bb_logmax_7 <- subset(popmeanfit_bb,para=="logmax"&lag==7)
pf_bb_logmax_14 <- subset(popmeanfit_bb,para=="logmax"&lag==14)

d_bp_logmax_0 <- subset(selection_bp_d,para=="logmax"&lag==0)
d_bp_logmax_7 <- subset(selection_bp_d,para=="logmax"&lag==7)
d_bp_logmax_14 <- subset(selection_bp_d,para=="logmax"&lag==14)

pf_bp_logmax_0 <- subset(popmeanfit_bp,para=="logmax"&lag==0)
pf_bp_logmax_7 <- subset(popmeanfit_bp,para=="logmax"&lag==7)
pf_bp_logmax_14 <- subset(popmeanfit_bp,para=="logmax"&lag==14)

X <- d_bb_logsigma_0
cater_fitsel_dif[1,] <- c("bb","logsig","dirsel","0",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bp_logsigma_0
cater_fitsel_dif[2,] <- c("bp","logsig","dirsel","0",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bb_logsigma_7
cater_fitsel_dif[3,] <- c("bb","logsig","dirsel","7",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bp_logsigma_7
cater_fitsel_dif[4,] <- c("bp","logsig","dirsel","7",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bb_logsigma_14
cater_fitsel_dif[5,] <- c("bb","logsig","dirsel","14",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bp_logsigma_14
cater_fitsel_dif[6,] <- c("bp","logsig","dirsel","14",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bb_logmax_0
cater_fitsel_dif[7,] <- c("bb","logmax","dirsel","0",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bp_logmax_0
cater_fitsel_dif[8,] <- c("bp","logmax","dirsel","0",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bb_logmax_7
cater_fitsel_dif[9,] <- c("bb","logmax","dirsel","7",
                          mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                          posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bp_logmax_7
cater_fitsel_dif[10,] <- c("bp","logmax","dirsel","7",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bb_logmax_14
cater_fitsel_dif[11,] <- c("bb","logmax","dirsel","14",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- d_bp_logmax_14
cater_fitsel_dif[12,] <- c("bp","logmax","dirsel","14",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bb_logsigma_0
cater_fitsel_dif[13,] <- c("bb","logsig","popfit","0",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bp_logsigma_0
cater_fitsel_dif[14,] <- c("bp","logsig","popfit","0",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bb_logsigma_7
cater_fitsel_dif[15,] <- c("bb","logsig","popfit","7",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bp_logsigma_7
cater_fitsel_dif[16,] <- c("bp","logsig","popfit","7",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bb_logsigma_14
cater_fitsel_dif[17,] <- c("bb","logsig","popfit","14",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bp_logsigma_14
cater_fitsel_dif[18,] <- c("bp","logsig","popfit","14",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bb_logmax_0
cater_fitsel_dif[19,] <- c("bb","logmax","popfit","0",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bp_logmax_0
cater_fitsel_dif[20,] <- c("bp","logmax","popfit","0",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bb_logmax_7
cater_fitsel_dif[21,] <- c("bb","logmax","popfit","7",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bp_logmax_7
cater_fitsel_dif[22,] <- c("bp","logmax","popfit","7",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bb_logmax_14
cater_fitsel_dif[23,] <- c("bb","logmax","popfit","14",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

X <- pf_bp_logmax_14
cater_fitsel_dif[24,] <- c("bp","logmax","popfit","14",
                           mean(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])),
                           posterior_interval(as.matrix(as.numeric(X[1,6:ncol(X)])-as.numeric(X[3,6:ncol(X)])), prob=0.95))

cater_fitsel_dif$mean <- round(as.numeric(cater_fitsel_dif$mean),4)
cater_fitsel_dif$lci <- round(as.numeric(cater_fitsel_dif$lci),4)
cater_fitsel_dif$uci <- round(as.numeric(cater_fitsel_dif$uci),4)

write.csv(cater_fitsel_dif,"Selection&PopFitDifByCater.csv",row.names = F)


## plot selection 

selection_logsigma <- subset(selection_cater,para=="logsigma")
selection_logmax <- subset(selection_cater,para=="logmax")

selection_logsigma$logsigmaQ <- c(0.1,0.5,0.9)
selection_logmax$logmaxQ <- c(0.1,0.5,0.9)

selection_logsigma$lag <- paste("Lag =",selection_logsigma$lag)
selection_logsigma$lag <- factor(selection_logsigma$lag, levels=c("Lag = 0","Lag = 7","Lag = 14"))

selection_logmax$lag <- paste("Lag =",selection_logmax$lag)
selection_logmax$lag <- factor(selection_logmax$lag, levels=c("Lag = 0","Lag = 7","Lag = 14"))

lm_cols <- c("#A6CFF5",  "#4289CC",  "#06529E")
ls_cols <- c("#A2E0AA",  "#32AB46", "#046E18")

(plot_bb_pf_ls <- ggplot(selection_logsigma,aes(as.factor(logsigmaQ),bb_pmf_mn,col=as.factor(logsigmaQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logsigmaQ),ymin=bb_pmf_lci,ymax=bb_pmf_uci),width=0.5)+
    guides(color = "none")+
    scale_colour_manual(values=ls_cols)+
    scale_y_continuous(limits=c(0,1),
                       breaks=c(0,0.25,0.5,0.75,1),
                       labels=c(" 0.00"," 0.25"," 0.50"," 0.75"," 1.00"))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylab("")+
    xlab(""))

(plot_bp_pf_ls <- ggplot(selection_logsigma,aes(as.factor(logsigmaQ),bp_pmf_mn,col=as.factor(logsigmaQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logsigmaQ),ymin=bp_pmf_lci,ymax=bp_pmf_uci),width=0.5)+
    guides(color = "none")+
    scale_colour_manual(values=ls_cols)+
    scale_y_continuous(limits=c(3,10),
                       breaks=c(4,6,8,10),
                       labels=c("     4","     6","     8","    10"))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylab("")+
    xlab(""))

(plot_bb_d_ls <- ggplot(selection_logsigma,aes(as.factor(logsigmaQ),bb_d_mn,col=as.factor(logsigmaQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logsigmaQ),ymin=bb_d_lci,ymax=bb_d_uci),width=0.5)+
    guides(color = "none")+
    scale_colour_manual(values=ls_cols)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          strip.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylim(-0.12,0.01)+
    ylab("")+
    xlab(""))

(plot_bp_d_ls <- ggplot(selection_logsigma,aes(as.factor(logsigmaQ),bp_d_mn,col=as.factor(logsigmaQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logsigmaQ),ymin=bp_d_lci,ymax=bp_d_uci),width=0.5)+
    guides(color = "none")+
    scale_colour_manual(values=ls_cols)+
    #scale_y_continuous(breaks=c(-0.01,-0.02,-0.03),
    #                   labels=c("  -0.01","  -0.02","  -0.03"))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          strip.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylim(-0.12,0.01)+
    ylab("")+
    xlab(""))



(plot_bb_pf_lm <- ggplot(selection_logmax,aes(as.factor(logmaxQ),bb_pmf_mn,col=as.factor(logmaxQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logmaxQ),ymin=bb_pmf_lci,ymax=bb_pmf_uci),width=0.5)+
    guides(color = "none")+
    annotate("segment", x = 1, xend = 3, y = 0.03, yend = 0.03,col="grey")+
    annotate("point", x = 2, y = 0.03,size=3,col="white")+
    annotate("text", label="*",x = 2, y = 0.012,size=5,col="grey")+
    scale_colour_manual(values=lm_cols)+
    scale_y_continuous(limits=c(0,1),
                       breaks=c(0,0.25,0.5,0.75,1),
                       labels=c(" 0.00"," 0.25"," 0.50"," 0.75"," 1.00"))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylab("")+
    xlab(""))

(plot_bp_pf_lm <- ggplot(selection_logmax,aes(as.factor(logmaxQ),bp_pmf_mn,col=as.factor(logmaxQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logmaxQ),ymin=bp_pmf_lci,ymax=bp_pmf_uci),width=0.5)+
    guides(color = "none")+
    annotate("segment", x = 1, xend = 3, y = 3.2, yend = 3.2,col="grey")+
    annotate("point", x = 2, y = 3.2,size=3,col="white")+
    annotate("text", label="*",x = 2, y = 3.05,size=5,col="grey")+
    scale_colour_manual(values=lm_cols)+
    scale_y_continuous(limits=c(3,10),
                       breaks=c(4,6,8,10),
                       labels=c("     4","     6","     8","    10"))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylab("")+
    xlab(""))

(plot_bb_d_lm <- ggplot(selection_logmax,aes(as.factor(logmaxQ),bb_d_mn,col=as.factor(logmaxQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logmaxQ),ymin=bb_d_lci,ymax=bb_d_uci),width=0.5)+
    guides(color = "none")+
    scale_colour_manual(values=lm_cols)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          strip.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylim(-0.12,0.01)+
    ylab("")+
    xlab(""))

(plot_bp_d_lm <- ggplot(selection_logmax,aes(as.factor(logmaxQ),bp_d_mn,col=as.factor(logmaxQ)))+
    geom_point()+
    geom_errorbar(aes(x=as.factor(logmaxQ),ymin=bp_d_lci,ymax=bp_d_uci),width=0.5)+
    guides(color = "none")+
    scale_colour_manual(values=lm_cols)+
    #scale_y_continuous(breaks=c(-0.01,-0.02,-0.03),
    #                   labels=c("  -0.01","  -0.02","  -0.03"))+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          strip.text.x = element_blank(),
          panel.border = element_rect(linewidth=1.2))+
    facet_wrap(~lag,nrow=1)+
    ylim(-0.12,0.01)+
    ylab("")+
    xlab(""))

space <- ggplot()+theme_void()

width_plot <- grid.arrange(space,space,space,space,plot_bb_pf_ls,plot_bp_pf_ls,space,plot_bb_d_ls,plot_bp_d_ls,space,ncol=3,nrow=4,heights=c(0.22,1,1,0.15),widths=c(0.08,1,1))
height_plot <- grid.arrange(space,space,space,space,plot_bb_pf_lm,plot_bp_pf_lm,space,plot_bb_d_lm,plot_bp_d_lm,space,ncol=3,nrow=4,heights=c(0.22,1,1,0.15),widths=c(0.08,1,1))
grid.arrange(height_plot,space,width_plot,nrow=3,heights=c(1,0.05,1)) #6.5*9"
grid.text(label=expression(bold("Probability of success")),x=(0.33), y=(0.97), rot=0)
grid.text(label=expression(bold("Number fledged")),x=(0.79), y=(0.97), rot=0)
grid.text(label="Direc. Selection",x=(0.04), y=(0.66), rot=90)
grid.text(label="Pop. mean fitness",x=(0.04), y=(0.85), rot=90)
grid.text(label="Caterpillar peak height (quantile)",x=(0.55), y=(0.53), rot=0) #5*6

grid.text(label=expression(bold("Probability of success")),x=(0.33), y=(0.46), rot=0)
grid.text(label=expression(bold("Number fledged")),x=(0.79), y=(0.46), rot=0)
grid.text(label="Direc. Selection",x=(0.04), y=(0.15), rot=90)
grid.text(label="Pop. mean fitness",x=(0.04), y=(0.33), rot=90)
grid.text(label="Caterpillar peak width (quantile)",x=(0.57), y=(0.02), rot=0) #5*6


### Differences across lag by caterpilalr metric

lagdif <- function(dat,para){
  dat$logsigma <- round(dat$logsigma,1)
  dat$logmax <- round(dat$logmax,1)
  
  dat$logsigma <- ifelse(dat$logsigma==max(dat$logsigma),"q0.9",
                         ifelse(dat$logsigma==min(dat$logsigma),"q0.1","q0.5"))
  dat$logmax <- ifelse(dat$logmax==max(dat$logmax),"q0.9",
                       ifelse(dat$logmax==min(dat$logmax),"q0.1","q0.5"))
  term <- para
  dat_p <- subset(dat,para==term)
  dat_lg <- gather(dat_p,key="iteration",value="estimate",6:ncol(dat_p))
  dat_sp <- spread(dat_lg,key=term,value="estimate")
  dat_0 <- subset(dat_sp, lag==0)
  dat_7 <- subset(dat_sp, lag==7)
  dat_14 <- subset(dat_sp, lag==14)
  
  h_0.9 <- dat_0$q0.9-dat_14$q0.9
  h_0.1 <- dat_0$q0.1-dat_14$q0.1
  
  mn <- mean(as.numeric(h_0.9-h_0.1))
  cis <- posterior_interval(as.matrix(as.numeric(h_0.9-h_0.1), prob=0.95))
  
  return(c(mn,cis))
}

lagdif_df <- data.frame(model=c(rep(c("bb","bb","bp","bp"),2)),
                        coef=c(rep("meanfit",4),rep("dirsel",4)),
                        para=c(rep(c("logmax","logsigma"),4)),
                        mn=NA,
                        lci=NA,
                        uci=NA) 

lagdif_df[1,c(4:6)] <- lagdif(dat=popmeanfit_bb,para="logmax")
lagdif_df[2,c(4:6)] <- lagdif(dat=popmeanfit_bb,para="logsigma")
lagdif_df[3,c(4:6)] <- lagdif(dat=popmeanfit_bp,para="logmax")
lagdif_df[4,c(4:6)] <- lagdif(dat=popmeanfit_bp,para="logsigma")

lagdif_df[5,c(4:6)] <- lagdif(dat=selection_bb_d,para="logmax")
lagdif_df[6,c(4:6)] <- lagdif(dat=selection_bb_d,para="logsigma")
lagdif_df[7,c(4:6)] <- lagdif(dat=selection_bp_d,para="logmax")
lagdif_df[8,c(4:6)] <- lagdif(dat=selection_bp_d,para="logsigma")

write.csv(lagdif_df, "LagEffDifByCater.csv", row.names = F)


lageffect <- function(dat,para){
  dat$logsigma <- round(dat$logsigma,1)
  dat$logmax <- round(dat$logmax,1)
  
  dat$logsigma <- ifelse(dat$logsigma==max(dat$logsigma),"q0.9",
                         ifelse(dat$logsigma==min(dat$logsigma),"q0.1","q0.5"))
  dat$logmax <- ifelse(dat$logmax==max(dat$logmax),"q0.9",
                       ifelse(dat$logmax==min(dat$logmax),"q0.1","q0.5"))
  term <- para
  dat_p <- subset(dat,para==term)
  dat_lg <- gather(dat_p,key="iteration",value="estimate",6:ncol(dat_p))
  dat_sp <- spread(dat_lg,key=term,value="estimate")
  dat_0 <- subset(dat_sp, lag==0)
  dat_7 <- subset(dat_sp, lag==7)
  dat_14 <- subset(dat_sp, lag==14)
  
  h_0.9 <- dat_0$q0.9-dat_14$q0.9
  h_0.5 <- dat_0$q0.5-dat_14$q0.5
  h_0.1 <- dat_0$q0.1-dat_14$q0.1
  
  lagef <- data.frame(c=c(0.9,0.5,0.1),
                      mean=c(mean(h_0.9),mean(h_0.5),mean(h_0.1)),
                      lci=c(posterior_interval(as.matrix(h_0.9),prob=0.95)[1],posterior_interval(as.matrix(h_0.5),prob=0.95)[1],posterior_interval(as.matrix(h_0.1),prob=0.95)[1]),
                      uci=c(posterior_interval(as.matrix(h_0.9),prob=0.95)[2],posterior_interval(as.matrix(h_0.5),prob=0.95)[2],posterior_interval(as.matrix(h_0.1),prob=0.95)[2]))
  
  
  return(lagef)
}

lageffect(dat=popmeanfit_bb,para="logmax")
lageffect(dat=popmeanfit_bb,para="logsigma")
lageffect(dat=popmeanfit_bp,para="logmax")
lageffect(dat=popmeanfit_bp,para="logsigma")

lageffect(dat=selection_bb_d,para="logmax")
lageffect(dat=selection_bb_d,para="logsigma")
lageffect(dat=selection_bp_d,para="logmax")
lageffect(dat=selection_bp_d,para="logsigma")

