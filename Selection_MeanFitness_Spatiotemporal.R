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
library(brms)
library(gridExtra)
library(grid)
library(cowplot)
library(viridis)
library(ggh4x)
library(grid)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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


## Hatch date model - run in Lag_Prediction
load(file="hd_mod")


### Site year mean hatch dates
hd_ef <- posterior_samples(hd_mod,pars=c("year","site"))
s_hd_ef <- hd_ef[,c(384:427)]
y_hd_ef <- hd_ef[,c(372:382)]
sy_hd_ef <- hd_ef[,c(3:371)]
int_post <- posterior_samples(hd_mod,pars=c("b_Intercept"))

sy_hd_post <- data.frame(matrix(NA,ncol=ncol(sy_hd_ef),nrow=nrow(sy_hd_ef))) 

for(i in 1:ncol(sy_hd_post)){
  sy_hd_post[,i] <- int_post[,1]+sy_hd_ef[,i]+y_hd_ef[,SY_id$year_id[i]]+s_hd_ef[,SY_id$site_id[i]]
}

### SD in hatch dates within site year
hd_sd_post <- posterior_samples(hd_mod,pars=c("sigma"))

### Site year optimum hatch dates
opt_bb <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy)))
opt_bp <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy)))

for(i in 1:nrow(SY_id)){
  opt_bb[,i] <- th_bb[,1]+mu_th_bb[,1]*cater_mu_bb[,SY_id$siteyear_id[i]]+bb_s[,SY_id$site_id[i]]+bb_y[,SY_id$year_id[i]]+bb_sy[,SY_id$siteyear_id[i]]
}
for(i in 1:nrow(SY_id_p)){
  opt_bp[,i] <- th_bp[,1]+mu_th_bp[,1]*cater_mu_bp[,SY_id_p$siteyear_id[i]]+bp_s[,SY_id_p$site_id[i]]+bp_y[,SY_id_p$year_id[i]]+bp_sy[,SY_id_p$siteyear_id[i]]
}

#Just caterpillar
opt_bb_c <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy)))
opt_bp_c <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy)))

for(i in 1:nrow(SY_id)){
  opt_bb_c[,i] <- th_bb[,1]+mu_th_bb[,1]*cater_mu_bb[,SY_id$siteyear_id[i]]
}
for(i in 1:nrow(SY_id_p)){
  opt_bp_c[,i] <- th_bp[,1]+mu_th_bp[,1]*cater_mu_bp[,SY_id_p$siteyear_id[i]]
}

### Site year fitness function width
sy_lo_bb <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy)))
sy_lo_bp <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy)))

for(i in 1:nrow(SY_id)){
  sy_lo_bb[,i] <- lo_bb[,1]+ls_lo_bb[,1]*cater_ls_bb[,SY_id$siteyear_id[i]]+lm_lo_bb[,1]*cater_lm_bb[,SY_id$siteyear_id[i]]+bb_s[,SY_id$site_id[i]+44]+bb_y[,SY_id$year_id[i]+11]+bb_sy[,SY_id$siteyear_id[i]+369]
}
for(i in 1:nrow(SY_id_p)){
  sy_lo_bp[,i] <- lo_bp[,1]+ls_lo_bp[,1]*cater_ls_bp[,SY_id_p$siteyear_id[i]]+lm_lo_bp[,1]*cater_lm_bp[,SY_id_p$siteyear_id[i]]+bp_s[,SY_id_p$site_id[i]+44]+bp_y[,SY_id_p$year_id[i]+11]+bp_sy[,SY_id_p$siteyear_id[i]+359]
}

# Just caterpillar
sy_lo_bb_c <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy)))
sy_lo_bp_c <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy)))

for(i in 1:nrow(SY_id)){
  sy_lo_bb_c[,i] <- lo_bb[,1]+ls_lo_bb[,1]*cater_ls_bb[,SY_id$siteyear_id[i]]+lm_lo_bb[,1]*cater_lm_bb[,SY_id$siteyear_id[i]]
}
for(i in 1:nrow(SY_id_p)){
  sy_lo_bp_c[,i] <- lo_bp[,1]+ls_lo_bp[,1]*cater_ls_bp[,SY_id_p$siteyear_id[i]]+lm_lo_bp[,1]*cater_lm_bp[,SY_id_p$siteyear_id[i]]
}

### Site year fitness function height
sy_lw_bb <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy)))
sy_lw_bp <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy)))

for(i in 1:nrow(SY_id)){
  sy_lw_bb[,i] <- lw_bb[,1]+lm_lw_bb[,1]*cater_lm_bb[,SY_id$siteyear_id[i]]+bb_s[,SY_id$site_id[i]+44+44]+bb_y[,SY_id$year_id[i]+11+11]+bb_sy[,SY_id$siteyear_id[i]+369+369]
}
for(i in 1:nrow(SY_id_p)){
  sy_lw_bp[,i] <- lw_bp[,1]+lm_lw_bp[,1]*cater_lm_bp[,SY_id_p$siteyear_id[i]]+bp_s[,SY_id_p$site_id[i]+44+44]+bp_y[,SY_id_p$year_id[i]+11+11]+bp_sy[,SY_id_p$siteyear_id[i]+359+359]
}

#Just caterpillars
sy_lw_bb_c <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy)))
sy_lw_bp_c <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy)))

for(i in 1:nrow(SY_id)){
  sy_lw_bb_c[,i] <- lw_bb[,1]+lm_lw_bb[,1]*cater_lm_bb[,SY_id$siteyear_id[i]]
}
for(i in 1:nrow(SY_id_p)){
  sy_lw_bp_c[,i] <- lw_bp[,1]+lm_lw_bp[,1]*cater_lm_bp[,SY_id_p$siteyear_id[i]]
}

# for each iteration
# draw 1000 hds from sy mean and hd_sd
# fitness from SY specific function
# calculate selection and mean fitness

sy_d_bb <- data.frame(matrix(NA,nrow=nrow(sy_lw_bb),ncol=ncol(sy_lw_bb)))
sy_d_bb_c <- data.frame(matrix(NA,nrow=nrow(sy_lw_bb),ncol=ncol(sy_lw_bb)))
sy_d2_bb <- data.frame(matrix(NA,nrow=nrow(sy_lw_bb),ncol=ncol(sy_lw_bb)))
sy_d2_bb_c <- data.frame(matrix(NA,nrow=nrow(sy_lw_bb),ncol=ncol(sy_lw_bb)))
sy_mf_bb <- data.frame(matrix(NA,nrow=nrow(sy_lw_bb),ncol=ncol(sy_lw_bb)))
sy_mf_bb_c <- data.frame(matrix(NA,nrow=nrow(sy_lw_bb),ncol=ncol(sy_lw_bb)))

sy_d_bp <- data.frame(matrix(NA,nrow=nrow(sy_lw_bp),ncol=ncol(sy_lw_bp)))
sy_d_bp_c <- data.frame(matrix(NA,nrow=nrow(sy_lw_bp),ncol=ncol(sy_lw_bp)))
sy_d2_bp <- data.frame(matrix(NA,nrow=nrow(sy_lw_bp),ncol=ncol(sy_lw_bp)))
sy_d2_bp_c <- data.frame(matrix(NA,nrow=nrow(sy_lw_bp),ncol=ncol(sy_lw_bp)))
sy_mf_bp <- data.frame(matrix(NA,nrow=nrow(sy_lw_bp),ncol=ncol(sy_lw_bp)))
sy_mf_bp_c <- data.frame(matrix(NA,nrow=nrow(sy_lw_bp),ncol=ncol(sy_lw_bp)))

n.hd <- 1000

pb <- txtProgressBar(min = 1, max = nrow(sy_hd_post), style = 3)

for(i in 1:nrow(sy_hd_post)){
  
  EW.bb <- data.frame(matrix(NA,nrow=n.hd,ncol=nrow(SY_id)))
  EW.bb.c <- data.frame(matrix(NA,nrow=n.hd,ncol=nrow(SY_id)))
  EW.bp <- data.frame(matrix(NA,nrow=n.hd,ncol=nrow(SY_id_p)))
  EW.bp.c <- data.frame(matrix(NA,nrow=n.hd,ncol=nrow(SY_id_p)))
  
  sim.hd.bb <- data.frame(matrix(NA,nrow=n.hd,ncol=ncol(sy_hd_post)))
  
  for(j in 1:n.hd){
    sim.hd.bb[j,] <- rnorm(ncol(sy_hd_post), mean=as.numeric(sy_hd_post[i,]), sd=hd_sd_post[i,1])
    sim.hd.bp <- sim.hd.bb[,SY_id_p$siteyear_id_c]
    fem.ef.bb <- rnorm(ncol(sy_hd_post), mean=0, sd=bb_fem_sd[i,1])
    fem.ef.bp <- rnorm(ncol(sim.hd.bp), mean=0, sd=bp_fem_sd[i,1])
    
    EW.bb[j,] <- plogis((fem.ef.bb+as.numeric(sy_lw_bb[i,]))-((as.numeric(sim.hd.bb[j,])-as.numeric(opt_bb[i,]))^2)/(2*exp(as.numeric(sy_lo_bb[i,]))^2))
    EW.bb.c[j,] <- plogis((fem.ef.bb+as.numeric(sy_lw_bb_c[i,]))-((as.numeric(sim.hd.bb[j,])-as.numeric(opt_bb_c[i,]))^2)/(2*exp(as.numeric(sy_lo_bb_c[i,]))^2))
    
    EW.bp[j,] <- exp((fem.ef.bp+as.numeric(sy_lw_bp[i,]))-((as.numeric(sim.hd.bp[j,])-as.numeric(opt_bp[i,]))^2)/(2*exp(as.numeric(sy_lo_bp[i,]))^2))
    EW.bp.c[j,] <- exp((fem.ef.bp+as.numeric(sy_lw_bp_c[i,]))-((as.numeric(sim.hd.bp[j,])-as.numeric(opt_bp_c[i,]))^2)/(2*exp(as.numeric(sy_lo_bp_c[i,]))^2))
    
    rm(fem.ef.bb,fem.ef.bp)
  }
  
  EW.rel.bb <- t(t(EW.bb)/colMeans(EW.bb))
  EW.rel.bb.c <- t(t(EW.bb.c)/colMeans(EW.bb.c))
  EW.rel.bp <- t(t(EW.bp)/colMeans(EW.bp))
  EW.rel.bp.c <- t(t(EW.bp.c)/colMeans(EW.bp.c))
  
  ## Storing hatch dates 
  
  sy_mf_bb[i,] <- colMeans(EW.bb)
  sy_mf_bb_c[i,] <- colMeans(EW.bb.c)
  
  sy_mf_bp[i,] <- colMeans(EW.bp)
  sy_mf_bp_c[i,] <- colMeans(EW.bp.c)
  
  for(k in 1:ncol(EW.rel.bb)){
    sim.hd.cent <- sim.hd.bb[,k]-mean(sim.hd.bb[,k])
    mod.bb <- lm(EW.rel.bb[,k]~sim.hd.cent+I(sim.hd.cent^2))
    mod.bb.c <- lm(EW.rel.bb.c[,k]~sim.hd.cent+I(sim.hd.cent^2))
    
    sy_d_bb[i,k] <- as.numeric(coefficients(mod.bb)[2])
    sy_d_bb_c[i,k] <- as.numeric(coefficients(mod.bb.c)[2])
    sy_d2_bb[i,k] <- as.numeric(coefficients(mod.bb)[3])
    sy_d2_bb_c[i,k] <- as.numeric(coefficients(mod.bb.c)[3])
    
    rm(mod.bb,mod.bb.c,sim.hd.cent)
  }
  
  
  for(l in 1:ncol(EW.rel.bp)){
    sim.hd.cent <- sim.hd.bp[,l]-mean(sim.hd.bp[,l])
    mod.bp <- lm(EW.rel.bp[,l]~sim.hd.cent+I(sim.hd.cent^2))
    mod.bp.c <- lm(EW.rel.bp.c[,l]~sim.hd.cent+I(sim.hd.cent^2))
    
    sy_d_bp[i,l] <- as.numeric(coefficients(mod.bp)[2])
    sy_d_bp_c[i,l] <- as.numeric(coefficients(mod.bp.c)[2])
    sy_d2_bp[i,l] <- as.numeric(coefficients(mod.bp)[3])
    sy_d2_bp_c[i,l] <- as.numeric(coefficients(mod.bp.c)[3])
    
    rm(mod.bp,mod.bp.c,sim.hd.cent)
  }
  
  rm(EW.bb,EW.bb.c,EW.bp,EW.bp.c,sim.hd.bb,EW.rel.bb,EW.rel.bb.c,EW.rel.bp,EW.rel.bp.c,sim.hd.bb,sim.hd.bp)
  
  setTxtProgressBar(pb, i)
} 

close(pb)  


write.csv(sy_d_bb,"sy_d_bb.csv",row.names=F)
write.csv(sy_d_bb_c,"sy_d_bb_c.csv",row.names=F)
write.csv(sy_d_bp,"sy_d_bp.csv",row.names=F)
write.csv(sy_d_bp_c,"sy_d_bp_c.csv",row.names=F)

write.csv(sy_d2_bb,"sy_d2_bb.csv",row.names=F)
write.csv(sy_d2_bb_c,"sy_d2_bb_c.csv",row.names=F)
write.csv(sy_d2_bp,"sy_d2_bp.csv",row.names=F)
write.csv(sy_d2_bp_c,"sy_d2_bp_c.csv",row.names=F)

write.csv(sy_mf_bb,"sy_mf_bb.csv",row.names=F)
write.csv(sy_mf_bb_c,"sy_mf_bb_c.csv",row.names=F)
write.csv(sy_mf_bp,"sy_mf_bp.csv",row.names=F)
write.csv(sy_mf_bp_c,"sy_mf_bp_c.csv",row.names=F)

# sy_d_bb <- read.csv("sy_d_bb.csv")
# sy_d_bb_c <- read.csv("sy_d_bb_c.csv")
# sy_d_bp <- read.csv("sy_d_bp.csv")
# sy_d_bp_c <- read.csv("sy_d_bp_c.csv")
# 
# sy_d2_bb <- read.csv("sy_d2_bb.csv")
# sy_d2_bb_c <- read.csv("sy_d2_bb_c.csv")
# sy_d2_bp <- read.csv("sy_d2_bp.csv")
# sy_d2_bp_c <- read.csv("sy_d2_bp_c.csv")
# 
# sy_mf_bb <- read.csv("sy_mf_bb.csv")
# sy_mf_bb_c <- read.csv("sy_mf_bb_c.csv")
# sy_mf_bp <- read.csv("sy_mf_bp.csv")
# sy_mf_bp_c <- read.csv("sy_mf_bp_c.csv")

#### Variance decomposition in space and time ####

mf_d_var <- data.frame(matrix(NA,nrow=nrow(opt_bb),ncol=8*3))
colnames(mf_d_var) <- c("s_mf_bb","s_mf_bb_c","s_mf_bp","s_mf_bp_c","s_d_bb","s_d_bb_c","s_d_bp","s_d_bp_c",
                        "y_mf_bb","y_mf_bb_c","y_mf_bp","y_mf_bp_c","y_d_bb","y_d_bb_c","y_d_bp","y_d_bp_c",
                        "sy_mf_bb","sy_mf_bb_c","sy_mf_bp","sy_mf_bp_c","sy_d_bb","sy_d_bb_c","sy_d_bp","sy_d_bp_c")

for(i in 1:nrow(opt_bb)){
  
  SY_id$mf <- as.numeric(sy_mf_bb[i,])
  SY_id_p$mf <- as.numeric(sy_mf_bp[i,])
  SY_id$d <- as.numeric(sy_d_bb[i,])
  SY_id_p$d <- as.numeric(sy_d_bp[i,])
  
  SY_id$mf_c <- as.numeric(sy_mf_bb_c[i,])
  SY_id_p$mf_c <- as.numeric(sy_mf_bp_c[i,])
  SY_id$d_c <- as.numeric(sy_d_bb_c[i,])
  SY_id_p$d_c <- as.numeric(sy_d_bp_c[i,])
  
  mod_mf_bb <- lmer(mf~1+(1|site)+(1|year),data=SY_id)
  mod_mf_bp <- lmer(mf~1+(1|site)+(1|year),data=SY_id_p)
  mod_d_bb <- lmer(d~1+(1|site)+(1|year),data=SY_id)
  mod_d_bp <- lmer(d~1+(1|site)+(1|year),data=SY_id_p)
  
  mod_mf_bb_c <- lmer(mf_c~1+(1|site)+(1|year),data=SY_id)
  mod_mf_bp_c <- lmer(mf_c~1+(1|site)+(1|year),data=SY_id_p)
  mod_d_bb_c <- lmer(d_c~1+(1|site)+(1|year),data=SY_id)
  mod_d_bp_c <- lmer(d_c~1+(1|site)+(1|year),data=SY_id_p)
  
  mf_d_var$s_mf_bb[i] <- as.data.frame(VarCorr(mod_mf_bb))$vcov[1]
  mf_d_var$y_mf_bb[i] <- as.data.frame(VarCorr(mod_mf_bb))$vcov[2]
  mf_d_var$sy_mf_bb[i] <- as.data.frame(VarCorr(mod_mf_bb))$vcov[3]
  
  mf_d_var$s_mf_bb_c[i] <- as.data.frame(VarCorr(mod_mf_bb_c))$vcov[1]
  mf_d_var$y_mf_bb_c[i] <- as.data.frame(VarCorr(mod_mf_bb_c))$vcov[2]
  mf_d_var$sy_mf_bb_c[i] <- as.data.frame(VarCorr(mod_mf_bb_c))$vcov[3]
  
  mf_d_var$s_mf_bp[i] <- as.data.frame(VarCorr(mod_mf_bp))$vcov[1]
  mf_d_var$y_mf_bp[i] <- as.data.frame(VarCorr(mod_mf_bp))$vcov[2]
  mf_d_var$sy_mf_bp[i] <- as.data.frame(VarCorr(mod_mf_bp))$vcov[3]
  
  mf_d_var$s_mf_bp_c[i] <- as.data.frame(VarCorr(mod_mf_bp_c))$vcov[1]
  mf_d_var$y_mf_bp_c[i] <- as.data.frame(VarCorr(mod_mf_bp_c))$vcov[2]
  mf_d_var$sy_mf_bp_c[i] <- as.data.frame(VarCorr(mod_mf_bp_c))$vcov[3] 
  
  mf_d_var$s_d_bb[i] <- as.data.frame(VarCorr(mod_d_bb))$vcov[1]
  mf_d_var$y_d_bb[i] <- as.data.frame(VarCorr(mod_d_bb))$vcov[2]
  mf_d_var$sy_d_bb[i] <- as.data.frame(VarCorr(mod_d_bb))$vcov[3]
  
  mf_d_var$s_d_bb_c[i] <- as.data.frame(VarCorr(mod_d_bb_c))$vcov[1]
  mf_d_var$y_d_bb_c[i] <- as.data.frame(VarCorr(mod_d_bb_c))$vcov[2]
  mf_d_var$sy_d_bb_c[i] <- as.data.frame(VarCorr(mod_d_bb_c))$vcov[3]
  
  mf_d_var$s_d_bp[i] <- as.data.frame(VarCorr(mod_d_bp))$vcov[1]
  mf_d_var$y_d_bp[i] <- as.data.frame(VarCorr(mod_d_bp))$vcov[2]
  mf_d_var$sy_d_bp[i] <- as.data.frame(VarCorr(mod_d_bp))$vcov[3]
  
  mf_d_var$s_d_bp_c[i] <- as.data.frame(VarCorr(mod_d_bp_c))$vcov[1]
  mf_d_var$y_d_bp_c[i] <- as.data.frame(VarCorr(mod_d_bp_c))$vcov[2]
  mf_d_var$sy_d_bp_c[i] <- as.data.frame(VarCorr(mod_d_bp_c))$vcov[3]
  
  print(i)   
}

mf_d_var$tot_mf_bb <- mf_d_var$s_mf_bb+mf_d_var$y_mf_bb+mf_d_var$sy_mf_bb
mf_d_var$tot_mf_bb_c<- mf_d_var$s_mf_bb_c+mf_d_var$y_mf_bb_c+mf_d_var$sy_mf_bb_c
mf_d_var$tot_mf_bp <- mf_d_var$s_mf_bp+mf_d_var$y_mf_bp+mf_d_var$sy_mf_bp
mf_d_var$tot_mf_bp_c<- mf_d_var$s_mf_bp_c+mf_d_var$y_mf_bp_c+mf_d_var$sy_mf_bp_c

mf_d_var$tot_d_bb <- mf_d_var$s_d_bb+mf_d_var$y_d_bb+mf_d_var$sy_d_bb
mf_d_var$tot_d_bb_c<- mf_d_var$s_d_bb_c+mf_d_var$y_d_bb_c+mf_d_var$sy_d_bb_c
mf_d_var$tot_d_bp <- mf_d_var$s_d_bp+mf_d_var$y_d_bp+mf_d_var$sy_d_bp
mf_d_var$tot_d_bp_c<- mf_d_var$s_d_bp_c+mf_d_var$y_d_bp_c+mf_d_var$sy_d_bp_c

mf_d_var$prop_cat_mf_bb <- mf_d_var$tot_mf_bb_c/mf_d_var$tot_mf_bb
mf_d_var$prop_cat_mf_bp <- mf_d_var$tot_mf_bp_c/mf_d_var$tot_mf_bp
mf_d_var$prop_cat_d_bb <- mf_d_var$tot_d_bb_c/mf_d_var$tot_d_bb
mf_d_var$prop_cat_d_bp <- mf_d_var$tot_d_bp_c/mf_d_var$tot_d_bp

mf_d_var$prop_site_mf_bb <- mf_d_var$s_mf_bb/mf_d_var$tot_mf_bb
mf_d_var$prop_site_mf_bp <- mf_d_var$s_mf_bp/mf_d_var$tot_mf_bp
mf_d_var$prop_site_d_bb <- mf_d_var$s_d_bb/mf_d_var$tot_d_bb
mf_d_var$prop_site_d_bp <- mf_d_var$s_d_bp/mf_d_var$tot_d_bp

mf_d_var$prop_year_mf_bb <- mf_d_var$y_mf_bb/mf_d_var$tot_mf_bb
mf_d_var$prop_year_mf_bp <- mf_d_var$y_mf_bp/mf_d_var$tot_mf_bp
mf_d_var$prop_year_d_bb <- mf_d_var$y_d_bb/mf_d_var$tot_d_bb
mf_d_var$prop_year_d_bp <- mf_d_var$y_d_bp/mf_d_var$tot_d_bp

mf_d_var$prop_siteyear_mf_bb <- mf_d_var$sy_mf_bb/mf_d_var$tot_mf_bb
mf_d_var$prop_siteyear_mf_bp <- mf_d_var$sy_mf_bp/mf_d_var$tot_mf_bp
mf_d_var$prop_siteyear_d_bb <- mf_d_var$sy_d_bb/mf_d_var$tot_d_bb
mf_d_var$prop_siteyear_d_bp <- mf_d_var$sy_d_bp/mf_d_var$tot_d_bp

mf_d_var$prop_site_mf_bb_c <- mf_d_var$s_mf_bb_c/mf_d_var$tot_mf_bb_c
mf_d_var$prop_site_mf_bp_c <- mf_d_var$s_mf_bp_c/mf_d_var$tot_mf_bp_c
mf_d_var$prop_site_d_bb_c <- mf_d_var$s_d_bb_c/mf_d_var$tot_d_bb_c
mf_d_var$prop_site_d_bp_c <- mf_d_var$s_d_bp_c/mf_d_var$tot_d_bp_c

mf_d_var$prop_year_mf_bb_c <- mf_d_var$y_mf_bb_c/mf_d_var$tot_mf_bb_c
mf_d_var$prop_year_mf_bp_c <- mf_d_var$y_mf_bp_c/mf_d_var$tot_mf_bp_c
mf_d_var$prop_year_d_bb_c <- mf_d_var$y_d_bb_c/mf_d_var$tot_d_bb_c
mf_d_var$prop_year_d_bp_c <- mf_d_var$y_d_bp_c/mf_d_var$tot_d_bp_c

mf_d_var$prop_siteyear_mf_bb_c <- mf_d_var$sy_mf_bb_c/mf_d_var$tot_mf_bb_c
mf_d_var$prop_siteyear_mf_bp_c <- mf_d_var$sy_mf_bp_c/mf_d_var$tot_mf_bp_c
mf_d_var$prop_siteyear_d_bb_c <- mf_d_var$sy_d_bb_c/mf_d_var$tot_d_bb_c
mf_d_var$prop_siteyear_d_bp_c <- mf_d_var$sy_d_bp_c/mf_d_var$tot_d_bp_c

mf_d_var$prop_site_mf_bb_c_tot <- (mf_d_var$s_mf_bb_c/mf_d_var$tot_mf_bb_c)*mf_d_var$prop_cat_mf_bb
mf_d_var$prop_site_mf_bp_c_tot <- (mf_d_var$s_mf_bp_c/mf_d_var$tot_mf_bp_c)*mf_d_var$prop_cat_mf_bp
mf_d_var$prop_site_d_bb_c_tot <- (mf_d_var$s_d_bb_c/mf_d_var$tot_d_bb_c)*mf_d_var$prop_cat_d_bb
mf_d_var$prop_site_d_bp_c_tot <- (mf_d_var$s_d_bp_c/mf_d_var$tot_d_bp_c)*mf_d_var$prop_cat_d_bp

mf_d_var$prop_year_mf_bb_c_tot <- (mf_d_var$y_mf_bb_c/mf_d_var$tot_mf_bb_c)*mf_d_var$prop_cat_mf_bb
mf_d_var$prop_year_mf_bp_c_tot <- (mf_d_var$y_mf_bp_c/mf_d_var$tot_mf_bp_c)*mf_d_var$prop_cat_mf_bp
mf_d_var$prop_year_d_bb_c_tot <- (mf_d_var$y_d_bb_c/mf_d_var$tot_d_bb_c)*mf_d_var$prop_cat_d_bb
mf_d_var$prop_year_d_bp_c_tot <- (mf_d_var$y_d_bp_c/mf_d_var$tot_d_bp_c)*mf_d_var$prop_cat_d_bp

mf_d_var$prop_siteyear_mf_bb_c_tot <- (mf_d_var$sy_mf_bb_c/mf_d_var$tot_mf_bb_c)*mf_d_var$prop_cat_mf_bb
mf_d_var$prop_siteyear_mf_bp_c_tot <- (mf_d_var$sy_mf_bp_c/mf_d_var$tot_mf_bp_c)*mf_d_var$prop_cat_mf_bp
mf_d_var$prop_siteyear_d_bb_c_tot <- (mf_d_var$sy_d_bb_c/mf_d_var$tot_d_bb_c)*mf_d_var$prop_cat_d_bb
mf_d_var$prop_siteyear_d_bp_c_tot <- (mf_d_var$sy_d_bp_c/mf_d_var$tot_d_bp_c)*mf_d_var$prop_cat_d_bp

sy_prop <- expand.grid(model=c("Probability of success","Number fledged"),
                       coef=c("Mean fitness","Directional selection"),
                       term=c("Site","Year","Site-year"))

sy_prop$mean <- c(mean(mf_d_var$prop_site_mf_bb),
                  mean(mf_d_var$prop_site_mf_bp),
                  mean(mf_d_var$prop_site_d_bb),
                  mean(mf_d_var$prop_site_d_bp),
                  mean(mf_d_var$prop_year_mf_bb),
                  mean(mf_d_var$prop_year_mf_bp),
                  mean(mf_d_var$prop_year_d_bb),
                  mean(mf_d_var$prop_year_d_bp),
                  mean(mf_d_var$prop_siteyear_mf_bb),
                  mean(mf_d_var$prop_siteyear_mf_bp),
                  mean(mf_d_var$prop_siteyear_d_bb),
                  mean(mf_d_var$prop_siteyear_d_bp))

sy_prop$lci <- c(posterior_interval(matrix(mf_d_var$prop_site_mf_bb), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_site_mf_bp), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_site_d_bb), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_site_d_bp), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_year_mf_bb), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_year_mf_bp), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_year_d_bb), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_year_d_bp), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bb), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bp), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_d_bb), prob=0.95)[1],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_d_bp), prob=0.95)[1])

sy_prop$uci <- c(posterior_interval(matrix(mf_d_var$prop_site_mf_bb), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_site_mf_bp), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_site_d_bb), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_site_d_bp), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_year_mf_bb), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_year_mf_bp), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_year_d_bb), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_year_d_bp), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bb), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bp), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_d_bb), prob=0.95)[2],
                 posterior_interval(matrix(mf_d_var$prop_siteyear_d_bp), prob=0.95)[2])

sy_prop_c_tot <- expand.grid(model=c("Probability of success","Number fledged"),
                             coef=c("Mean fitness","Directional selection"),
                             term=c("Site","Year","Site-year"))

sy_prop_c_tot$mean <- c(mean(mf_d_var$prop_site_mf_bb_c_tot),
                        mean(mf_d_var$prop_site_mf_bp_c_tot),
                        mean(mf_d_var$prop_site_d_bb_c_tot),
                        mean(mf_d_var$prop_site_d_bp_c_tot),
                        mean(mf_d_var$prop_year_mf_bb_c_tot),
                        mean(mf_d_var$prop_year_mf_bp_c_tot),
                        mean(mf_d_var$prop_year_d_bb_c_tot),
                        mean(mf_d_var$prop_year_d_bp_c_tot),
                        mean(mf_d_var$prop_siteyear_mf_bb_c_tot),
                        mean(mf_d_var$prop_siteyear_mf_bp_c_tot),
                        mean(mf_d_var$prop_siteyear_d_bb_c_tot),
                        mean(mf_d_var$prop_siteyear_d_bp_c_tot))

sy_prop_c_tot$lci <- c(posterior_interval(matrix(mf_d_var$prop_site_mf_bb_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_site_mf_bp_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_site_d_bb_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_site_d_bp_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_year_mf_bb_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_year_mf_bp_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_year_d_bb_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_year_d_bp_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bb_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bp_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_d_bb_c_tot), prob=0.95)[1],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_d_bp_c_tot), prob=0.95)[1])

sy_prop_c_tot$uci <- c(posterior_interval(matrix(mf_d_var$prop_site_mf_bb_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_site_mf_bp_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_site_d_bb_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_site_d_bp_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_year_mf_bb_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_year_mf_bp_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_year_d_bb_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_year_d_bp_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bb_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_mf_bp_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_d_bb_c_tot), prob=0.95)[2],
                       posterior_interval(matrix(mf_d_var$prop_siteyear_d_bp_c_tot), prob=0.95)[2])


sy_prop$Var <- "Total" 
sy_prop_c_tot$Var <- "Caterpillar" 
all_prop <- rbind(sy_prop,sy_prop_c_tot)
all_prop$Var <- factor(all_prop$Var,levels=c("Total", "Caterpillar"))

(full_plot <- ggplot(all_prop, aes(alpha=term, fill=Var, y=mean, x=Var))+ 
    geom_bar(position="stack", stat="identity")+
    scale_fill_viridis(discrete=T,option="A",begin=0.45,end=0.65, name="Component")+
    facet_nested_wrap(coef~model, ncol=4)+
    xlab("")+
    ylab("")+
    ggtitle("")+
    coord_cartesian(y=c(0.02,1))+
    theme_classic()+
    scale_alpha_discrete(range = c(0.3, 1.0), name="Term")+
    theme(axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          strip.text.x = element_blank(),
          plot.margin = margin(t = 0.5,  # Top margin
                               r = 0.8,  # Right margin
                               b = 0.8,  # Bottom margin
                               l = 0.8,  # Left margin
                               unit = "cm")))

grid.text("Mean fitness", x=0.28,y=0.93,gp=gpar(fontface="bold"))
grid.text("Directional selection", x=0.63,y=0.93,gp=gpar(fontface="bold"))
grid.text("a", x=0.14,y=0.93,gp=gpar(cex=0.9))
grid.text("b", x=0.48,y=0.93,gp=gpar(cex=0.9))
grid.text("Proportion of total variance",rot=90,x=0.04,y=0.5,gp=gpar(cex=0.9))
grid.text("Probability of\nsuccess", x=0.195,y=0.06,gp=gpar(cex=0.9))
grid.text("Number\nfledged", x=0.365,y=0.06,gp=gpar(cex=0.9))
grid.text("Probability of\nsuccess", x=0.535,y=0.06,gp=gpar(cex=0.9))
grid.text("Number\nfledged", x=0.71,y=0.06,gp=gpar(cex=0.9))
grid.lines(x=unit(c(0.45,0.45),"npc"), y=unit(c(0.03,0.95),"npc"),gp=gpar(lty="dashed",width=0.5))
#7*4


# Violin plots
prop_cat_violin <- mf_d_var[,c("prop_cat_mf_bb","prop_cat_mf_bp","prop_cat_d_bb","prop_cat_d_bp" )]
for(i in 1:ncol(prop_cat_violin)){
  prop_cat_violin[,i] <- ifelse(prop_cat_violin[,i]>=posterior_interval(matrix(prop_cat_violin[,i]),prob=0.95)[1]&prop_cat_violin[,i]<=posterior_interval(matrix(prop_cat_violin[,i]),prob=0.95)[2],prop_cat_violin[,i],NA)
}
prop_cat_violin <- gather(prop_cat_violin,value="PseudoR2",key="Term",1:4)
prop_cat_violin$Model <- substr(prop_cat_violin$Term,nchar(prop_cat_violin$Term)-1,nchar(prop_cat_violin$Term))
prop_cat_violin$Coef <- substr(prop_cat_violin$Term,nchar(prop_cat_violin$Term)-4,nchar(prop_cat_violin$Term)-2)

prop_cat <- expand.grid(Model=c("bb","bp"),
                        Coef=c("mf_","_d_"))

prop_cat$Mean <- c(mean(mf_d_var$prop_cat_mf_bb),
                   mean(mf_d_var$prop_cat_mf_bp),
                   mean(mf_d_var$prop_cat_d_bb),
                   mean(mf_d_var$prop_cat_d_bp))

prop_cat$lci <- c(posterior_interval(matrix(mf_d_var$prop_cat_mf_bb), prob=0.95)[1],
                  posterior_interval(matrix(mf_d_var$prop_cat_mf_bp), prob=0.95)[1],
                  posterior_interval(matrix(mf_d_var$prop_cat_d_bb), prob=0.95)[1],
                  posterior_interval(matrix(mf_d_var$prop_cat_d_bp), prob=0.95)[1])

prop_cat$uci <- c(posterior_interval(matrix(mf_d_var$prop_cat_mf_bb), prob=0.95)[2],
                  posterior_interval(matrix(mf_d_var$prop_cat_mf_bp), prob=0.95)[2],
                  posterior_interval(matrix(mf_d_var$prop_cat_d_bb), prob=0.95)[2],
                  posterior_interval(matrix(mf_d_var$prop_cat_d_bp), prob=0.95)[2])

supp.labs <- c("Mean fitness","Directional selection")
names(supp.labs) <- c("mf_","_d_")

prop_cat$Coef <- factor(prop_cat$Coef, levels=c("mf_","_d_"))
prop_cat_violin$Coef <- factor(prop_cat_violin$Coef, levels=c("mf_","_d_"))

ggplot()+
  geom_violin(data=prop_cat_violin,aes(Model,PseudoR2, fill=Model))+
  scale_fill_discrete(palette=c("grey90","grey50"))+
  geom_point(data=prop_cat,aes(Model,Mean))+
  facet_wrap(~Coef, labeller=labeller(Coef=supp.labs))+
  theme_classic()+
  ylab("Prop. of variance from caterpillar peak")+
  xlab("Response")+
  scale_x_discrete(labels = c("Prob. of success","Number fledged"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) #5*5

grid.text(label="a",x=(0.13), y=(0.965), rot=0, gp=gpar(fontsize=10))
grid.text(label="b",x=(0.57), y=(0.965), rot=0, gp=gpar(fontsize=10))


