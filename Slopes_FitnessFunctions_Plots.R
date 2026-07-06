rm(list=ls())
setwd("")
source("OrganisingDataframes.R")
source("ResultsFunctions.R")
library(rstan)
library(rstanarm)
library(MCMCglmm)
library(dplyr)
library(forcats)
library(gridExtra)
library(grid)
library(cowplot)
library(viridis)

mod <- readRDS("FullMMH_mod.rds")
#rstan::get_num_divergent(mod)


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


### Slope mean and CIs ###

slopes <- stanpost(model=mod, parameters=c("mu_theta_bb","ls_logomega_bb","lm_logomega_bb","lm_logWmax_bb","mu_theta_bp","ls_logomega_bp","lm_logomega_bp","lm_logWmax_bp"))
betas <- as.data.frame(posterior_interval(as.matrix(slopes), prob=0.95))
betas$mean <- colMeans(slopes)
for(i in 1:ncol(slopes)){
  betas$PropPos[i] <- length(which(slopes[,i]>0))/nrow(slopes)
}

betas <- as.data.frame(as.matrix(round(betas,2)))
#write.csv(betas,"slope_coefs.csv")

### Site-year predictions for fitness function parameters

SY_id$mu <- colMeans(cater_mu_bb)+8+146
SY_id$ls <- exp(colMeans(cater_ls_bb)+2.5)
SY_id$lm <- exp(colMeans(cater_lm_bb)-2.5)

SY_id_p$mu <- colMeans(cater_mu_bp)+8+146
SY_id_p$ls <- exp(colMeans(cater_ls_bp)+2.5)
SY_id_p$lm <- exp(colMeans(cater_lm_bp)-2.5)

for(i in 1:nrow(SY_id)){
  SY_id$th_bb[i] <- mean(th_bb[,1] + 
                           mu_th_bb[,1]*cater_mu_bb[,c(SY_id$siteyear_id[i])]+
                           bb_s[,c(SY_id$site_id[i])]+
                           bb_y[,c(SY_id$year_id[i])]+
                           bb_sy[,c(SY_id$siteyear_id[i])])+146
  
  SY_id$lo_bb[i] <- mean(exp(lo_bb[,1] + 
                               ls_lo_bb[,1]*cater_ls_bb[,c(SY_id$siteyear_id[i])]+
                               lm_lo_bb[,1]*cater_lm_bb[,c(SY_id$siteyear_id[i])]+
                               bb_s[,c(SY_id$site_id[i]+44)]+
                               bb_y[,c(SY_id$year_id[i]+11)]+
                               bb_sy[,c(SY_id$siteyear_id[i]+369)]))
  
  mean_bb <- lw_bb[,1] + 
    lm_lw_bb[,1]*cater_lm_bb[,c(SY_id$siteyear_id[i])]+
    bb_s[,c(SY_id$site_id[i]+44+44)]+
    bb_y[,c(SY_id$year_id[i]+11+11)]+
    bb_sy[,c(SY_id$siteyear_id[i]+369+369)]
  var_bb <- bb_fem_sd[,1]^2
  
  SY_id$lw_bb[i] <- mean(diggle_approx(mean=mean_bb,var=var_bb))
  
  rm(mean_bb,var_bb)
  
}


for(i in 1:nrow(SY_id_p)){  
  SY_id_p$th_bp[i] <- mean(th_bp[,1] + 
                             mu_th_bp[,1]*cater_mu_bp[,c(SY_id_p$siteyear_id[i])]+
                             bp_s[,c(SY_id_p$site_id[i])]+
                             bp_y[,c(SY_id_p$year_id[i])]+
                             bp_sy[,c(SY_id_p$siteyear_id[i])])+146
  
  SY_id_p$lo_bp[i] <- mean(exp(lo_bp[,1] + 
                                 ls_lo_bp[,1]*cater_ls_bp[,c(SY_id_p$siteyear_id[i])]+
                                 lm_lo_bp[,1]*cater_lm_bp[,c(SY_id_p$siteyear_id[i])]+
                                 bp_s[,c(SY_id_p$site_id[i]+44)]+
                                 bp_y[,c(SY_id_p$year_id[i]+11)]+
                                 bp_sy[,c(SY_id_p$siteyear_id[i]+359)]))
  
  SY_id_p$lw_bp[i] <- mean(exp(lw_bp[,1] + 
                                 lm_lw_bp[,1]*cater_lm_bp[,c(SY_id_p$siteyear_id[i])]+
                                 bp_s[,c(SY_id_p$site_id[i]+44+44)]+
                                 bp_y[,c(SY_id_p$year_id[i]+11+11)]+
                                 bp_sy[,c(SY_id_p$siteyear_id[i]+359+359)]+
                                 (bp_fem_sd[,1]^2)/2))
}

### Pseudo-R^2 for each slope
# R^2(z by x) = beta^2*var(x)/(beta^2*var(x)+var(z:s)+var(z:y)+var(z:sy))

r2 <- data.frame(parameter=c("mu_th_bb","ls_lo_bb","lm_lo_bb","lm_lw_bb","mu_th_bp","ls_lo_bp","lm_lo_bp","lm_lw_bp"),
                 mean=NA,
                 loci=NA,
                 upci=NA)

pseudo_r2 <- function(beta,var_x,other_var){
  p_r2 <- (var_x*beta^2)/(var_x*beta^2+other_var)
  summary <- c(mean(p_r2),
               posterior_interval(as.matrix(p_r2), prob=0.95)[1],
               posterior_interval(as.matrix(p_r2), prob=0.95)[2])
  return(list(p_r2,summary))
}

r2[1,c(2:4)] <- pseudo_r2(beta=mu_th_bb[,1],
                          var_x=mean(apply(cater_mu_bb,1,var)),
                          other_var=bb_s_sd[,1]^2+bb_y_sd[,1]^2+bb_sy_sd[,1]^2)[[2]]

r2[2,c(2:4)] <- pseudo_r2(beta=ls_lo_bb[,1],
                          var_x=mean(apply(cater_ls_bb,1,var)),
                          other_var=bb_s_sd[,2]^2+bb_y_sd[,2]^2+bb_sy_sd[,2]^2)[[2]]

r2[3,c(2:4)] <- pseudo_r2(beta=lm_lo_bb[,1],
                          var_x=mean(apply(cater_ls_bb,1,var)),
                          other_var=bb_s_sd[,2]^2+bb_y_sd[,2]^2+bb_sy_sd[,2]^2)[[2]]

r2[4,c(2:4)] <- pseudo_r2(beta=lm_lw_bb[,1],
                          var_x=mean(apply(cater_lm_bb,1,var)),
                          other_var=bb_s_sd[,3]^2+bb_y_sd[,3]^2+bb_sy_sd[,3]^2)[[2]]

r2[5,c(2:4)] <- pseudo_r2(beta=mu_th_bp[,1],
                          var_x=mean(apply(cater_mu_bp,1,var)),
                          other_var=bp_s_sd[,1]^2+bp_y_sd[,1]^2+bp_sy_sd[,1]^2)[[2]]

r2[6,c(2:4)] <- pseudo_r2(beta=ls_lo_bp[,1],
                          var_x=mean(apply(cater_ls_bp,1,var)),
                          other_var=bp_s_sd[,2]^2+bp_y_sd[,2]^2+bp_sy_sd[,2]^2)[[2]]

r2[7,c(2:4)] <- pseudo_r2(beta=lm_lo_bp[,1],
                          var_x=mean(apply(cater_ls_bp,1,var)),
                          other_var=bp_s_sd[,2]^2+bp_y_sd[,2]^2+bp_sy_sd[,2]^2)[[2]]

r2[8,c(2:4)] <- pseudo_r2(beta=lm_lw_bp[,1],
                          var_x=mean(apply(cater_lm_bp,1,var)),
                          other_var=bp_s_sd[,3]^2+bp_y_sd[,3]^2+bp_sy_sd[,3]^2)[[2]]

r2[,2:4] <- round(r2[,2:4]*100,2)
#write.csv(r2,"slope_r2s.csv",row.names = F)

### Quantiles and colours

### Plot of change in optimum with change in caterpillar timing      
# theta_bb_obs = theta_bb+mu_theta_bb*c_mu_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,1]+year_effs_bb[year_id_bb,1]+siteyear_effs_bb[siteyear_id_bb,1];
# theta_bp_obs = theta_bp+mu_theta_bp*c_mu_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,1]+year_effs_bp[year_id_bp,1]+siteyear_effs_bp[siteyear_id_bp,1];

mu_cent <- seq(min(colMeans(cater_mu_bb)),max(colMeans(cater_mu_bb)),0.2)
mu_uncent <- mu_cent+8+146

thta_pred <- data.frame(mu_cent=mu_cent,
                        mu_uncent=mu_uncent,
                        bb_mean=NA,
                        bb_loci=NA,
                        bb_upci=NA,
                        bp_mean=NA,
                        bp_loci=NA,
                        bp_upci=NA)

for(i in 1:nrow(thta_pred)){
  X_bb <- th_bb+mu_th_bb*thta_pred$mu_cent[i]
  X_bp <- th_bp+mu_th_bp*thta_pred$mu_cent[i]
  
  thta_pred$bb_mean[i] <- mean(X_bb[,1])
  thta_pred$bb_loci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[1]
  thta_pred$bb_upci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[2]
  
  thta_pred$bp_mean[i] <- mean(X_bp[,1])
  thta_pred$bp_loci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[1]
  thta_pred$bp_upci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[2]
  
  rm(X_bb,X_bp)
}

thta_pred$bb_mean_uncent <- thta_pred$bb_mean+146
thta_pred$bb_loci_uncent <- thta_pred$bb_loci+146
thta_pred$bb_upci_uncent <- thta_pred$bb_upci+146

thta_pred$bp_mean_uncent <- thta_pred$bp_mean+146
thta_pred$bp_loci_uncent <- thta_pred$bp_loci+146
thta_pred$bp_upci_uncent <- thta_pred$bp_upci+146



### Plot of change in max fitness with change in caterpillar height
# logWmax_bb_obs = logWmax_bb + lm_logWmax_bb*c_lm_bb[siteyear_id_bb] + site_effs_bb[site_id_bb,3] + year_effs_bb[year_id_bb,3] + siteyear_effs_bb[siteyear_id_bb,3] + fem_effs_bb[fem_id_bb];
# logWmax_bp_obs = logWmax_bp + lm_logWmax_bp*c_lm_bp[siteyear_id_bp] + site_effs_bp[site_id_bp,3] + year_effs_bp[year_id_bp,3] + siteyear_effs_bp[siteyear_id_bp,3] + fem_effs_bp[fem_id_bp];

lm_cent <- seq(min(colMeans(cater_lm_bb)),max(colMeans(cater_lm_bb)),0.05)
lm_uncent <- lm_cent-2.5

wmax_pred <- data.frame(lm_cent=lm_cent,
                        lm_uncent=lm_uncent,
                        max_uncent=exp(lm_uncent),
                        bb_mean=NA,
                        bb_loci=NA,
                        bb_upci=NA,
                        bp_mean=NA,
                        bp_loci=NA,
                        bp_upci=NA,
                        bb_omga_mean=NA,
                        bb_omga_loci=NA,
                        bb_omga_upci=NA,
                        bp_omga_mean=NA,
                        bp_omga_loci=NA,
                        bp_omga_upci=NA)

for(i in 1:nrow(wmax_pred)){
  mn_bb <- (lw_bb+lm_lw_bb*wmax_pred$lm_cent[i])[,1]
  var_bb <- bb_s_sd[,3]^2+bb_y_sd[,3]^2+bb_sy_sd[,3]^2+bb_fem_sd[,1]^2
  X_bb <- diggle_approx(mean=mn_bb, var=var_bb)
  X_bp <- exp(lw_bp+lm_lw_bp*wmax_pred$lm_cent[i]+(bp_s_sd[,3]^2+bp_y_sd[,3]^2+bp_sy_sd[,3]^2+bp_fem_sd[,1]^2)/2)
  
  wmax_pred$bb_mean[i] <- mean(X_bb)
  wmax_pred$bb_loci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[1]
  wmax_pred$bb_upci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[2]
  
  wmax_pred$bp_mean[i] <- mean(X_bp[,1])
  wmax_pred$bp_loci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[1]
  wmax_pred$bp_upci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[2]
  
  rm(mn_bb,var_bb,X_bb,X_bp)
}


### Plot of change in width of fitness function with change in caterpillar width
# omega_bb_obs = exp(logomega_bb+ls_logomega_bb*c_ls_bb[siteyear_id_bb]+lm_logomega_bb*c_lm_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,2]+year_effs_bb[year_id_bb,2]+siteyear_effs_bb[siteyear_id_bb,2]);
# omega_bp_obs = exp(logomega_bp+ls_logomega_bp*c_ls_bp[siteyear_id_bp]+lm_logomega_bp*c_lm_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,2]+year_effs_bp[year_id_bp,2]+siteyear_effs_bp[siteyear_id_bp,2]);

ls_cent <- seq(min(colMeans(cater_ls_bb)),max(colMeans(cater_ls_bb)),0.01)
ls_uncent <- ls_cent+2.5

omga_pred <- data.frame(ls_cent=ls_cent,
                        ls_uncent=ls_uncent,
                        sig_uncent=exp(ls_uncent),
                        bb_mean=NA,
                        bb_loci=NA,
                        bb_upci=NA,
                        bp_mean=NA,
                        bp_loci=NA,
                        bp_upci=NA)

for(i in 1:nrow(omga_pred)){
  X_bb <- exp(lo_bb+ls_lo_bb*omga_pred$ls_cent[i]+lm_lo_bb*mean(colMeans(cater_lm_bb))+(bb_s_sd[,2]^2+bb_y_sd[,2]^2+bb_sy_sd[,2]^2)/2) #
  X_bp <- exp(lo_bp+ls_lo_bp*omga_pred$ls_cent[i]+lm_lo_bp*mean(colMeans(cater_lm_bb))+(bp_s_sd[,2]^2+bp_y_sd[,2]^2+bp_sy_sd[,2]^2)/2) #
  
  omga_pred$bb_mean[i] <- mean(X_bb[,1])
  omga_pred$bb_loci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[1]
  omga_pred$bb_upci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[2]
  
  omga_pred$bp_mean[i] <- mean(X_bp[,1])
  omga_pred$bp_loci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[1]
  omga_pred$bp_upci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[2]
  
  rm(X_bb,X_bp)
}


### Plot of change in width of fitness function with change in caterpillar height
# omega_bb_obs = exp(logomega_bb+ls_logomega_bb*c_ls_bb[siteyear_id_bb]+lm_logomega_bb*c_lm_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,2]+year_effs_bb[year_id_bb,2]+siteyear_effs_bb[siteyear_id_bb,2]);
# omega_bp_obs = exp(logomega_bp+ls_logomega_bp*c_ls_bp[siteyear_id_bp]+lm_logomega_bp*c_lm_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,2]+year_effs_bp[year_id_bp,2]+siteyear_effs_bp[siteyear_id_bp,2]);

for(i in 1:nrow(wmax_pred)){
  X_bb <- exp(lo_bb+ls_lo_bb*mean(colMeans(cater_ls_bb))+lm_lo_bb*wmax_pred$lm_cent[i]+(bb_s_sd[,2]^2+bb_y_sd[,2]^2+bb_sy_sd[,2]^2)/2) #
  X_bp <- exp(lo_bp+ls_lo_bp*mean(colMeans(cater_ls_bb))+lm_lo_bp*wmax_pred$lm_cent[i]+(bp_s_sd[,2]^2+bp_y_sd[,2]^2+bp_sy_sd[,2]^2)/2) #
  
  wmax_pred$bb_omga_mean[i] <- mean(X_bb[,1])
  wmax_pred$bb_omga_loci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[1]
  wmax_pred$bb_omga_upci[i] <- posterior_interval(as.matrix(X_bb), prob=0.95)[2]
  
  wmax_pred$bp_omga_mean[i] <- mean(X_bp[,1])
  wmax_pred$bp_omga_loci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[1]
  wmax_pred$bp_omga_upci[i] <- posterior_interval(as.matrix(X_bp), prob=0.95)[2]
  
  rm(X_bb,X_bp)
}

lo_col <- "#148F2B"
th_col <- "#9D39BF" 
lw_col <- "#216DB8" 

(plot_th_mu_bb <- ggplot()+
    geom_abline(intercept=-10,slope=1,linetype="dashed",col="grey")+
    geom_point(data=SY_id,aes(mu,th_bb),alpha=0.2,size=0.5,col=th_col)+
    geom_ribbon(data=thta_pred,aes(x=mu_uncent,ymin=bb_loci_uncent,ymax=bb_upci_uncent),alpha=0.3,fill=th_col)+
    geom_line(data=thta_pred,aes(mu_uncent,bb_mean_uncent),col=th_col)+
    coord_cartesian(ylim=c(116,158))+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    theme_classic()+ 
    theme(axis.text.x = element_blank()))

(plot_th_mu_bp <- ggplot()+
    geom_abline(intercept=-10,slope=1,linetype="dashed",col="grey")+
    geom_point(data=SY_id_p,aes(mu,th_bp),alpha=0.2,size=0.5,col=th_col)+
    geom_ribbon(data=thta_pred,aes(x=mu_uncent,ymin=bp_loci_uncent,ymax=bp_upci_uncent),alpha=0.3,fill=th_col)+
    geom_line(data=thta_pred,aes(mu_uncent,bp_mean_uncent),col=th_col)+
    coord_cartesian(ylim=c(116,158))+
    guides(color = "none")+
    xlab("Caterpillar mean timing")+
    ylab("")+
    theme_classic())

(plot_lw_lm_bb <- ggplot()+
    geom_point(data=SY_id,aes(lm,lw_bb),alpha=0.2,size=0.5,col=lw_col)+
    geom_ribbon(data=wmax_pred,aes(x=max_uncent,ymin=bb_loci,ymax=bb_upci),alpha=0.3,fill=lw_col)+
    geom_line(data=wmax_pred,aes(max_uncent,bb_mean),col=lw_col)+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    theme_classic()+
    theme(axis.text.x = element_blank()))

(plot_lw_lm_bp <- ggplot()+
    geom_point(data=SY_id_p,aes(lm,lw_bp),alpha=0.2,size=0.5,col=lw_col)+
    geom_ribbon(data=wmax_pred,aes(x=max_uncent,ymin=bp_loci,ymax=bp_upci),alpha=0.3,fill=lw_col)+
    geom_line(data=wmax_pred,aes(max_uncent,bp_mean),col=lw_col)+
    scale_y_continuous(breaks=c(4,6,8,10,12,14),
                       labels=c(" 4"," 6"," 8"," 10","12","14"))+
    guides(color = "none")+
    xlab("Caterpillar height")+
    ylab("")+
    theme_classic())

(plot_lw_lm_bb_zoom <- ggplot()+
    geom_point(data=SY_id,aes(lm,lw_bb),alpha=0.2,size=0.5,col=lw_col)+
    geom_ribbon(data=wmax_pred,aes(x=max_uncent,ymin=bb_loci,ymax=bb_upci),alpha=0.3,fill=lw_col)+
    geom_line(data=wmax_pred,aes(max_uncent,bb_mean),col=lw_col)+
    coord_cartesian(xlim=c(0,0.7),ylim=c(0.518,1))+ # only 2% (9) siteyears have height>0.5
    guides(color = "none")+
    xlab("")+
    ylab("Max. probability of success")+
    theme_classic()+
    theme(axis.text.x = element_blank()))

(plot_lw_lm_bp_zoom <- ggplot()+
    geom_point(data=SY_id_p,aes(lm,lw_bp),alpha=0.2,size=0.5,col=lw_col)+
    geom_ribbon(data=wmax_pred,aes(x=max_uncent,ymin=bp_loci,ymax=bp_upci),alpha=0.3,fill=lw_col)+
    geom_line(data=wmax_pred,aes(max_uncent,bp_mean),col=lw_col)+
    coord_cartesian(xlim=c(0,0.7),ylim=c(4.13,10.5))+ # only 2% (9) siteyears have height>0.5
    scale_y_continuous(breaks=c(4,6,8,10,12,14),
                       labels=c(" 4"," 6"," 8"," 10","12","14"))+
    guides(color = "none")+
    xlab("Caterpillar height")+
    ylab("Max. number fledged")+
    theme_classic())

(plot_lo_ls_bb <- ggplot()+
    geom_point(data=SY_id,aes(ls,lo_bb),alpha=0.2,size=0.5,col=lo_col)+
    geom_ribbon(data=omga_pred,aes(x=sig_uncent,ymin=bb_loci,ymax=bb_upci),alpha=0.3,fill=lo_col)+
    geom_line(data=omga_pred,aes(sig_uncent,bb_mean),col=lo_col)+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    theme_classic()+ 
    theme(axis.text.x = element_blank()))

(plot_lo_ls_bp <- ggplot()+
    geom_point(data=SY_id_p,aes(ls,lo_bp),alpha=0.2,size=0.5,col=lo_col)+
    geom_ribbon(data=omga_pred,aes(x=sig_uncent,ymin=bp_loci,ymax=bp_upci),alpha=0.3,fill=lo_col)+
    geom_line(data=omga_pred,aes(sig_uncent,bp_mean),col=lo_col)+
    scale_y_continuous(breaks=c(20,40,60),
                       labels=c("  20","  40","  60"))+
    guides(color = "none")+
    xlab("Caterpillar width")+
    ylab("")+
    theme_classic())

(plot_lo_lm_bb <- ggplot()+
    geom_point(data=SY_id,aes(lm,lo_bb),alpha=0.2,size=0.5)+
    geom_ribbon(data=wmax_pred,aes(x=max_uncent,ymin=bb_omga_loci,ymax=bb_omga_upci),alpha=0.3,fill=1)+
    geom_line(data=wmax_pred,aes(max_uncent,bb_omga_mean))+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    theme_classic()+ 
    theme(axis.text.x = element_blank()))

(plot_lo_lm_bp <- ggplot()+
    geom_point(data=SY_id_p,aes(lm,lo_bp),alpha=0.2,size=0.5)+
    geom_ribbon(data=wmax_pred,aes(x=max_uncent,ymin=bp_omga_loci,ymax=bp_omga_upci),alpha=0.3,fill=1)+
    geom_line(data=wmax_pred,aes(max_uncent,bp_omga_mean))+
    scale_y_continuous(breaks=c(20,40,60),
                       labels=c("  20","  40","  60"))+
    guides(color = "none")+
    xlab("Caterpillar height")+
    ylab("")+
    theme_classic())



#### Average fitness function plots ####
bb_pred <- read.csv("bb_pred_function.csv")
bp_pred <- read.csv("bp_pred_function.csv")

bb_pred_0.5 <- subset(bb_pred, round(mu,2)==-0.18&para=="mu")
bp_pred_0.5 <- subset(bp_pred, round(mu,2)==-0.18&para=="mu")

(plot_av_bb <- ggplot()+
    geom_ribbon(data=bb_pred_0.5,aes(x=dat,ymin=lwci,ymax=upci),alpha=0.3)+
    geom_line(data=bb_pred_0.5,aes(dat,mean),lwd=1)+
    ylim(0,1)+
    xlab("")+
    ylab("")+
    theme_classic()+
    scale_x_continuous(breaks=c(130,150,170))+ 
    theme(axis.text.x = element_blank()))

(plot_av_bp <- ggplot()+
    geom_ribbon(data=bp_pred_0.5,aes(x=dat,ymin=lwci,ymax=upci),alpha=0.3)+
    geom_line(data=bp_pred_0.5,aes(dat,mean),lwd=1)+
    guides(color = "none")+
    xlab("Hatch date (1 = 1st Jan)")+
    ylab("")+
    theme_classic()+
    scale_x_continuous(breaks=c(130,150,170))+ 
    scale_y_continuous(breaks=c(2,4,6,8),
                       labels=c("     2","     4","     6","     8"),
                       limits = c(1,8.6)))

# Fitness function by caterpeak values - bernoulli
bb_pred_mu <- subset(bb_pred, para=="mu")
bb_pred_ls <- subset(bb_pred, para=="logsigma")
bb_pred_lm <- subset(bb_pred, para=="logmax")

textx <- 178
bb_texty <- 0.97
bp_texty <- 8.3
textsize <- 4


# Plot peaks for abundance estimate by date at quantiles of mu
(plot_mu_bb <- ggplot()+
    geom_line(data=bb_pred_mu, aes(dat, mean, linetype=as.factor(mu)),lwd=0.9, col="#9D39BF")+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"))+
    ylim(0,1)+
    theme_classic()+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+ 
    theme(axis.text.x = element_blank(),
          legend.position = "none"))

# Plot peaks for abundance estimate by date at quantiles of sigma
(plot_ls_bb <- ggplot()+
    geom_line(data=bb_pred_ls, aes(dat, mean, linetype=as.factor(logsigma)),lwd=0.9,col=lo_col)+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"))+
    ylim(0,1)+
    theme_classic()+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+ 
    theme(axis.text.x = element_blank(),
          legend.position = "none"))

# Plot peaks for abundance estimate by date at quantiles of max
(plot_lm_bb <- ggplot()+
    geom_line(data=bb_pred_lm, aes(dat, mean, linetype=as.factor(logmax)),lwd=0.9,col=lw_col)+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"))+
    ylim(0,1)+
    theme_classic()+
    guides(color = "none")+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+
    theme(legend.position = "none"))


# Fitness function by caterpeak values - poisson
bp_pred_mu <- subset(bp_pred, para=="mu")
bp_pred_ls <- subset(bp_pred, para=="logsigma")
bp_pred_lm <- subset(bp_pred, para=="logmax")


# Plot peaks for abundance estimate by date at quantiles of mu
(plot_mu_bp <- ggplot()+
    geom_line(data=bp_pred_mu, aes(dat, mean, linetype=as.factor(mu)),lwd=0.9,col="#9D39BF")+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"))+
    ylim(1,8.6)+
    guides(color = "none")+
    theme_classic()+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+ 
    theme(axis.text.x = element_blank(),
          legend.position = "none"))

# Plot peaks for abundance estimate by date at quantiles of sigma
(plot_ls_bp <- ggplot()+
    geom_line(data=bp_pred_ls, aes(dat, mean, linetype=as.factor(logsigma)),lwd=0.9,col=lo_col)+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"))+
    ylim(1,8.6)+
    guides(color = "none")+
    theme_classic()+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+ 
    theme(axis.text.x = element_blank(),
          legend.position = "none"))

# Plot peaks for abundance estimate by date at quantiles of max
(plot_lm_bp <- ggplot()+
    geom_line(data=bp_pred_lm, aes(dat, mean, linetype=as.factor(logmax)),lwd=0.9,col=lw_col)+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"))+
    ylim(1,8.6)+
    guides(color = "none")+
    theme_classic()+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+
    theme(legend.position = "none"))

bp_pred_lm$quantile <- ifelse(round(bp_pred_lm$logmax,1)==-2.3,"0.025",
                              ifelse(round(bp_pred_lm$logmax,1)==-1.2,"0.25",
                                     ifelse(round(bp_pred_lm$logmax,1)==-0.4,"0.5",
                                            ifelse(round(bp_pred_lm$logmax,1)==0.3,"0.75",
                                                   ifelse(round(bp_pred_lm$logmax,1)==1.9,"0.975",NA)))))
(plot_legend <- ggplot()+
    geom_line(data=bp_pred_lm, aes(dat, mean, linetype=as.factor(quantile)),lwd=0.9,col=1)+
    scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash","solid"), name="Caterpillar\nmetric quantile")+
    ylim(1,8.6)+
    guides(color = "none")+
    theme_classic()+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks=c(130,150,170))+
    theme(legend.key.width = unit(40, "pt")))

legd <- get_legend(plot_legend)


#### Organise figures ####

space <- ggplot()+ theme_void()


#Figure 2: Average fitness function and peaks
plot_row_bb <- grid.arrange(space,plot_av_bb,space,plot_th_mu_bb,space,plot_lw_lm_bb,space,plot_lo_ls_bb,ncol=8,widths=c(0.1,1,0.05,1,0.05,1,0.05,1))
plot_row_bp <- grid.arrange(space,plot_av_bp,space,plot_th_mu_bp,space,plot_lw_lm_bp,space,plot_lo_ls_bp,ncol=8,widths=c(0.1,1,0.05,1,0.05,1,0.05,1))
full_plot <- grid.arrange(plot_row_bb,plot_row_bp,nrow=2,heights=c(1,1.05))
grid.text(label="Optimum hatch date",x=(0.277), y=(0.52), rot=90)
grid.text(label="Max. fledging success (prob. or no.)",x=(0.53), y=(0.52), rot=90)
grid.text(label="Fitness function width",x=(0.77), y=(0.52), rot=90)
grid.text(label=expression(bold("Probability of success")),x=(0.025), y=(0.77), rot=90, gp=gpar(fontsize=14))
grid.text(label=expression(bold("Number fledged")),x=(0.025), y=(0.28), rot=90, gp=gpar(fontsize=14))
grid.text(label="a",x=(0.095), y=(0.98), rot=0, gp=gpar(fontsize=10))
grid.text(label="b",x=(0.34), y=(0.98), rot=0, gp=gpar(fontsize=10))
grid.text(label="c",x=(0.585), y=(0.98), rot=0, gp=gpar(fontsize=10))
grid.text(label="d",x=(0.83), y=(0.98), rot=0, gp=gpar(fontsize=10))
grid.text(label="e",x=(0.095), y=(0.49), rot=0, gp=gpar(fontsize=10))
grid.text(label="f",x=(0.34), y=(0.49), rot=0, gp=gpar(fontsize=10))
grid.text(label="g",x=(0.585), y=(0.49), rot=0, gp=gpar(fontsize=10))
grid.text(label="h",x=(0.83), y=(0.49), rot=0, gp=gpar(fontsize=10))
#5.5*9.5. 


# Supplementary figure: slopes of fitness width by caterpillar height
supp_lo_lm <- grid.arrange(plot_lo_lm_bb,plot_lo_lm_bp)
grid.text(label="Width of Fitness Function",x=(0.03), y=(0.52), rot=90) 
grid.text(label="a",x=(0.22), y=(0.98), rot=0, gp=gpar(fontsize=10))
grid.text(label="b",x=(0.22), y=(0.48), rot=0, gp=gpar(fontsize=10)) #5.5*3

# Supplementary figure: slopes of max fitness  by caterpillar height zoomed in
supp_lw_lm <- grid.arrange(plot_lw_lm_bb_zoom,plot_lw_lm_bp_zoom)
grid.text(label="a",x=(0.22), y=(0.98), rot=0, gp=gpar(fontsize=10))
grid.text(label="b",x=(0.22), y=(0.48), rot=0, gp=gpar(fontsize=10)) #5.5*3

# Supplementary figure: fitness function by caterpilla peak timing height and width
plot_col_bb <- grid.arrange(plot_mu_bb,plot_lm_bb,plot_ls_bb,space,nrow=4,heights=c(1,1,1,0.1))
plot_col_bp <- grid.arrange(plot_mu_bp,plot_lm_bp,plot_ls_bp,space,nrow=4,heights=c(1,1,1,0.1))
full_plot <- grid.arrange(space,plot_col_bb,space,plot_col_bp,legd,ncol=5,widths=c(0.05,1,0.15,1,0.7))
grid.text(label="Probability of success (min. 1 fledge)",x=(0.02), y=(0.5), rot=90)
grid.text(label="No. fledged (cond. min. 1 fledge)",x=(0.42), y=(0.5), rot=90)
grid.text(label="Hatch date (1 = 1 Jan)",x=(0.43), y=(0.03), rot=0)
grid.text(label="a",x=(0.105), y=(0.985), rot=0, gp=gpar(fontsize=10))
grid.text(label="b",x=(0.48), y=(0.985), rot=0, gp=gpar(fontsize=10))
grid.text(label="c",x=(0.105), y=(0.66), rot=0, gp=gpar(fontsize=10))
grid.text(label="d",x=(0.48), y=(0.66), rot=0, gp=gpar(fontsize=10))
grid.text(label="e",x=(0.105), y=(0.34), rot=0, gp=gpar(fontsize=10))
grid.text(label="f",x=(0.48), y=(0.34), rot=0, gp=gpar(fontsize=10))
#8*8


### Optimum relative to caterpillars for early and late caterpillar peaks
mu_0.1 <- 140-8-146
mu_0.9 <- 170-8-146
mu_0.1_uncent <- 140
mu_0.9_uncent <- 170

thta_0.1_bb <- th_bb+mu_th_bb*mu_0.1+146
thta_0.1_bp <- th_bp+mu_th_bp*mu_0.1+146
thta_0.9_bb <- th_bb+mu_th_bb*mu_0.9+146
thta_0.9_bp <- th_bp+mu_th_bp*mu_0.9+146

dif_0.1_bb <- thta_0.1_bb - (mu_0.1_uncent)
dif_0.1_bp <- thta_0.1_bp - (mu_0.1_uncent)
dif_0.9_bb <- thta_0.9_bb - (mu_0.9_uncent) 
dif_0.9_bp <- thta_0.9_bp - (mu_0.9_uncent) 

par(mfrow=c(2,2))
hist(dif_0.1_bb[,1],100,xlim=c(-50,0))  
abline(v=0,col=2)
hist(dif_0.1_bp[,1],100,xlim=c(-50,0))  
abline(v=0,col=2)
hist(dif_0.9_bb[,1],100,xlim=c(-50,0))  
abline(v=0,col=2)
hist(dif_0.9_bp[,1],100,xlim=c(-50,0))  
abline(v=0,col=2)

mean(dif_0.1_bb[,1])
posterior_interval(as.matrix(dif_0.1_bb), prob=0.95)

hist(dif_0.1_bb[,1]-dif_0.9_bb[,1],100)
hist(dif_0.1_bp[,1]-dif_0.9_bp[,1],100)
