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


#### Model SY mean hatch dates ####

hd_mod <- brm(hd_cent~1+(1|site)+(1|year)+(1|siteyear)+(1|fem), data=SYBhd, warmup = 1000, iter=12000, thin=10, chains=4)
plot(hd_mod)
summary(hd_mod)
save(hd_mod, file="hd_mod")
#load(file="hd_mod")

hd_ef <- posterior_samples(hd_mod,pars=c("year","site"))
s_hd_ef <- hd_ef[,c(384:427)]
y_hd_ef <- hd_ef[,c(372:382)]
sy_hd_ef <- hd_ef[,c(3:371)]
int_post <- posterior_samples(hd_mod,pars=c("b_Intercept"))

sy_hd_post <- data.frame(matrix(NA,ncol=ncol(sy_hd_ef),nrow=nrow(sy_hd_ef))) 

for(i in 1:ncol(sy_hd_post)){
  sy_hd_post[,i] <- int_post[,1]+sy_hd_ef[,i]+y_hd_ef[,SY_id$year_id[i]]+s_hd_ef[,SY_id$site_id[i]]
}


#### Model SY optimum hatch dates ####
opt_bb <-  data.frame(matrix(NA,ncol=nrow(SY_id),nrow=nrow(bb_sy))) 
opt_bp <-  data.frame(matrix(NA,ncol=nrow(SY_id_p),nrow=nrow(bp_sy))) 

for(i in 1:nrow(SY_id)){
  opt_bb[,i] <- th_bb[,1]+mu_th_bb[,1]*cater_mu_bb[,SY_id$siteyear_id[i]]+bb_s[,SY_id$site_id[i]]+bb_y[,SY_id$year_id[i]]+bb_sy[,SY_id$siteyear_id[i]]
}
for(i in 1:nrow(SY_id_p)){
  opt_bp[,i] <- th_bp[,1]+mu_th_bp[,1]*cater_mu_bp[,SY_id_p$siteyear_id[i]]+bp_s[,SY_id_p$site_id[i]]+bp_y[,SY_id_p$year_id[i]]+bp_sy[,SY_id_p$siteyear_id[i]]
}

#### Lag for each SY ####
lag_bb <- as.data.frame(as.matrix(sy_hd_post)-as.matrix(opt_bb))
lag_bp <- as.data.frame(as.matrix(sy_hd_post[,c(SY_id_p$siteyear_id_c)])-as.matrix(opt_bp))

#### Lag for each nest ####
SYBhd$lag_bb <- SYBhd$hd_cent-colMeans(opt_bb)[pmatch(SYBhd$siteyear,SY_id$siteyear,duplicates.ok = T)]
SYBhd_p$lag_bp <- SYBhd_p$hd_cent-colMeans(opt_bp)[pmatch(SYBhd_p$siteyear,SY_id_p$siteyear,duplicates.ok = T)]

## Lag histograms
xlims <- c(-11,33)
ylims <- c(1.5,45)
(plot_lag_sy_bb <- ggplot()+
    geom_histogram(aes(colMeans(lag_bb)),fill="grey",col="black", binwidth=1)+
    coord_cartesian(xlim = xlims,ylim=ylims)+
    xlab("Lag (site-year)")+
    ylab("Frequency")+
    theme_classic())

(plot_lag_sy_bp <- ggplot()+
    geom_histogram(aes(colMeans(lag_bp)),fill="grey",col="black", binwidth=1)+
    coord_cartesian(xlim = xlims,ylim=ylims)+
    xlab("Lag (site-year)")+
    ylab("")+
    theme_classic())

ylims <- c(1.5,105)

(plot_lag_n_bb <- ggplot()+
    geom_histogram(data=SYBhd,aes(lag_bb),fill="grey",col="black", binwidth=1)+
    coord_cartesian(xlim = xlims,ylim=ylims)+
    xlab("Lag (nest)")+
    ylab("Frequency")+
    theme_classic())

(plot_lag_n_bp <- ggplot()+
    geom_histogram(data=SYBhd_p,aes(lag_bp),fill="grey",col="black", binwidth=1)+
    coord_cartesian(xlim = xlims,ylim=ylims)+
    xlab("Lag (nest)")+
    ylab("")+
    theme_classic())


#### Cater and mean lag ####

SY_id$mu <- colMeans(cater_mu_bb)+8+146
SY_id_p$mu <- colMeans(cater_mu_bp)+8+146
SY_id$opt_bb <- colMeans(opt_bb)+146
SY_id_p$opt_bp <- colMeans(opt_bp)+146
SY_id$lag_bb <- colMeans(lag_bb)
SY_id_p$lag_bp <- colMeans(lag_bp)

(plot_lag_cater_bb <- ggplot(SY_id,aes(mu,lag_bb,col=year))+
    geom_point()+
    ylab("Lag (site-year)")+
    xlab("Caterpillar mean timing")+
    theme_classic()+
    ylim(-5.5,22)+
    xlim(135,180)+
    labs(col = "Year")+ theme(legend.position = "none"))

(plot_lag_cater_bp <- ggplot(SY_id_p,aes(mu,lag_bp,col=as.factor(year)))+
    geom_point()+
    ylab("")+
    xlab("Caterpillar mean timing")+
    theme_classic()+
    ylim(-5.5,22)+
    xlim(135,180)+
    labs(col = "Year")+ 
    scale_colour_discrete(labels=c("2014","","2016","","2018","","2020","","2022","","2024"))+
    theme(legend.key.height = unit(0, "pt"),
          legend.key.size = unit(0, 'pt'),
          legend.position.inside = c(0.95,0.72),
          legend.text = element_text(size=7))+
    guides(color = guide_legend(
      override.aes=list(shape = 15,size=1.8),
      position = "inside")))

### Arrange figure 3
space <- ggplot()+theme_void()
F3_col1 <- grid.arrange(space, plot_lag_sy_bb,plot_lag_n_bb, plot_lag_cater_bb,nrow=4,heights=c(0.2,1,1,1))
F3_col2 <- grid.arrange(space, plot_lag_sy_bp,plot_lag_n_bp, plot_lag_cater_bp,nrow=4,heights=c(0.2,1,1,1))

F3 <- grid.arrange(F3_col1,F3_col2, ncol=2)
grid.text(label=expression(bold("Probability of success")),x=(0.27), y=(0.96), rot=0, gp=gpar(fontsize=14))
grid.text(label=expression(bold("Number fledged")),x=(0.77), y=(0.96), rot=0, gp=gpar(fontsize=14))
grid.text(label="a",x=(0.1), y=(0.92), rot=0, gp=gpar(fontsize=10))
grid.text(label="b",x=(0.6), y=(0.92), rot=0, gp=gpar(fontsize=10))
grid.text(label="c",x=(0.1), y=(0.6), rot=0, gp=gpar(fontsize=10))
grid.text(label="d",x=(0.6), y=(0.6), rot=0, gp=gpar(fontsize=10))
grid.text(label="e",x=(0.1), y=(0.29), rot=0, gp=gpar(fontsize=10))
grid.text(label="f",x=(0.6), y=(0.29), rot=0, gp=gpar(fontsize=10))
#8.5*6.8
