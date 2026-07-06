rm(list=ls())
setwd("")
source("OrganisingDataframes.R")
source("ResultsFunctions.R")

library(rstan)
library(extraDistr)
library(gp)


# Load model
mod <- readRDS("FullMMH_estmissing_mod.rds")


## Structure
#Bird Bernoulli Gaussian parameter calculations
#theta_bb_obs = theta_bb+mu_theta_bb*c_mu_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,1]+year_effs_bb[year_id_bb,1]+siteyear_effs_bb[siteyear_id_bb,1];
#omega_bb_obs = exp(logomega_bb+ls_logomega_bb*c_ls_bb[siteyear_id_bb]+lm_logomega_bb*c_lm_bb[siteyear_id_bb]+site_effs_bb[site_id_bb,2]+year_effs_bb[year_id_bb,2]+siteyear_effs_bb[siteyear_id_bb,2]);
#logWmax_bb_obs = logWmax_bb + lm_logWmax_bb*c_lm_bb[siteyear_id_bb] + site_effs_bb[site_id_bb,3] + year_effs_bb[year_id_bb,3] + siteyear_effs_bb[siteyear_id_bb,3] + fem_effs_bb[fem_id_bb];

# theta_bb
# mu_theta_bb
# c_mu_bb[siteyear_id_bb]
# site_effs_bb[site_id_bb,1:3]
# year_effs_bb[year_id_bb,1:3]
# siteyear_effs_bb[siteyear_id_bb,1:3];
# logomega_bb
# ls_logomega_bb
# c_ls_bb[siteyear_id_bb]
# lm_logomega_bb
# logWmax_bb
# lm_logWmax_bb
# c_lm_bb[siteyear_id_bb]
# fem_effs_bb[fem_id_bb]


#Bird truncated-Poisson Gaussian parameter calculations
#theta_bp_obs = theta_bp+mu_theta_bp*c_mu_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,1]+year_effs_bp[year_id_bp,1]+siteyear_effs_bp[siteyear_id_bp,1];
#omega_bp_obs = exp(logomega_bp+ls_logomega_bp*c_ls_bp[siteyear_id_bp]+lm_logomega_bp*c_lm_bp[siteyear_id_bp]+site_effs_bp[site_id_bp,2]+year_effs_bp[year_id_bp,2]+siteyear_effs_bp[siteyear_id_bp,2]);
#logWmax_bp_obs = logWmax_bp + lm_logWmax_bp*c_lm_bp[siteyear_id_bp] + site_effs_bp[site_id_bp,3] + year_effs_bp[year_id_bp,3] + siteyear_effs_bp[siteyear_id_bp,3] + fem_effs_bp[fem_id_bp];
#
# theta_bp
# mu_theta_bp
# c_mu_bp[siteyear_id_bp]
# site_effs_bp[site_id_bp,1:3]
# year_effs_bp[year_id_bp,1:3]
# siteyear_effs_bp[siteyear_id_bp,1:3]
# logomega_bp
# ls_logomega_bp
# c_ls_bp[siteyear_id_bp]
# lm_logomega_bp
# logWmax_bp
# lm_logWmax_bp
# c_lm_bp[siteyear_id_bp]
# fem_effs_bp[fem_id_bp]
# 
# siteyear_missing_effs[bb_missing_in_bp,1:3]

#y_1bb = logWmax_bb_obs - square(date_bb-theta_bb_obs) ./ (2 * square(omega_bb_obs));
#y_bb ~ bernoulli_logit(y_1bb);

#y_1bp = exp(logWmax_bp_obs - square(date_bp-theta_bp_obs) ./ (2 * square(omega_bp_obs)));
#tgenpoiss_lpmf(y_bp[i] | y_1bp[i], lambda)


## Posteriors
cater_mu <- stanpost(model=mod, parameters=c("c_mu_bb"))
cater_ls <- stanpost(model=mod, parameters=c("c_ls_bb"))
cater_lm <- stanpost(model=mod, parameters=c("c_lm_bb"))

th_bb <- stanpost(model=mod, parameters=c("theta_bb"))
lo_bb <- stanpost(model=mod, parameters=c("logomega_bb"))
lw_bb <- stanpost(model=mod, parameters=c("logWmax_bb"))

mu_th_bb <- stanpost(model=mod, parameters=c("mu_theta_bb"))
ls_lo_bb <- stanpost(model=mod, parameters=c("ls_logomega_bb"))
lm_lo_bb<- stanpost(model=mod, parameters=c("lm_logomega_bb"))
lm_lw_bb<- stanpost(model=mod, parameters=c("lm_logWmax_bb"))

bb_fem_sd <- stanpost(model=mod, parameters=c("sd_fem_bb"))
bb_fem_scl <- stanpost(model=mod, parameters=c("fem_scaled_bb"))

th_bp <- stanpost(model=mod, parameters=c("theta_bp"))
lo_bp <- stanpost(model=mod, parameters=c("logomega_bp"))
lw_bp <- stanpost(model=mod, parameters=c("logWmax_bp"))

mu_th_bp <- stanpost(model=mod, parameters=c("mu_theta_bp"))
ls_lo_bp <- stanpost(model=mod, parameters=c("ls_logomega_bp"))
lm_lo_bp <- stanpost(model=mod, parameters=c("lm_logomega_bp"))
lm_lw_bp <- stanpost(model=mod, parameters=c("lm_logWmax_bp"))

bp_fem_sd <- stanpost(model=mod, parameters=c("sd_fem_bp"))
bp_fem_scl <- stanpost(model=mod, parameters=c("fem_scaled_bp"))

bb_s <- stanpost(model=mod, parameters=c("site_effs_bb"))
bb_y <- stanpost(model=mod, parameters=c("year_effs_bb"))
bb_sy <- stanpost(model=mod, parameters=c("siteyear_effs_bb"))

bp_s <- stanpost(model=mod, parameters=c("site_effs_bp"))
bp_y <- stanpost(model=mod, parameters=c("year_effs_bp"))
bp_sy <- stanpost(model=mod, parameters=c("siteyear_effs_bp"))

missing <- stanpost(model=mod, parameters=c("siteyear_missing_effs"))
missingrows <- as.numeric(setdiff(SY_id_p$siteyear_id,SY_id_p$siteyear_id_c))

lambda <- stanpost(model=mod, parameters=c("lambda"))

# Lining up missing SY in bp
SY_id$siteyear_id_bp <- SY_id_p$siteyear_id[pmatch(SY_id$siteyear,SY_id_p$siteyear)]
SY_id$siteyear_bp_01 <- ifelse(is.na(SY_id$siteyear_id_bp),0,1) 

for(i in 1:nrow(SY_id)){
  if(SY_id$siteyear_bp_01[i]==0){
    SY_id$siteyear_id_bp[i] <- which(missingrows==i)+359
  } else {
    SY_id$siteyear_id_bp[i] <- SY_id$siteyear_id_bp[i]
  }
}

# ordering posterior columns so missing theta are with main effects theta etc
bp_sy_incmissing <- cbind(bp_sy[,1:359],missing[,1:10],bp_sy[,360:718],missing[,11:20],bp_sy[,719:1077],missing[21:30])

# Idenitfy females without BP estimate
fembp <- data.frame(fem=levels(as.factor(SYBhd_p$fem)),
                    femid=as.numeric(as.factor(levels(as.factor(SYBhd_p$fem)))))

SYBhd$femid_bb <- as.numeric(as.factor(SYBhd$fem))
SYBhd$femid_bp <- fembp$femid[pmatch(SYBhd$fem,fembp$fem,duplicates.ok = T)]
SYBhd$femid_bp <- ifelse(is.na(SYBhd$femid_bp),881,SYBhd$femid_bp) #all females with out BP as 881
bp_fem_scl$no_bp <- 0 # 881 'female' effect size as 0 - no effect estimated, assume average

# reduce to relevant data
sim_df <- data.frame(year=SYBhd$year,
                     year_id=as.numeric(as.factor(SYBhd$year)),
                     site=SYBhd$site,
                     site_id=as.numeric(as.factor(SYBhd$site)),
                     siteyear=SYBhd$siteyear,
                     siteyear_id=as.numeric(as.factor(SYBhd$siteyear)),
                     siteyear_id_bp=SY_id$siteyear_id_bp[pmatch(SYBhd$siteyear,SY_id$siteyear,duplicates.ok = T)],
                     fem_id_bb=SYBhd$femid_bb,  
                     fem_id_bp=SYBhd$femid_bp,  
                     date=SYBhd$hd_cent,
                     y_bb=SYBhd$fledged.binom,
                     y_bp=SYBhd$fledged)

# set up data frames - run through every iteration
iter <- 1:nrow(th_bb)
sim_dat_bb <- data.frame(matrix(NA,nrow=nrow(sim_df),ncol=length(iter)))
sim_dat_bp <- data.frame(matrix(NA,nrow=nrow(sim_df),ncol=length(iter)))

pb <- txtProgressBar(min = 1, max = length(iter), style = 3)

for(i in 1:length(iter)){ # for a given iteration
  sim_df$th_bb <- th_bb[iter[i],1] + 
    mu_th_bb[iter[i],1]*as.numeric(cater_mu[iter[i],c(sim_df$siteyear_id)])+
    as.numeric(bb_s[iter[i],c(sim_df$site_id)])+
    as.numeric(bb_y[iter[i],c(sim_df$year_id)])+
    as.numeric(bb_sy[iter[i],c(sim_df$siteyear_id)]) # predicted theta bb
  
  sim_df$lo_bb <- lo_bb[iter[i],1] + 
    ls_lo_bb[iter[i],1]*as.numeric(cater_ls[iter[i],c(sim_df$siteyear_id)])+
    lm_lo_bb[iter[i],1]*as.numeric(cater_lm[iter[i],c(sim_df$siteyear_id)])+
    as.numeric(bb_s[iter[i],c(sim_df$site_id+44)])+
    as.numeric(bb_y[iter[i],c(sim_df$year_id+11)])+
    as.numeric(bb_sy[iter[i],c(sim_df$siteyear_id+369)]) # predicted logomega bb
  
  sim_df$lw_bb <- lw_bb[iter[i],1] + 
    lm_lw_bb[iter[i],1]*as.numeric(cater_lm[iter[i],c(sim_df$siteyear_id)])+
    as.numeric(bb_s[iter[i],c(sim_df$site_id+44+44)])+
    as.numeric(bb_y[iter[i],c(sim_df$year_id+11+11)])+
    as.numeric(bb_sy[iter[i],c(sim_df$siteyear_id+369+369)])+
    bb_fem_sd[iter[1],1]*as.numeric(bb_fem_scl[iter[i],c(sim_df$fem_id_bb)]) #predicted logWmax bb
  
  
  sim_df$th_bp <- th_bp[iter[i],1] + 
    mu_th_bp[iter[i],1]*as.numeric(cater_mu[iter[i],c(sim_df$siteyear_id)])+
    as.numeric(bp_s[iter[i],c(sim_df$site_id)])+
    as.numeric(bp_y[iter[i],c(sim_df$year_id)])+
    as.numeric(bp_sy_incmissing[iter[i],c(sim_df$siteyear_id_bp)]) # predicted theta bp
  
  sim_df$lo_bp <- lo_bp[iter[i],1] + 
    ls_lo_bp[iter[i],1]*as.numeric(cater_ls[iter[i],c(sim_df$siteyear_id)])+
    lm_lo_bp[iter[i],1]*as.numeric(cater_lm[iter[i],c(sim_df$siteyear_id)])+
    as.numeric(bp_s[iter[i],c(sim_df$site_id+44)])+
    as.numeric(bp_y[iter[i],c(sim_df$year_id+11)])+
    as.numeric(bp_sy_incmissing[iter[i],c(sim_df$siteyear_id_bp+369)]) #predicted logomega bp
  
  sim_df$lw_bp <- lw_bp[iter[i],1] + 
    lm_lw_bp[iter[i],1]*as.numeric(cater_lm[iter[i],c(sim_df$siteyear_id)])+
    as.numeric(bp_s[iter[i],c(sim_df$site_id+44+44)])+
    as.numeric(bp_y[iter[i],c(sim_df$year_id+11+11)])+
    as.numeric(bp_sy_incmissing[iter[i],c(sim_df$siteyear_id_bp+369+369)])+
    bp_fem_sd[iter[1],1]*as.numeric(bp_fem_scl[iter[i],c(sim_df$fem_id_bp)]) # predicted logWmax bp
  
  sim_df$mean_bb <- plogis(sim_df$lw_bb - (sim_df$date-sim_df$th_bb)^2 / (2 * exp(sim_df$lo_bb)^2)) # predicted probability of success
  sim_df$mean_bp <- exp(sim_df$lw_bp - (sim_df$date-sim_df$th_bp)^2 / (2 * exp(sim_df$lo_bp)^2)) # predicted number fledged
  
  sim_dat_bb[,i] <- rbinom(nrow(sim_df),1,prob=sim_df$mean_bb) # random draw from binomial dist
  for(j in 1:nrow(sim_df)){
    sim_dat_bp[j,i] <- tryCatch(rgp(5,theta=sim_df$mean_bp[j],lambda=lambda[i,1],method="Inversion")[1],error=function(e){NA})
  } # random draw from generalised poisson dist 
  
  sim_df[,c(13:ncol(sim_df))] <- NULL
  
  setTxtProgressBar(pb, i)
}

close(pb)

#Count NAs from estimates not fitting generalise poisson generator
store <- c()
for(k in 1:ncol(sim_dat_bp)){
  store[k] <- length(which(is.na(sim_dat_bp[,k])==T)) 
}
rm(k)
hist(store)
(mean(store)/nrow(sim_dat_bp))*100 # mean % NA

# dataframe for estimates combining bb and bp
sim_dat_pred <- data.frame(matrix(NA,nrow=nrow(sim_df),ncol=length(iter)))

pb <- txtProgressBar(min = 1, max = ncol(sim_dat_bb), style = 3)

for(l in 1:ncol(sim_dat_bb)){
  for(k in 1:nrow(sim_dat_bb)){
    if(sim_dat_bb[k,l]==0){ # if bb is 0 then zero
      sim_dat_pred[k,l] <- 0
    } else if(is.na(sim_dat_bp[k,l])){ # if bp is NA then NA
      sim_dat_pred[k,l] <- NA
    } else if(sim_dat_bb[k,l]==1&sim_dat_bp[k,l]>0){ # if bb is 1 and bp is >0 then bp
      sim_dat_pred[k,l] <- sim_dat_bp[k,l]
    } else if(sim_dat_bb[k,l]==1&sim_dat_bp[k,l]==0){ # if bb is 1 and bp is 0 then NA (truncated)
      sim_dat_pred[k,l] <- NA
    } else {
      sim_dat_pred[k,l] <- NA
    }
  }
  setTxtProgressBar(pb, l)
}

close(pb)

# recalculate NAs
store2 <- c()
for(m in 1:ncol(sim_dat_bp)){
  store2[m] <- length(which(is.na(sim_dat_pred[,m])==T))
}
rm(m)
hist(store2)
(mean(store2)/nrow(sim_dat_bp))*100 # mean % NA

# Save simulation
#write.csv(sim_dat_pred,"simulation_periter.csv",row.names=F)
#sim_dat_pred <- read.csv("simulation_periter.csv")

#Compare total fledged in the dataset to total predicted
hist(colSums(sim_dat_pred,na.rm=T),100)
abline(v=sum(SYBhd$fledged),col=2)


# Compare proportion zero in prediction to in data
propzero <- c()
for(i in 1:ncol(sim_dat_pred)){
  propzero[i] <- length(which(sim_dat_pred[,i]==0))/nrow(sim_dat_pred)
}

hist(propzero,100)
abline(v=length(which(SYBhd$fledged.binom==0))/nrow(SYBhd),col=2)

library(ggplot2)
(fledge <- ggplot()+
    geom_histogram(aes(colSums(sim_dat_pred,na.rm=T)),fill="lightgrey",col=1,linewidth = 0.5)+
    geom_vline(aes(xintercept=sum(SYBhd$fledged)),col=2,linetype="longdash",linewidth=1.5)+
    xlab("Total fledged")+
    ylab("Frequency")+
    theme_classic())
(fledgevar <- ggplot()+
    geom_histogram(aes(apply(sim_dat_pred,2,function(x) var(x,na.rm = T))),fill="lightgrey",col=1,linewidth = 0.5)+
    geom_vline(aes(xintercept=var(SYBhd$fledged)),col=2,linetype="longdash",linewidth=1.5)+
    xlab("Variance in no. fledged")+
    ylab("Frequency")+
    theme_classic())
(zeros <- ggplot()+
    geom_histogram(aes(propzero),fill="lightgrey",col=1,linewidth = 0.5)+
    geom_vline(aes(xintercept=length(which(SYBhd$fledged.binom==0))/nrow(SYBhd)),col=2,linetype="longdash",linewidth=1.5)+
    xlab("Proportion failed")+
    ylab("Frequency")+
    theme_classic())
library(gridExtra)
library(grid)
space <- ggplot()+theme_void()
grid.arrange(fledge,space,fledgevar,space,zeros,nrow=5,heights=c(1,0.1,1,0.1,1))
grid.text("a",x=0.17,y=0.98)
grid.text("b",x=0.17,y=0.6)
grid.text("c",x=0.17,y=0.27) #8*5"