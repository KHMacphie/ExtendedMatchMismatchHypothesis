## Fitness Model Table ##

# Load model
mod <- readRDS("FullMMH_mod.rds")
rstan::get_num_divergent(mod)

sum_names <- data.frame(name=rownames(summary(mod)$summary))
keepeffs <- c("theta_bb","logomega_bb","logWmax_bb",
              "mu_theta_bb","ls_logomega_bb","lm_logomega_bb","lm_logWmax_bb",
              "sd_site_bb[1]","sd_site_bb[2]","sd_site_bb[3]",
              "sd_year_bb[1]","sd_year_bb[2]","sd_year_bb[3]",
              "sd_siteyear_bb[1]","sd_siteyear_bb[2]","sd_siteyear_bb[3]",
              "sd_fem_bb",
              "omega_bb[1,1]","omega_bb[2,1]","omega_bb[3,1]","omega_bb[1,2]","omega_bb[2,2]","omega_bb[3,2]","omega_bb[1,3]","omega_bb[2,3]","omega_bb[3,3]",
              "theta_bp","logomega_bp","logWmax_bp",
              "mu_theta_bp","ls_logomega_bp","lm_logomega_bp","lm_logWmax_bp",
              "sd_site_bp[1]","sd_site_bp[2]","sd_site_bp[3]",
              "sd_year_bp[1]","sd_year_bp[2]","sd_year_bp[3]",
              "sd_siteyear_bp[1]","sd_siteyear_bp[2]","sd_siteyear_bp[3]",
              "sd_fem_bp", "lambda",
              "omega_bp[1,1]","omega_bp[2,1]","omega_bp[3,1]","omega_bp[1,2]","omega_bp[2,2]","omega_bp[3,2]","omega_bp[1,3]","omega_bp[2,3]","omega_bp[3,3]",
              "bbp_reg[1]","bbp_reg[2]",
              "mu","logsigma","logmax",
              "sd_site_c[1]","sd_site_c[2]","sd_site_c[3]",
              "sd_year_c[1]","sd_year_c[2]","sd_year_c[3]",
              "sd_siteyear_c[1]","sd_siteyear_c[2]","sd_siteyear_c[3]",
              "sd_siteday","sd_treeID","sd_rec","sd_resid_c",
              "omega_c[1,1]","omega_c[2,1]","omega_c[3,1]","omega_c[1,2]","omega_c[2,2]","omega_c[3,2]","omega_c[1,3]","omega_c[2,3]","omega_c[3,3]")

keepfixed <- as.data.frame(summary(mod)$summary[keepeffs,c(1,4,6,8:10)])
keepfixed <- round(keepfixed,3)
colnames(keepfixed) <- c("Posterior Mean","2.5% CI","Posterior Median","97.5% CI","Effective Sample Size","Rhat") 
keepfixed <- keepfixed[,c(1,3,2,4,5,6)]
keepfixed$`Effective Sample Size` <- round(keepfixed$`Effective Sample Size`,0)
keepfixed$`Prior Standard Deviation or Scale` <- 
  as.character(c(20,5,sqrt(pi^2/3),5,5,5,5,10,3,3,10,3,3,10,3,3,10,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                 20,5,10,5,5,5,5,10,3,3,10,3,3,10,3,3,10,10,NA,NA,NA,NA,NA,NA,NA,NA,NA,10,10,
                 20,5,10,10,10,10,10,10,10,10,10,10,10,10,10,10,NA,NA,NA,NA,NA,NA,NA,NA,NA))

write.csv(keepfixed,"FitnessModelTable.csv")
