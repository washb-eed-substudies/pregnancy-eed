rm(list=ls())

# Load libraries and configuration
source(here::here("src/0-config.R"))
source(here::here("functions/cowboy_glm.R"))
library(here)
library(VIM)

# Construct the path to the RDS file directly within the project folder
#data_path <- here("bangladesh-cleaned-master-data.RDS")
data_path <- "C:/Users/andre/Documents/EE/eed-substudy-data/bangladesh-cleaned-master-data.RDS"

# Read the data from the RDS file
d <- readRDS(data_path)


#-------------------------------------------------------------------------------
# Hypothesis 1
#-------------------------------------------------------------------------------
Xvars <- c("ln_preg_cort")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", 
           "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", 
           "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

#Fit models

res_adj <- washb_glm_lasso_boot(d=d, X=Xvars, Y=Yvars[1],  W=NULL)


H1_unadj_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    
    res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j,  res_adj$TR)
    H1_unadj_res <- bind_rows(H1_unadj_res, res)
  }
}


# Save results
saveRDS(H1_unadj_res, here("results/unadjusted/H1_unadj_res.RDS"))


#-------------------------------------------------------------------------------
# Hypothesis 2
#-------------------------------------------------------------------------------
Xvars <- c("agp", "crp", "sumscore_t0_mom_Z")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", 
           "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", 
           "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

H2_unadj_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    
    res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j,  res_adj$TR)
    H2_unadj_res <- bind_rows(H2_unadj_res, res)
  }
}

# Save results
saveRDS(H2_unadj_res, here("results/unadjusted/H2_unadj_res.RDS"))

#-------------------------------------------------------------------------------
# Hypothesis 3
#-------------------------------------------------------------------------------
Xvars <- c("ln_preg_estri")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", 
           "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", 
           "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

H3_unadj_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    
    res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j,  res_adj$TR)
    H3_unadj_res <- bind_rows(H3_unadj_res, res)
  }
}

# Save results
saveRDS(H3_unadj_res, here("results/unadjusted/H3_unadj_res.RDS"))
