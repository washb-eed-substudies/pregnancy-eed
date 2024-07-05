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

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu","gest_age_weeks",
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth_scaled",
         "tr", "life_viol_any_t3_cat", "viol_any_preg_cat")

#Add in time varying covariates:
Wvars1 <- c(Wvars, c("ageday_st1", "month_blood_t0", "month_bt1"))
Wvars2 <- c(Wvars, c("ageday_st2", "month_blood_t0", "month_bt2"))
Wvars3 <- c(Wvars, c("ageday_st3", "month_blood_t0", "month_bt3"))

#Add time of day for blood collection (for cortisol measurement)
d$time_of_day_cort<-as.numeric(hms(d$time_of_day_cort))

pick_covariates <- function(j){
  # j is outcome as string
  # choose correct adjustment set based on outcome
  if(grepl("t1", j)){Wset = Wvars1}
  if(grepl("t2", j)){Wset = Wvars2}
  if(grepl("t3", j)){Wset = Wvars3}
  else{Wset = Wvars3}
  
  if(grepl("cort", j)){Wset = c(Wset, "time_of_day_cort")}
  return(Wset)
}


#-------------------------------------------------------------------------------
# Hypothesis 1
#-------------------------------------------------------------------------------
Xvars <- c("ln_preg_cort")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", 
           "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", 
           "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

#Fit models

i=Xvars
j=Yvars[1]
res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=Wvars)

H1_adj_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j,  res_adj$TR)
    H1_adj_res <- bind_rows(H1_adj_res, res)
  }
}


# Save results
saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))


#-------------------------------------------------------------------------------
# Hypothesis 2
#-------------------------------------------------------------------------------
Xvars <- c("agp", "crp", "sumscore_t0_mom_Z")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", 
           "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", 
           "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

H2_adj_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j,  res_adj$TR)
    H2_adj_res <- bind_rows(H2_adj_res, res)
  }
}

# Save results
saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))

#-------------------------------------------------------------------------------
# Hypothesis 3
#-------------------------------------------------------------------------------
Xvars <- c("ln_preg_estri")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", 
           "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", 
           "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

H3_adj_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- washb_glm_lasso_boot(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j,  res_adj$TR)
    H3_adj_res <- bind_rows(H3_adj_res, res)
  }
}

# Save results
saveRDS(H3_adj_res, here("results/adjusted/H3_adj_res.RDS"))
