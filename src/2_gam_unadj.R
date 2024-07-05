rm(list=ls())

# Load libraries and configuration
source(here::here("src/0-config.R"))
library(here)
library(VIM)

# Construct the path to the RDS file directly within the project folder
#data_path <- here("bangladesh-cleaned-master-data.RDS")
data_path <- "C:/Users/andre/Documents/EE/eed-substudy-data/bangladesh-cleaned-master-data.RDS"

# Read the data from the RDS file
d <- readRDS(data_path)

# Define exposure (Xvars) and outcome (Yvars) variables
Xvars <- c("agp", "crp", "sumscore_t0_mom_Z")
Yvars <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")

selected_columns <- c(Xvars, Yvars)

# Hypothesis testing for these variables
H2_models <- NULL

# Loop over exposure-outcome pairs
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H2_models <- bind_rows(H2_models, res)
  }
}


H2_res <- NULL

for (i in 1:nrow(H2_models)) {
  res <- data.frame(X = H2_models$X[i], Y = H2_models$Y[i])
  
  # Get predictions of differences from the 25th percentile of exposure
  preds <- predict_gam_diff(fit = H2_models$fit[i][[1]], d = H2_models$dat[i][[1]], quantile_diff = c(0.25, 0.75), Xvar = res$X, Yvar = res$Y)
  
  # Store primary contrasts
  H2_res <- bind_rows(H2_res, preds$res)
}

H2_plot_list <- NULL
H2_plot_data <- NULL

for (i in 1:nrow(H2_models)) {
  res <- data.frame(X = H2_models$X[i], Y = H2_models$Y[i])
  
  # Fit spline with simultaneous confidence intervals
  simul_plot <- gam_simul_CI(H2_models$fit[i][[1]], H2_models$dat[i][[1]], xlab = res$X, ylab = res$Y, title = "")
  
  H2_plot_list[[i]] <- simul_plot$p
  H2_plot_data <- rbind(H2_plot_data, data.frame(Xvar = res$X, Yvar = res$Y, adj = 0, simul_plot$pred))
}

# Save models
saveRDS(H2_models, file = here("H2_models.RDS"))

# Save results
saveRDS(H2_res, here("H2_res.RDS"))

# Save plot data
saveRDS(H2_plot_data, here("H2_plot_data.RDS"))


