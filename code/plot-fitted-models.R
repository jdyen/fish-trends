# set working directory
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")

# load and plot fitted models
mod_list <- dir("./outputs/fitted/")

# set up a summary output table and predictor names for labels
cov_names <- c("Mean annual flow",
               "CV of annual flow", 
               "CV of pre-spawning flow")
summary_table <- matrix(NA, nrow = length(mod_list), ncol = 8)
rownames(summary_table) <- seq_len(length(mod_list))
for (i in seq_along(mod_list)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list[i])))
  # plot fitted trends
  pdf(paste0("./outputs/plots/fitted_", mod$spp, "_", mod$resp, 
             ifelse(!is.null(mod$mod_sum$cov_plot_vals), "_covar", ""), ".pdf"))
  plot_fitted(mod$mod_sum,
              mod$bugsdata,
              mod$sp_names,
              mod$spp, mod$resp)
  dev.off()

  # plot relative trends
  pdf(paste0("./outputs/plots/trends_", mod$spp, "_", mod$resp, 
             ifelse(!is.null(mod$mod_sum$cov_plot_vals), "_covar", ""), ".pdf"))
  plot_trend(mod$mod_sum,
             mod$bugsdata,
             mod$sp_names,
             mod$spp, mod$resp)
  dev.off()
  
  # plot covariate associations
  if (!is.null(mod$mod_sum$cov_plot_vals)) {
    pdf(paste0("./outputs/plots/", mod$spp, "_flow_effects_", mod$resp, ".pdf"))
    plot_covars(mod$mod_sum,
                mod$bugsdata,
                resp = mod$resp,
                cov_names = cov_names,
                covar_std = mod$covar_std)
    dev.off()
  }
  
  if (!is.na(mod$mod_sum$cov_inc[1])) {
    summary_table[i, ] <- round(c(mod$mod_sum$r2, mean(mod$mod_sum$ppps), mod$mod_sum$cov_inc, mod$mod_sum$cov_or), 2)
  } else {
    summary_table[i, ] <- round(c(mod$mod_sum$r2, mean(mod$mod_sum$ppps), rep(NA, 6)), 2)
  }
  rownames(summary_table)[i] <- substr(mod_list[i], 1, nchar(mod_list[i]) - 6)
  
}
colnames(summary_table) <- c("r2", "PPP", rep(cov_names, 2))
