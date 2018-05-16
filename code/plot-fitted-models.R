# set working directory
setwd("~/Dropbox/research/fish-trends/")

# load packages
library(ggplot2)
library(gridExtra)
# library(scales)

# load and plot fitted models
mod_list <- dir("./outputs/fitted/")

for (i in seq_along(mod_list)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list[i])))
  # plot fitted trends
  pdf(paste0("./outputs/plots/", mod$spp, "_", mod$resp, 
             ifelse(!is.null(mod$mod_sum$cov_plot_vals), "_covar", ""), ".pdf"))
  plot_fitted(mod$mod_sum,
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
                cov_names = c("Mean pre-spawning flow",
                              "CV of annual flow", 
                              "CV of pre-spawning flow"))
    dev.off()
  }
  
}