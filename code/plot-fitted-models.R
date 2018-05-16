# set working directory
setwd("~/Dropbox/research/fish-trends/")

# load packages
library(ggplot2)
library(gridExtra)
# library(scales)

# load and plot fitted models
mod_list <- dir("./outputs/fitted/")

for (i in seq_along(mod_list)) {
  
  mod <- get(load(mod_list[i]))
  # plot fitted trends
  pdf(paste0("./outputs/plots/", spp, "_", resp, "_", save.name, ".pdf"))
  plot_fitted(mod$mod_sum,
              mod$bugsdata,
              mod$sp_names,
              mod$spp, resp)
  dev.off()
  
  # plot covariate associations
  if (!is.null(mod$mod_sum$cov_plot_vals)) {
    pdf(paste0("./outputs/plots/", spp, "_flow_effects_", resp, ".pdf"))
    plot_covars(mod$mod_sum,
                mod$bugsdata,
                resp = mod$resp,
                cov_names = c("Mean pre-spawning flow",
                              "CV of annual flow", 
                              "CV of pre-spawning flow"))
    dev.off()
  }
  
}