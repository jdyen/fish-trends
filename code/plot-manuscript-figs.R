# set working directory
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")

# load and plot fitted models
mod_list <- dir("./outputs/fitted/")

# set up a summary output table and predictor names for labels
cov_names <- c("Mean pre-spawning flow",
               "CV of annual flow", 
               "CV of pre-spawning flow")

# choose combos of plots to include in panels
#
# plot_fitted(mod$mod_sum,
#             mod$bugsdata,
#             mod$sp_names,
#             mod$spp, mod$resp)
# plot_trend(mod$mod_sum,
#            mod$bugsdata,
#            mod$sp_names,
#            mod$spp, mod$resp)
# plot_covars(mod$mod_sum,
#             mod$bugsdata,
#             resp = mod$resp,
#             cov_names = cov_names,
#             covar_std = mod$covar_std)
