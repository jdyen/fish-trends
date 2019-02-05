# set working directory
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")
source("./code/plot-helpers2.R")

# load and plot fitted models
mod_list_all <- dir("./outputs/fitted/")

# subset to the trend and covar_trend models for plotting
mod_list <- mod_list_all[-grep("covar.rds", mod_list_all)]
mod_list <- mod_list[-grep("int_re.rds", mod_list)]

# set up a summary output table and predictor names for labels
cov_names <- c("Mean daily flow (ML)",
               "CV of daily flows")
sp_name <- c('Murray cod', 'Trout cod', 'Golden perch',
             'Silver perch', 'Murray river rainbowfish',
             'Common carp')
sp_name_covs_vol <- c('Trout cod', 'Golden perch', 'Golden perch', 'Silver perch')
sp_name_covs_var <- c('Trout cod', rep('Silver perch', 3))

# trends plot
mod_list_abund <- mod_list[grep('abundance.rds', mod_list)][c(3, 6, 2, 5, 4, 1)]
mod_list_biom <- mod_list[grep('weight.rds', mod_list)][c(3, 6, 2, 5, 4, 1)]
mod_list_covs_vol <- mod_list_all[c(46, 10, 14, 34)]
vars_vol <- rep(1, 4)
mod_list_covs_var <- mod_list_all[c(42, 46, 34, 38)]
vars_var <- rep(2, 4)

# plot fitted abundance trends
pdf(file = "./outputs/plots/Fig1.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 4.5, 2.1, 1.1))
for (i in seq_along(mod_list_abund)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_abund[i])))
  
  plot_trend(mod$mod_sum,
             mod$bugsdata,
             mod$sp_names,
             mod$spp, mod$resp)
  
  mtext(sp_name[i], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot fitted biomass trends
pdf(file = "./outputs/plots/Fig2.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 4.5, 2.1, 1.1))
for (i in seq_along(mod_list_biom)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_biom[i])))
  
  plot_trend(mod$mod_sum,
             mod$bugsdata,
             mod$sp_names,
             mod$spp, mod$resp)
  
  mtext(sp_name[i], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot covariate effects (volume)
pdf(file = "./outputs/plots/Fig3.pdf", width = 7, height = 7)
par(mfrow = c(2, 2), mar = c(4.5, 4.3, 2.1, 1.1))
for (i in seq_along(mod_list_covs_vol)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_covs_vol[i])))
  
  plot_covars_single(mod$mod_sum,
                     mod$bugsdata,
                     resp = mod$resp,
                     cov_names = cov_names,
                     covar_std = mod$covar_std,
                     subset = vars_vol[i])
  
  mtext(sp_name_covs_vol[i], side = 3, adj = 1, line = 0, cex = 1)
  
}
dev.off()

# plot covariate effects (variability)
pdf(file = "./outputs/plots/Fig4.pdf", width = 7, height = 7)
par(mfrow = c(2, 2), mar = c(4.5, 4.3, 2.1, 1.1))
for (i in seq_along(mod_list_covs_var)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_covs_var[i])))
  
  plot_covars_single(mod$mod_sum,
                     mod$bugsdata,
                     resp = mod$resp,
                     cov_names = cov_names,
                     covar_std = mod$covar_std,
                     subset = vars_var[i])
  
  mtext(sp_name_covs_var[i], side = 3, adj = 1, line = 0, cex = 1)
  
}
dev.off()

# plot fitted abundance values
pdf(file = "./outputs/plots/FigS1.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.5, 2.1, 1.1))
for (i in seq_along(mod_list_abund)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_abund[i])))
  
  plot_fitted2(mod$mod_sum,
               mod$bugsdata,
               mod$sp_names,
               mod$spp, mod$resp,
               system = length(mod$mod_sum$sysnames),
               'Abundance CPUE')
  
  mtext(sp_name[i], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot fitted biomass values
pdf(file = "./outputs/plots/FigS2.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.5, 2.1, 1.1))
for (i in seq_along(mod_list_biom)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_biom[i])))
  
  plot_fitted2(mod$mod_sum,
               mod$bugsdata,
               mod$sp_names,
               mod$spp, mod$resp,
               system = length(mod$mod_sum$sysnames),
               'Biomass CPUE')
  
  mtext(sp_name[i], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot covariate effects (all species, all predictors)
sp_list_ordered <- c("murraycod", "troutcod", "goldenperch",
                     "silverperch", "murrayriverrainbowfish",
                     "commoncarp")
for (i in seq_along(sp_list_ordered)) {
  
  mod_list_sub <- mod_list_all[grep(sp_list_ordered[i], mod_list_all)][c(2, 6)]
  print(mod_list_sub)
  
  pdf(file = paste0("./outputs/plots/FigS", i + 2, ".pdf"), width = 6, height = 8)
  par(mfrow = c(2, 2), mar = c(4.5, 4.3, 2.1, 1.1))
  for (j in seq_along(mod_list_sub)) {
    
    mod <- get(load(paste0("./outputs/fitted/", mod_list_sub[j])))
    
    for (k in seq_len(mod$bugsdata$Q)) {
      
      plot_covars_single(mod$mod_sum,
                         mod$bugsdata,
                         resp = mod$resp,
                         cov_names = cov_names,
                         covar_std = mod$covar_std,
                         subset = k)
    }
    
  }
  
  dev.off()

}
