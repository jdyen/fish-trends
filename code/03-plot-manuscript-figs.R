# set working directory
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")
source("./code/plot-helpers2.R")

# load and plot fitted models
mod_list_all <- dir("./outputs/fitted/")
mod_list_all <- mod_list_all[-grep("old", mod_list_all)]

# subset to the trend and covar_trend models for plotting
mod_list <- mod_list_all[-grep("covar.rds", mod_list_all)]
mod_list <- mod_list[-grep("int_re.rds", mod_list)]

# set up a summary output table and predictor names for labels
cov_names <- c("spawning_variability", "prop_spring_lt", "prop_max_antecedent_lt",
               "prop_summer_lt", "prop_winter_lt", "spawning_temp",
               "number_low_days")
sp_name <- c("maccullochellapeelii" = "Murray cod", 
             "maccullochellamacquariensis" = "Trout cod",
             "macquariaambigua" = "Golden perch",
             "bidyanusbidyanus" = "Silver perch",
             "melanotaeniafluviatilis" = "Murray-Darling rainbowfish",
             "cyprinuscarpio" = "Common carp")

# trends plot
mod_list_abund <- mod_list[grep("abundance.rds", mod_list)][c(4, 3, 5, 1, 6, 2)]
mod_list_biom <- mod_list[grep("weight_g.rds", mod_list)][c(4, 3, 5, 1, 6, 2)]
mod_list_covs <- mod_list_all[grep("covar_trend", mod_list_all)]
mod_list_covs_marginal <- mod_list_all[grep("covar.rds", mod_list_all)]

# plot fitted abundance trends
jpeg(file = "./outputs/plots/Fig1.jpg", width = 6, height = 7, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(4.1, 4.5, 2.1, 1.1))
for (i in seq_along(mod_list_abund)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_abund[i])))
  
  plot_trend(mod$mod_sum,
             mod$bugsdata,
             mod$sp_names,
             mod$spp, mod$resp)
  
  mtext(sp_name[mod$spp], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot fitted biomass trends
jpeg(file = "./outputs/plots/Fig2.jpg", width = 6, height = 7, units = "in", res = 300)
# pdf(file = "./outputs/plots/Fig2_x.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 4.5, 2.1, 1.1))
for (i in seq_along(mod_list_biom)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_biom[i])))
  
  plot_trend(mod$mod_sum,
             mod$bugsdata,
             mod$sp_names,
             mod$spp, mod$resp)
  
  mtext(sp_name[mod$spp], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# fig 3: spring flow effects on abundance of all species (6 panels)
# plot covariate effects (all species, all predictors)
sp_list_ordered <- names(sp_name)

jpeg(file = "outputs/plots/Fig3.jpg", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(3.5, 4, 2.1, 0.5))
k <- 2
for (i in seq_along(sp_list_ordered)) {
  
  mod_list_sub <- mod_list_all[grep(sp_list_ordered[i], mod_list_all)][2]

  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_sub)))
  
  plot_covars_single(mod$mod_sum,
                     mod$bugsdata,
                     resp = mod$resp,
                     cov_names = "",
                     covar_std = mod$covar_std,
                     subset = k)

  if (i %in% c(5, 6))
    mtext("Proportional spring flow", side = 1, line = 2.1, cex = 0.8)  
  if (i %in% c(1, 3, 5))
    mtext("Effect on abundance", side = 2, line = 2.6, cex = 0.8)  
  mtext(sp_name[sp_list_ordered[i]], side = 3, line = 0.4, adj = 0, xpd = TRUE)
  
}
dev.off()

# fig 4: low flow effects on abundance of all species (6 panels)
jpeg(file = "outputs/plots/Fig4.jpg", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(3.5, 4, 2.1, 0.5))
k <- 7
for (i in seq_along(sp_list_ordered)) {
  
  mod_list_sub <- mod_list_all[grep(sp_list_ordered[i], mod_list_all)][2]
  
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_sub)))
  
  plot_covars_single(mod$mod_sum,
                     mod$bugsdata,
                     resp = mod$resp,
                     cov_names = "",
                     covar_std = mod$covar_std,
                     subset = k)
  
  if (i %in% c(5, 6))
    mtext("Number of low-flow days", side = 1, line = 2.1, cex = 0.8)  
  if (i %in% c(1, 3, 5))
    mtext("Effect on abundance", side = 2, line = 2.6, cex = 0.8)  
  mtext(sp_name[sp_list_ordered[i]], side = 3, line = 0.4, adj = 0, xpd = TRUE)
  
}
dev.off()

# plot fitted abundance values
jpeg(file = "./outputs/plots/FigS1.jpg", width = 6, height = 7, units = "in", res = 300)
# pdf(file = "./outputs/plots/FigS1_x.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.5, 2.1, 1.1))
for (i in seq_along(mod_list_abund)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_abund[i])))
  
  plot_fitted2(mod$mod_sum,
               mod$bugsdata,
               mod$sp_names,
               mod$spp, mod$resp,
               system = length(mod$mod_sum$sysnames),
               "Abundance CPUE")
  
  mtext(sp_name[i], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot fitted biomass values
jpeg(file = "./outputs/plots/FigS2.jpg", width = 6, height = 7, units = "in", res = 300)
# pdf(file = "./outputs/plots/FigS2_x.pdf", width = 6, height = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.5, 2.1, 1.1))
for (i in seq_along(mod_list_biom)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_biom[i])))
  
  plot_fitted2(mod$mod_sum,
               mod$bugsdata,
               mod$sp_names,
               mod$spp, mod$resp,
               system = length(mod$mod_sum$sysnames),
               "Biomass CPUE")
  
  mtext(sp_name[i], side = 3, adj = 1, line = 0.5, cex = 1)
  
}
dev.off()

# plot covariate effects (all species, all predictors)
sp_list_ordered <- names(sp_name)
plot_labels <- c("Variability in spawning flow",
                 "Proportional spring flow",
                 "Proportional maximum antecedent flow",
                 "Proportional summer flow",
                 "Proportional winter flow",
                 "Water temperature during spawning period",
                 "Number of low-flow days")
resp_name <- c("weight_g" = "Biomass", "abundance" = "Abundance")
for (i in seq_along(sp_list_ordered)) {
  
  mod_list_sub <- mod_list_all[grep(sp_list_ordered[i], mod_list_all)][c(2, 6)]
  # mod_list_sub <- mod_list_all[grep(sp_list_ordered[i], mod_list_all)][c(1, 5)]

  jpeg(file = paste0("./outputs/plots/FigS", i + 2, ".jpg"), width = 8, height = 8, units = "in", res = 300)
  
  par(mfrow = c(4, 4), mar = c(3.5, 4, 2.1, 0.5))
  
  for (j in seq_along(mod_list_sub)) {
    
    mod <- get(load(paste0("./outputs/fitted/", mod_list_sub[j])))
  
    for (k in seq_len(mod$bugsdata$Q)) {
      
      plot_covars_single(mod$mod_sum,
                         mod$bugsdata,
                         resp = mod$resp,
                         cov_names = plot_labels,
                         covar_std = mod$covar_std,
                         subset = k)
      
      if (k == 4 & j == 1)
        mtext(sp_name[sp_list_ordered[i]], side = 3, line = 0.1, adj = 1, xpd = TRUE)
      if (k == 1)
        mtext(resp_name[mod$resp], side = 3, line = 0.1, adj = 0, xpd = TRUE)

    }
    
    plot(rnorm(100), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    
  }
  
  dev.off()
  
}

# load and plot fitted models by river
mod_list_all <- dir("./outputs/fitted/")
mod_list_all <- mod_list_all[-grep("old", mod_list_all)]

# subset to the trend and covar_trend models for plotting
mod_list <- mod_list_all[-grep("covar.rds", mod_list_all)]
mod_list <- mod_list[-grep("int_re.rds", mod_list)]

sp_list <- unique(sapply(strsplit(mod_list_all, "_"), function(x) x[1]))
sp_names <- c("maccullochellapeelii" = "Murray cod", 
              "maccullochellamacquariensis" = "Trout cod",
              "macquariaambigua" = "Golden perch",
              "bidyanusbidyanus" = "Silver perch",
              "melanotaeniafluviatilis" = "Murray-Darling rainbowfish",
              "cyprinuscarpio" = "Common carp")

mod_list <- mod_list[-grep('covar_trend.rds', mod_list)]
for (i in seq_along(mod_list)) {
  
  mod <- get(load(paste0("outputs/fitted/", mod_list[i])))
  
  
  spname <- sapply(strsplit(mod_list[i], '_'), function(x) x[1])
  type <- ifelse(sapply(strsplit(mod_list[i], '_'), function(x) x[2]) == "abundance.rds",
                 "Abundance", "Biomass")
  plot_label <- paste0(spname, ".jpg")
  
  # setup plot to catch abundance and biomass
  if (i %% 2 == 1) {
    jpeg(file = paste0("outputs/plots/by_system/", plot_label), width = 8, height = 8, units = "in", res = 300, pointsize = 9)
    nplot <- mod$bugsdata$Nsystem + 1
    par(mfrow = c(4, 4), mar = c(4.1, 4.5, 2.8, 1.5))
  }
  
  # plot fitted trends
  for (j in seq_len(nplot)) {
    
    
    sysname <- mod$dat$system[mod$bugsdata$system == j][1]
    if (is.na(sysname)) {
      sysname <- "average"
    }
    
    sys_label <- paste0(toupper(substr(sysname, 1, 1)), substr(sysname, 2, nchar(sysname)))
    
    plot_trend(mod$mod_sum,
               mod$bugsdata,
               mod$sp_names,
               mod$spp, mod$resp,
               system = j)
    
    mtext(sys_label, side = 3, line = 0.1, adj = 1)
    
    if (j == 1)
      mtext(type, side = 3, adj = 0, line = 0.7, cex = 1.4)
    
  }
  
  # fill spaces with blank plots if needed
  if (nplot < 8) {
    for (j in (nplot + 1):8) {
      plot(c(1, 2, 3) ~ 1,
           type = "n", xaxt = "n",
           yaxt = "n", bty = "n",
           xlab = "", ylab = "")
    }
  }
  
  if (i %% 2 == 0)  
    dev.off()
  
}
