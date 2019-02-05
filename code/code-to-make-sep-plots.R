
# set working directory
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")
source("./code/plot-helpers2.R")

# load and plot fitted models
mod_list_all <- dir("./outputs/fitted/")

# subset to the trend and covar_trend models for plotting
mod_list <- mod_list_all[-grep("covar.rds", mod_list_all)]
mod_list <- mod_list[-grep("int.rds", mod_list)]
mod_list <- mod_list[-grep("int_re.rds", mod_list)]

sp_list <- unique(sapply(strsplit(mod_list_all, "_"), function(x) x[1]))
sp_names <- c('Australian smelt', 'Common carp', 'Golden perch',
              'Murray cod', 'Murray river rainbowfish', 'River blackfish',
              'Silver perch', 'Trout cod')

mod_list <- mod_list[-grep('covar_trend.rds', mod_list)]
for (i in seq_along(mod_list)) {
  
  mod <- get(load(paste0("~/Dropbox/research/fish-trends/outputs/fitted/", mod_list[i])))
  
  
  # plot fitted trends
  for (j in seq_len(mod$bugsdata$Nsystem + 1)) {
    
    spname <- sp_names[sp_list == sapply(strsplit(mod_list[i], '_'), function(x) x[1])]
    type <- ifelse(sapply(strsplit(mod_list[i], '_'), function(x) x[2]) == 'abundance.rds',
                   'abundance', 'biomass')
    
    plot_label <- paste0(spname, '_', type)
    
    sysname <- mod$dat$system[mod$bugsdata$system == j][1]
    if (is.na(sysname)) {
      sysname <- 'ALL'
    }
    
    plot_label <- paste0(plot_label, '_', sysname, '_corrected', '.pdf')

    pdf(file = plot_label, height = 7, width = 7)

    plot_trend(mod$mod_sum,
               mod$bugsdata,
               mod$sp_names,
               mod$spp, mod$resp,
               system = j)


    dev.off()

    pdf(file = paste0("fitted_", plot_label))
    
    plot_fitted2(mod$mod_sum,
                 mod$bugsdata,
                 mod$sp_names,
                 mod$spp, mod$resp,
                 system = j)
    
    dev.off()

  }
  

}
