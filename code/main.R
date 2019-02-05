# R script to fit Bayesian non-linear trend model
# WinBUGS code by Jim Thomson: jim.thomsson@gmail.com
# Modified by Jian Yen: jdl.yen@gmail.com
# Refer to "Abalone Trend Model: User Notes" for instructions

# clear workspace
rm(list = ls())

# load packages
library(R2WinBUGS)
library(lubridate)
library(reshape2)

# set file paths
setwd("C:/Users/Jian/Dropbox/research/fish-trends/")
# setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/length_to_mass_calculations.R")
source("./code/length_to_mass_calculations_sra.R")
source("./code/helpers.R")

# load data
source("./code/load-survey-data.R")
source("./code/load-flow-data-by-site.R")
source("./code/filter-survey-data.R")

# set bugs directories
# bugs_dir <- "/Applications/WinBUGS.app/drive_c/Program Files/WinBUGS14/"
bugs_dir <- "C:/Users/Jian/Documents/WinBUGS14/"
# WINE <- "/Applications/Wine.app/Contents/Resources/bin/wine"
# WINEPATH <- "/Applications/Wine.app/Contents/Resources/bin/winepath"

# model settings
# resp_all <- c("abundance", "weight")
resp_all <- c("abundance")
nits <- 50000
nburn <- 20000
nchain <- 4
debug <- FALSE
mod_type <- c("int_re", "trend", "covar", "covar_trend")

# expand records
alldat2 <- alldat[!is.na(alldat$recorded), ]
alldat2 <- alldat2[rep(seq_len(nrow(alldat2)), times = alldat2$recorded), ]
alldat2$recorded <- rep(1, nrow(alldat2))
alldat2$abundance <- alldat2$recorded

for (resp in resp_all) {
  
  for (mod_set in mod_type) {

    # loop through species
    # for (spp in levels(alldat2$species)[c(2:5, 7:8)]) {
    for (spp in levels(alldat2$species)[c(8)]) {
        
      bugs_set <- prepare_bugs_data(alldat2,
                                    spp,
                                    resp,
                                    mod_type = mod_set,
                                    nbreak = 10,
                                    pc = c(0.5, 0.3, 0.1, 0.1))
      dat <- bugs_set$dat
      bugsdata <- bugs_set$bugsdata

      make.model.file.hier(bugs_set$filename,
                           mod_type = mod_set,
                           bugsdata = bugsdata,
                           cont = TRUE)  
      
      setwd("./code/temp_files")
      fit <- bugs(data = bugsdata,
                  inits = bugs_set$inits,
                  parameters.to.save = bugs_set$params,
                  model.file = bugs_set$file_tmp,
                  n.chains = nchain,
                  n.iter = nits,
                  n.burnin = nburn,
                  debug = debug,
                  bugs.directory = bugs_dir) # ,
                  # useWINE = TRUE,
                  # WINE = WINE,
                  # WINEPATH = WINEPATH)
      setwd("../..")
      
      mod_summary <- summarise_fitted(dat = dat, fit = fit,
                                      bugsdata = bugsdata,
                                      mod_type = mod_set)
      
      mod_all <- list(mod_sum = mod_summary,
                      bugsdata = bugsdata,
                      dat = mod_summary$dat,
                      sp_names = sp_names,
                      spp = spp,
                      resp = resp,
                      covar_std = bugs_set$covar_std,
                      mod_type = mod_set)
      
      save_name <- paste0("./outputs/fitted/", spp, "_", resp)
      if (mod_set == "int")
        save_name <- paste0(save_name, "_int")
      if (mod_set == "int_re")
        save_name <- paste0(save_name, "_int_re")
      if (mod_set == "covar")
        save_name <- paste0(save_name, "_covar")
      if (mod_set == "covar_trend")
        save_name <- paste0(save_name, "_covar_trend")
      save_name <- paste0(save_name, ".rds")
      save(mod_all, file = save_name)
      
    } 
    
  } 
  
}  
