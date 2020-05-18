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
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")

# load data
# (loaded in code/prepare-data-object.R)
survey_data <- readRDS("data/data-loaded-Nov19.rds")
flow_data <- readRDS("data/flow-data-loaded-Nov19.rds")

# set bugs directories
bugs_dir <- "/Applications/WinBUGS.app/drive_c/Program Files/WinBUGS14/"
WINE <- "/Applications/Wine.app/Contents/Resources/bin/wine"
WINEPATH <- "/Applications/Wine.app/Contents/Resources/bin/winepath"

# model settings
resp_all <- c("weight_g")
# resp_all <- c("abundance", "weight_g")
nits <- 10000
nburn <- 5000
nchain <- 3
debug <- FALSE
mod_type <- c("int_re", "trend", "covar", "covar_trend")

# expand records
expanded_rows <- rep(seq_len(nrow(survey_data)), times = survey_data$abundance)
survey_data <- survey_data[expanded_rows, ]
flow_data <- flow_data[expanded_rows, ]

# fill abundance with ones now that counts are expanded out (one row per fish)
survey_data$abundance <- rep(1, nrow(survey_data))

# kill off some unneeded columns
cols_to_remove <- c("no_collected",
                    "no_observed",
                    "event_date",
                    "gear_type",
                    "scientific_name",
                    "total_length_mm",
                    "temp_c",
                    "do_mgl",
                    "ph",
                    "turbidity_ntu",
                    "secchi_depth_m",
                    "dataset",
                    "common_name",
                    "imputed")
survey_data[, cols_to_remove] <- NULL

# some data sets are missing site IDs
idx <- is.na(survey_data$site)
survey_data$site[idx] <- 
  paste(survey_data$system[idx], survey_data$reach[idx], sep = "_")

# which species?
sp_list <- unique(survey_data$spp_formatted)
# remove Macquarie perch?
sp_list <- sp_list[sp_list != "macquariaaustralasica"]

# which flow variables?
vars_to_include <- c("spawning_variability", "prop_spring_lt", "prop_max_antecedent_lt",
                     "prop_summer_lt", "prop_winter_lt", "spawning_temp",
                     "number_low_days")
flow_tmp <- flow_data[, vars_to_include]

# flow columns are lists for some reason
flow_tmp <- apply(flow_tmp, 2, unlist)
flow_tmp <- as.data.frame(flow_tmp)

for (resp in resp_all) {
  
  for (mod_set in mod_type) {

    # loop through species
    for (spp in sp_list) {

      # prepare data in a BUGS-friendly format
      bugs_set <- prepare_bugs_data(surveys = survey_data,
                                    flow = flow_tmp,
                                    species = spp,
                                    response = resp,
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
                  bugs.directory = bugs_dir,
                  useWINE = TRUE,
                  WINE = WINE,
                  WINEPATH = WINEPATH)

      setwd("../..")
      
      mod_summary <- summarise_fitted(dat = dat, fit = fit,
                                      bugsdata = bugsdata,
                                      mod_type = mod_set)
      
      mod <- list(mod_sum = mod_summary,
                  bugsdata = bugsdata,
                  dat = mod_summary$dat,
                  sp_names = sp_list,
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
      save(mod, file = save_name)
      
    } 
    
  } 
  
}  
