# set working directory
setwd("~/Dropbox/research/fish-trends/")

# source helper functions
source("./code/helpers.R")

# load and plot fitted models
mod_list_all <- dir("./outputs/fitted/")

# 
mod_list_all <- mod_list_all[grep("old", mod_list_all, invert = TRUE)]

# subset to the trend and covar_trend models for plotting
mod_list <- mod_list_all[-grep("covar_trend.rds", mod_list_all)]
mod_list <- mod_list[-grep("int_re.rds", mod_list)]

# pull out ppp values for all models
ppp_table <- matrix(NA, nrow = length(mod_list_all), ncol = 2)
for (i in seq_along(mod_list_all)) {
  mod <- get(load(paste0("./outputs/fitted/", mod_list_all[i])))
  ppp_table[i, ] <- round(c(mod$mod_sum$r2, mean(mod$mod_sum$ppps)), 2)
}
rownames(ppp_table) <- mod_list_all

# set up a summary output table and predictor names for labels
cov_names <- c("spawning_variability", "prop_spring_lt", "prop_max_antecedent_lt",
               "prop_summer_lt", "prop_winter_lt", "spawning_temp",
               "number_low_days")
summary_table <- matrix(NA, nrow = length(mod_list), ncol = (2 + 2 * length(cov_names)))
rownames(summary_table) <- seq_len(length(mod_list))
for (i in seq_along(mod_list)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list[i])))
  
  if (!is.na(mod$mod_sum$cov_inc[1])) {
    summary_table[i, ] <- round(c(mod$mod_sum$r2, mean(mod$mod_sum$ppps), mod$mod_sum$cov_inc, mod$mod_sum$cov_or), 2)
  } else {
    summary_table[i, ] <- round(c(mod$mod_sum$r2, mean(mod$mod_sum$ppps), rep(NA, 2 * length(cov_names))), 2)
  }
  rownames(summary_table)[i] <- substr(mod_list[i], 1, nchar(mod_list[i]) - 6)
  
}
colnames(summary_table) <- c("r2", "PPP", rep(cov_names, 2))

sp_list <- unique(sapply(strsplit(mod_list_all, "_"), function(x) x[1]))
hp_out <- vector("list", length = length(sp_list))
r2_out <- vector("list", length = length(sp_list))
for (i in seq_along(sp_list)) {
  
  mod_sub <- mod_list_all[grep(sp_list[i], mod_list_all)]
  mod_load <- vector("list", length = length(mod_sub))
  for (j in seq_along(mod_sub))
    mod_load[[j]] <- get(load(paste0("./outputs/fitted/", mod_sub[j])))
  
  r2_sub <- sapply(mod_load, function(x) x$mod_sum$r2)
  names(r2_sub) <- paste(sapply(mod_load, function(x) x$resp),
                         sapply(mod_load, function(x) x$mod_type),
                         sep = "_")
  
  hp.abund <- hier.part::partition(r2_sub[c("abundance_int_re",
                                            "abundance_trend",
                                            "abundance_covar",
                                            "abundance_covar_trend")],
                                   pcan = 2,
                                   var.names = c("trend", "flow"))
  hp.bioms <- hier.part::partition(r2_sub[c("weight_g_int_re",
                                            "weight_g_trend",
                                            "weight_g_covar",
                                            "weight_g_covar_trend")],
                                   pcan = 2,
                                   var.names = c("trend", "flow"))
  
  r2_out[[i]] <- r2_sub
  hp_out[[i]] <- list(abundance = hp.abund,
                      biomass = hp.bioms)
  
}
names(hp_out) <- names(r2_out) <- sp_list

mod_list_sub <- mod_list_all[grep("covar_trend", mod_list_all)]
conditional_pinc <- matrix(NA, nrow = length(mod_list_sub), ncol = 2 * length(cov_names))
for (i in seq_along(mod_list_sub)) {
  
  mod <- get(load(paste0("./outputs/fitted/", mod_list_sub[i])))

  if (!is.na(mod$mod_sum$cov_inc[1])) {
    conditional_pinc[i, ] <- round(c(mod$mod_sum$cov_inc, mod$mod_sum$cov_or), 2)
  }
  
}
colnames(conditional_pinc) <- rep(cov_names, 2)
rownames(conditional_pinc) <- sapply(strsplit(mod_list_sub, split = "_"),
                                     function(x) paste(x[1], x[2], sep = "_"))
