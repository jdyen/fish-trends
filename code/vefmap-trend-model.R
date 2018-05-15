# R script to fit Bayesian non-linear trend model
# Code by Jim Thomson: jim.thomsson@gmail.com
# Refer to "Abalone Trend Model: User Notes" for instructions

# clear workspace
rm(list = ls())

# load packages
library(R2WinBUGS)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(reshape2)
library(scales)

# set file paths
setwd("~/Dropbox/research/ari-links/")

# source helper functions
source("./code/length_to_mass_calculations.R")
source("./code/length_to_mass_calculations_sra.R")
source("./code/makeBUGSfile_hier.R")

# set bugs directories
bugs_dir <- "/Applications/WinBUGS.app/drive_c/Program Files/WinBUGS14/"
WINE <- "/Applications/Wine.app/Contents/Resources/bin/wine"
WINEPATH <- "/Applications/Wine.app/Contents/Resources/bin/winepath"

# model settings
resp_all <- c("abundance", "biomass")           # response variable
nits <- 10000                                   # number of iterations
debug <- FALSE                                  # debug model
covar_all <- c(TRUE, FALSE)                     # include covariates in model?
plot_outputs <- TRUE                            # plot the fitted trends?

for (banana in seq_along(resp_all)) {
  for (chocolate in seq_along(covar_all)) {

    # load data
    alldat <- read.csv("./data/VEFMAP_FISH_20171024.csv")
    
    # load MC conversion from length to mass
    mc_conv <- read.csv("./data/murray_cod_length_to_mass.csv")
    mc_conv$length_mm <- mc_conv$length_cm * 10
    mc_conv$weight_g <- mc_conv$weight_kg * 1000
    
    # load snags data set
    snags_data <- read.csv("./data/SNAGS_FISH_20171205.csv")
    snags_data$date_new <- format(dmy_hms(snags_data$surveydate), format = "%d/%m/%Y")
    snags_data$YEAR <- sapply(strsplit(snags_data$date_new, "/"),
                              function(x) x[3])
    snags_data$taxonname <- as.character(snags_data$taxonname)
    snags_data$taxonname <- ifelse(snags_data$taxonname == "Yellowbelly",
                                   "Golden perch",
                                   snags_data$taxonname)
    snags_data$taxonname <- factor(snags_data$taxonname)
    snags_data2 <- data.frame(SYSTEM = rep("LOWERMURRAY", nrow(snags_data)),
                              SITE_CODE = paste0("Lm", snags_data$idsite),
                              Reach = rep(1, nrow(snags_data)),
                              geartype = factor(rep("EF/Boat"), nrow(snags_data)),
                              Event_Date = snags_data$date_new,
                              Pass.No = rep(1, nrow(snags_data)),
                              total_no_passes = rep(1, nrow(snags_data)),
                              seconds = snags_data$seconds,
                              Common.Name = snags_data$taxonname,
                              Scientific.Name = snags_data$Scientific.Name,
                              totallength = snags_data$totallength,
                              WEIGHT = snags_data$weight,
                              Total.Sampled = rep(1, nrow(snags_data)),
                              VEFMAP.Stage = rep(NA, nrow(snags_data)),
                              YEAR = as.integer(snags_data$YEAR))
      
    # load ovens data and combine with alldat
    ovens_data <- read.table("./data/vba_ovens_2008_2017.csv", sep = "\t", header = TRUE)
    ovens_data$date_new <- format(dmy(ovens_data$date), format = "%d/%m/%Y")
    ovens_data$YEAR <- sapply(strsplit(ovens_data$date_new, "/"),
                              function(x) x[3])
    ovens_data$species <- as.character(ovens_data$species)
    ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii ",
                                 "Maccullochella peelii peelii",
                                 ovens_data$species)
    ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii",
                                 "Maccullochella peelii peelii",
                                 ovens_data$species)
    ovens_data$species <- factor(ovens_data$species)
    ovens_data$common_name <- alldat$Common.Name[match(ovens_data$species, alldat$Scientific.Name)]
    ovens_data2 <- data.frame(SYSTEM = rep("OVENS", nrow(ovens_data)),
                              SITE_CODE = paste0("Ov", ovens_data$site),
                              Reach = rep(1, nrow(ovens_data)),
                              geartype = ovens_data$gear_type,
                              Event_Date = ovens_data$date_new,
                              Pass.No = rep(1, nrow(ovens_data)),
                              total_no_passes = rep(1, nrow(ovens_data)),
                              seconds = ovens_data$electro_seconds,
                              Common.Name = ovens_data$common_name,
                              Scientific.Name = ovens_data$species,
                              totallength = ovens_data$total_length_mm,
                              WEIGHT = ovens_data$weight_g,
                              Total.Sampled = ovens_data$no_collected,
                              VEFMAP.Stage = rep(NA, nrow(ovens_data)),
                              YEAR = as.integer(ovens_data$YEAR))
    alldat <- rbind(alldat, ovens_data2, snags_data2)

    # load flow data
    source("./code/load-flow-data.R")
    
    resp <- resp_all[banana]
    covar <- covar_all[chocolate]
    
    # set file save name
    save.name <- ifelse(covar, "covar", "trend_only")
    
    # set covariate names
    cov_name <- c("Mean pre-spawning flow", "CV of annual flow", "CV of pre-spawning flow")
    
    # set systems of interest
    system_sub <- c("BROKEN",
                    "LOWERMURRAY",
                    "CAMPASPE",
                    "GLENELG",
                    "GOULBURN",
                    "LODDON",
                    "THOMSON",
                    "OVENS")
    alldat <- alldat[-which(is.na(match(alldat$SYSTEM, system_sub))), ]
    
    # clean up common names
    alldat$Common.Name <- tolower(alldat$Common.Name)
    alldat$Common.Name <- gsub(" ", "", alldat$Common.Name)
    alldat$Common.Name <- gsub("-[[:digit:]]*", "", alldat$Common.Name)
    alldat$Common.Name <- gsub("/", "", alldat$Common.Name)
    alldat$Common.Name <- gsub("sp\\.", "", alldat$Common.Name)
    alldat$Common.Name <- gsub("flatheadedgudgeon", "flatheadgudgeon", alldat$Common.Name)
    alldat$Common.Name <- gsub("europeancarp", "carp", alldat$Common.Name)
    alldat$Common.Name <- gsub("redfinperch", "redfin", alldat$Common.Name)
    alldat$Common.Name <- ifelse(alldat$Common.Name == "weatherloach",
                                 "orientalweatherloach",
                                 alldat$Common.Name)
    alldat$Common.Name <- ifelse(alldat$Common.Name == "rainbowfish",
                                 "murrayriverrainbowfish",
                                 alldat$Common.Name)
    alldat$Common.Name <- ifelse(alldat$Common.Name == "hardyhead",
                                 "unspeckedhardyhead",
                                 alldat$Common.Name)
    alldat$Common.Name <- ifelse(alldat$Common.Name == "blackfish",
                                 "riverblackfish",
                                 alldat$Common.Name)
    
    # remove spp not of interest
    sp_to_rm <- c("carpgoldfishhybrid",
                  "carpgudgeoncomplex",
                  "commonyabby",
                  "easternsnakeneckedturtle",
                  "freshwatershrimp",
                  "glenelgspinyfreshwatercrayfish",
                  "goldfish",
                  "hypseleostris",
                  "murrayspinycrayfish",
                  "nofish",
                  "unidentifiedcod",
                  "unknown",
                  "",
                  "codspp.",
                  "yabbie",
                  "longneckturtle",
                  "gudgeons",
                  "murraycray",
                  "murraycodtroutcodhybrid")
    alldat <- alldat[-which(!is.na(match(alldat$Common.Name, sp_to_rm))), ]
    if (any(is.na(alldat$Common.Name)))
      alldat <- alldat[-which(is.na(alldat$Common.Name)), ]
    
    # create SRA version of weights
    alldat$WEIGHT_SRA <- alldat$WEIGHT
    alldat$IMPUTED <- rep(FALSE, nrow(alldat))
    
    # fill missing weights
    sp_tmp <- unique(alldat$Common.Name)
    for (i in seq_along(sp_tmp)) {
      
      # subset to single species
      dat_tmp <- alldat[which(alldat$Common.Name == sp_tmp[i]), ]

      # check to make sure some weights are missing      
      if (any(is.na(dat_tmp$WEIGHT))) {
        
        # use species-specific equation if it exists, generic otherwise
        if (length_weight_conversion[[sp_tmp[i]]]$n) {
          coefs <- length_weight_conversion[[sp_tmp[i]]]$coef
        } else {
          coefs <- length_weight_conversion$generic
        }
        
        if (length(length_weight_conversion_sra[[sp_tmp[i]]])) {
          coefs_sra <- length_weight_conversion_sra[[sp_tmp[i]]]
        } else {
          coefs_sra <- length_weight_conversion_sra$generic
        }
        
        # subset NA observations
        na_sub <- which(is.na(dat_tmp$WEIGHT))
        
        # estimate weight from length
        dat_tmp$WEIGHT[na_sub] <- exp(coefs[1] + coefs[2] * log(dat_tmp$totallength[na_sub]))
        dat_tmp$WEIGHT_SRA[na_sub] <- 10 ** (coefs_sra["intercept"] +
                                               coefs_sra["slope"] * log((dat_tmp$totallength[na_sub] / coefs_sra["fork_scale_fac"]),
                                                                        base = 10))
        dat_tmp$IMPUTED[na_sub] <- TRUE
        
        # return estimated data to full data set
        alldat[which(alldat$Common.Name == sp_tmp[i]), ] <- dat_tmp
        
      }
    }
    
    # covariate settings
    nbreak <- 10                 # number of spline breakpoints for covariate model
    pc <- c(0.5, 0.3, 0.1, 0.1)  # prior probabilities of no effect, none, one or two breakpoints
    
    # reformat dates
    alldat$Date <- dmy(alldat$Event_Date)
    
    # arrange alldat to have more consistent set of variables
    alldat <- data.frame(date = alldat$Date,
                         year = alldat$YEAR,
                         site = alldat$SITE_CODE,
                         system = alldat$SYSTEM,
                         reach = alldat$Reach,
                         species = alldat$Common.Name,
                         biomass = alldat$WEIGHT,
                         abundance = alldat$Total.Sampled,
                         intensity = (alldat$total_no_passes * alldat$seconds))
    
    alldat$reach_alt <- alldat$reach
    alldat$reach_alt <- ifelse(alldat$system == "BROKEN",
                               ifelse(alldat$reach_alt == 4, 5, alldat$reach_alt), # downstream
                               alldat$reach_alt)
    alldat$reach_alt <- ifelse(alldat$system == "THOMSON",
                               ifelse(alldat$reach_alt == 2, 3,             # downstream
                                      ifelse(alldat$reach_alt == 6, 5,      # upstream
                                             alldat$reach_alt)),
                               alldat$reach_alt)
    alldat$reach_alt <- ifelse(alldat$system == "LODDON",
                               ifelse(alldat$reach_alt == 2, 3,             # downstream
                                      ifelse(alldat$reach_alt == 5, 4,      # upstream
                                             alldat$reach_alt)),
                               alldat$reach_alt)
    
    mean_ann_flow_vec <- mean_spr_flow_vec <- mean_sum_flow_vec <- rep(NA, nrow(alldat))
    maxm_ann_flow_vec <- covn_ann_flow_vec <- rep(NA, nrow(alldat))
    mean_spwn_flow_vec <- cov_spwn_flow_vec <- rep(NA, nrow(alldat))
    for (i in seq_len(nrow(alldat))) {
      row_tmp <- which((mean_flow_stats$system == alldat$system[i]) &
                         (mean_flow_stats$reach == alldat$reach[i]) &
                         (mean_flow_stats$year == alldat$year[i]))
      if (!length(row_tmp))
        row_tmp <- which((mean_flow_stats$system == alldat$system[i]) &
                           (mean_flow_stats$reach == alldat$reach_alt[i]) &
                           (mean_flow_stats$year == alldat$year[i]))
      if (length(row_tmp)) {
        mean_ann_flow_vec[i] <- mean_flow_stats$mean_ann_flow[row_tmp]
        mean_spr_flow_vec[i] <- mean_flow_stats$mean_spr_flow[row_tmp]
        mean_sum_flow_vec[i] <- mean_flow_stats$mean_sum_flow[row_tmp]
        maxm_ann_flow_vec[i] <- mean_flow_stats$max_ann_flow[row_tmp]
        covn_ann_flow_vec[i] <- mean_flow_stats$cov_ann_flow[row_tmp]
        mean_spwn_flow_vec[i] <- mean_flow_stats$mean_spwn_flow[row_tmp]
        cov_spwn_flow_vec[i] <- mean_flow_stats$cov_spwn_flow[row_tmp]
      }
    }
    alldat$mannf <- mean_ann_flow_vec
    alldat$msprf <- mean_spr_flow_vec
    alldat$msumf <- mean_sum_flow_vec
    alldat$maxaf <- maxm_ann_flow_vec
    alldat$covaf <- covn_ann_flow_vec
    alldat$mspwn <- mean_spwn_flow_vec
    alldat$cspwn <- cov_spwn_flow_vec
    
    # set up species names for plots
    sp_names <- data.frame(full = c("Australian bass",
                                    "Australian grayling",
                                    "Australian smelt",
                                    "Black bream",
                                    "Bony bream",
                                    "Brown trout",
                                    "Carp",
                                    "Carp gudgeon",
                                    "Common galaxias",
                                    "Dwarf flathead gudgeon",
                                    "Eastern gambusia",
                                    "Eel",
                                    "Estuary perch",
                                    "Flathead galaxias",
                                    "Flathead gudgeon",
                                    "Golden perch",
                                    "Long-finned eel",
                                    "Luderick",
                                    "Mountain galaxias",
                                    "Murray cod",
                                    "Murray river rainbowfish",
                                    "Obscure galaxias",
                                    "Oriental weather loach",
                                    "Pouched lamprey",
                                    "Pygmy perch",
                                    "Rainbow trout",
                                    "Redfin",
                                    "River blackfish",
                                    "River garfish",
                                    "Roach",
                                    "Sea mullet",
                                    "Short finned eel",
                                    "Short headed lamprey",
                                    "Silver perch",
                                    "Southern pygmy perch",
                                    "Tench",
                                    "Trout cod",
                                    "Tupong",
                                    "Two spined blackfish",
                                    "Unspecked hardyhead",
                                    "Variegated pygmy perch",
                                    "Western carp gudgeon",
                                    "Yarra pygmy perch",
                                    "Yellow eyed mullet"),
                           code = as.character(levels(alldat$species)))
    
    # clean up NAs in reach column
    if (any(is.na(alldat$reach))) {
      alldat <- alldat[-which(is.na(alldat$reach)), ]
    }
    alldat$reach_alt <- NULL
    
    # set up outputs
    r2_all <- NULL
    cov_res_all <- NULL
    
    # set species subset
    species_sub <- c("goldenperch", "murraycod", "murrayriverrainbowfish",
                     "silverperch", "troutcod")
    
    # loop through spatial management units
    for (spp in levels(alldat$species)[match(species_sub, levels(alldat$species))]) {
      
      row_sub <- which(alldat$species == spp)
      
      if (length(row_sub) > 10) {
        
        dat <- alldat[row_sub, ]
        
        dat_melt <- melt(dat, id.vars = c(colnames(dat)[c(1:6, 9)])) 
        row_sub1 <- which(dat_melt$variable == "biomass")
        dat_melt$value[row_sub1] <- dat_melt$value[row_sub1] / dat_melt$intensity[row_sub1]
        row_sub1 <- which(dat_melt$variable == "abundance")
        dat_melt$value[row_sub1] <- dat_melt$value[row_sub1] / dat_melt$intensity[row_sub1]
        
        dat <- dcast(dat_melt,
                     date + year + site + system + reach + species ~ variable,
                     sum)
        dat2 <- dcast(dat_melt,
                      date + year + site + system + reach + species ~ variable,
                      mean)
        dat$mannf <- dat2$mannf
        dat$msprf <- dat2$msprf
        dat$msumf <- dat2$msumf
        dat$maxaf <- dat2$maxaf
        dat$covaf <- dat2$covaf
        dat$mspwn <- dat2$mspwn
        dat$cspwn <- dat2$cspwn
        rm(dat2)
        
        y <- dat[, resp]
        
        if (any(is.na(y))) {
          dat <- dat[-which(is.na(y)), ]
          y <- y[-which(is.na(y))]
        }
        if (any(y == 0)) {
          dat <- dat[-which(y == 0), ]
          y <- y[-which(y == 0)]
        }
        
        if (nrow(dat) > 0) {
          year <- dat$year
          
          yearf <- as.numeric(as.factor(year))
          site <- as.numeric(as.factor(as.integer(dat$site)))
          system <- as.numeric(as.factor(as.integer(dat$system)))
          sysyear <- as.numeric(ordered(as.factor(system):as.factor(yearf)))
          reach <- as.numeric(as.factor(dat$reach))
          Nsystem <- max(system)
          components <- cbind(system, site, yearf, reach)
          
          components <- cbind(components, sysyear)
          Nlevels <- apply(components, 2, max)
          
          Nbatch <- length(Nlevels)
          N <- length(y)
          Nyear <- max(yearf)
          
          kmax <- max(floor(Nyear / 5), 1)
          
          meanyr <- min(dat$year)
          sdyr <- sd(dat$year)
          yr1 <- min(dat$year) - 1
          
          yearx <- c(1:Nyear)
          base.period <- c(min(year):max(year))[c(1, 3)]
          
          sysyr.ind <- matrix(max(components[, Nbatch]) + 1,
                              nrow = Nsystem + 1,
                              ncol = Nyear)
          for (i in 1:Nsystem) {
            for (j in 1:Nyear) {
              sysyr.ind[i, j] <- components[which((system == i) & (yearf == j)), Nbatch][1]
              if (is.na(sysyr.ind[i, j]))
                sysyr.ind[i, j] <- max(components[, Nbatch]) + 1
            }
          } 
          
          trend.period <- cbind(Nyear - c(5:2), rep(Nyear, 4))
          Ntps <- nrow(trend.period)
          
          if (!any(trend.period < 1)) {
            
            if (covar) {
              
              # extract plot variable from predictors
              Xcov <- cbind(dat$mspwn, dat$covaf, dat$cspwn)
              
              # interpolate NAs in covariate data
              Xcov <- apply(Xcov, 2, na_rm_fun)
              Xcov <- sweep(Xcov, 2, apply(Xcov, 2, mean), "-")
              Xcov <- sweep(Xcov, 2, apply(Xcov, 2, sd), "/")
              
              # create dummy matrix to give nonlinear effects
              if (!is.matrix(Xcov)) {
                Xcov <- matrix(Xcov, ncol = 1)
              }
              if (ncol(Xcov) > 1) {
                xdum <- apply(Xcov, 2, create_cov_matrix, nbreak = nbreak)
              } else {
                xdum <- create_cov_matrix(Xcov, nbreak = nbreak)
                xdum <- matrix(xdum, ncol = 1)
              }
              Q <- ncol(Xcov)
              Nplot <- 50
              Xplot <- NULL
              for (i in seq_len(Q)) {
                Xplot <- cbind(Xplot, seq(min(Xcov[, i]), max(Xcov[, i]), length = Nplot))
              }
            }
            
            bugsdata <- list(y = y,
                             yearf = yearf,
                             yearx = yearx,
                             N = N,
                             Nsystem = Nsystem,
                             system = system,
                             comp = components,
                             Nlevels = Nlevels, 
                             Nbatch = Nbatch,
                             Nyear = Nyear,
                             kmax = kmax,
                             b.st = 1,
                             b.end = 3,
                             ryind = sysyr.ind,
                             tpr = trend.period,
                             Ntps = Ntps)
            inits <- function(){list(sd.rand = rep(1.0, Nbatch),
                                     sd = rep(0.5, 3),
                                     alpha = rnorm(1),
                                     mu = rnorm(N),
                                     kmean = rep(1, Nsystem + 1))}
            params <- c("logmean", "sd",
                        "beta", "mu", "trend",
                        "pch", "pinc",
                        "y.sim",
                        "fitted", "fitted.yr",
                        "p.above", "p.ok",
                        "mean.trend", 
                        "pmt", "ppp")
            if (covar) {
              bugsdata <- c(bugsdata,
                            Q = Q,
                            Xcov = list(Xcov),
                            xdum = list(xdum),
                            nbreak = nbreak,
                            pc = list(pc),
                            Nplot = Nplot,
                            Xplot = list(Xplot))
              params <- c(params, "effect", "inc.cov", "beta.cov.all", "cov.plot")
            } 
            file_tmp <- paste0("trend_model_", spp, ".txt")
            filename <- paste0(getwd(), "/code/temp_files/", file_tmp)
            make.model.file.hier(filename, covar = covar, npred = Q, nbreak = nbreak)  
            
            setwd("./code/temp_files")
            fit <- bugs(data = bugsdata, 
                        inits = inits,
                        parameters.to.save = params,
                        model.file = file_tmp,
                        n.chains = 3,
                        n.iter = nits,
                        debug = debug,
                        bugs.directory = bugs_dir,
                        useWINE = TRUE,
                        WINE = WINE,
                        WINEPATH = WINEPATH)
            setwd("../..")
            
            dat$yadj <- y
            dat$qyear <- dmy(paste0("01/01/", dat$year))
            yr <- sort(unique(dat$year))
            yr.date <- dmy(paste0("01/06/", yr))
            
            real_fitted <- fit$summary[grep("y.sim", rownames(fit$summary)), ]
            r2_all <- c(r2_all, cor(real_fitted[, "mean"], y) ** 2)
            names(r2_all)[length(r2_all)] <- spp
            
            fitted <- data.frame(fit$summary[grep("fitted\\[", rownames(fit$summary)), ], yr, yr.date)
            fitted.yr <- data.frame(fit$summary[grep("fitted.yr\\[", rownames(fit$summary)), ], yr, yr.date)
            trends <- fit$summary[grep("trend", rownames(fit$summary)), ]
            trends[, c(1, 3:7)] <- trends[, c(1, 3:7)] - 1
            mean.trends <- data.frame(trends[grep("mean.trend", rownames(trends)), ])
            trends <- trends[-grep("mean.trend", rownames(trends)), ]
            trends <- data.frame(trends,round(fit$summary[grep("pch", rownames(fit$summary)), 1], 3),
                                 round(fit$summary[grep("pinc", rownames(fit$summary)), 1], 3),
                                 round(fit$summary[grep("p.above", rownames(fit$summary)), 1], 3),
                                 round(fit$summary[grep("p.ok", rownames(fit$summary)), 1], 3),
                                 yr, yr.date)
            names(fitted)[3:7] <- names(trends)[3:7] <- names(mean.trends)[3:7] <- paste0("pc", gsub("%", "", colnames(fit$summary)[3:7]))
            names(trends)[10:13] <- c("p.change", "p.increase", "p.above", "p.ok")
            mean.trends$pinc <- as.vector(t(round(fit$mean$pmt, 2)))
            
            prior.ch <- (kmax * 0.5) / Nyear
            trends$ORs <- trends$p.change * (1 - prior.ch) / ((1 - trends$p.change) * prior.ch)
            
            sysnames <- c(unique(as.character(dat$system))[order(unique(components[, "system"]))],
                          "All systems")
            
            ppps <- fit$mean$ppp
            names(ppps) <- sysnames
            
            sysnos <- as.numeric(gsub(",.*", "", gsub("mean.trend\\[", "", rownames(mean.trends))))
            yearnos <- as.numeric(gsub("\\]", "", gsub(".*,", "", gsub("mean.trend\\[", "", rownames(mean.trends)))))
            start.year <- trend.period[yearnos, 1] + max(year) - Nyear
            end.year <- trend.period[yearnos, 2] + max(year) - Nyear
            mean.trends <- data.frame(sysnames[sysnos], start.year, end.year, mean.trends)
            names(mean.trends)[1] <- c("System")
            
            sysnos <- as.numeric(gsub(",.*", "", gsub("trend\\[", "", rownames(trends))))
            yearnos <- as.numeric(gsub("\\]", "", gsub(".*,", "", gsub("trend\\[", "", rownames(trends)))))
            yearvalue <- yearnos + max(year) - Nyear
            trends <- data.frame(sysnames[sysnos], yearvalue, trends)
            names(trends)[1:2] <- c("System", "Year")
            
            if (covar) {
              cov.res <- fit$mean$inc.cov
              if (!is.null(colnames(Xcov)))
                names(cov.res) <- colnames(Xcov)
              cov.prior <- 1 - pc[1]
              OR.cov <- cov.res * (1 - cov.prior) / ((1 - cov.res) * cov.prior)
              cov.res <- data.frame(prob_inc = cov.res, or_inc = OR.cov)
              cov_res_all <- cbind(cov_res_all, cov.res$prob_inc)
              colnames(cov_res_all)[ncol(cov_res_all)] <- spp
            }
            
            # plot results
            if (plot_outputs) {
              
              pdf(paste0("./outputs/", spp, "_", resp, "_", save.name, ".pdf"))
              for (sysnum in seq_len(Nsystem + 1)) {
                
                fitted.r <- fitted[grep(paste0("\\[", sysnum, ","), rownames(fitted)), ]
                trends.r <- trends[grep(paste0("\\[", sysnum, ","), rownames(trends)), ]
                mean.trends.r <- mean.trends[grep(paste0("\\[", sysnum, ","), rownames(mean.trends)), ]
                dat.r <- dat[system == sysnum, ]
                if (length(unique(dat.r$system)) > 1)
                  stop("check data, >1 system code in subset")
                if(sysnum == (Nsystem + 1))
                  dat.r <- dat
                p <- ggplot() + geom_point(data = dat.r, aes(x = date, y = yadj))
                CIs <- ggplot() +
                  geom_ribbon(data = fitted.r,
                              aes(x = yr.date, ymin = pc2.5, ymax = pc97.5),
                              alpha = 0.5,
                              fill = "grey80",
                              color = NA) + 
                  geom_ribbon(data = fitted.r,
                              aes(x = yr.date, ymin = pc25, ymax = pc75),
                              fill = "grey",
                              color = NA) 
                count.plot <- CIs +
                  geom_point(data = dat.r,
                             aes(x = date, y = yadj),
                             color = "grey30") + 
                  geom_line(data = fitted.r,
                            aes(x = yr.date, y = mean)) +
                  ylab(paste0(sp_names$full[match(spp, sp_names$code)], " ", resp)) +
                  xlab("Year") +
                  ggtitle(paste0(sysnames[sysnum])) +
                  theme_bw()
                p2 <- ggplot(data = trends.r) +
                  geom_ribbon(aes(x = yr.date, ymin = pc2.5, ymax = pc97.5),
                              alpha = 0.5,
                              fill = "grey80",
                              color = NA) +
                  geom_ribbon(aes(x = yr.date, ymin = pc25, ymax = pc75),
                              fill = "grey",
                              color = NA) +
                  geom_line(aes(x = yr.date, y = mean), lwd = 2) 
                trend.plot <- p2 +
                  xlab("Year") +
                  ylab("Trend (proportional change per year)") +
                  geom_abline(intercept = 0, slope = 0) +
                  theme_bw()
                
                grid.arrange(count.plot, trend.plot, heights = c(1.5, 1))
                
              }
              
              dev.off()
              
              if (covar) {
                
                pdf(paste0("./outputs/", spp, "_flow_effects_", resp, ".pdf"))
                
                cov_plot_vals <- fit$summary[grep("cov.plot\\[", rownames(fit$summary)), ]
                cov_mean <- matrix(cov_plot_vals[, "mean"], ncol = Q, byrow = TRUE)
                cov_lower <- matrix(cov_plot_vals[, "2.5%"], ncol = Q, byrow = TRUE)
                cov_upper <- matrix(cov_plot_vals[, "97.5%"], ncol = Q, byrow = TRUE)
                
                par(mfrow = c(2, 2))
                for (i in seq_len(Q)) {
                  plot(cov_mean[, i] ~ Xplot[2:nrow(Xplot), i],
                       type = "l", bty = "l", las = 1,
                       xlab = cov_name[i], ylab = paste0("Effect on ", resp),
                       ylim = range(c(cov_mean, cov_lower, cov_upper)))
                  polygon(c(Xplot[2:nrow(Xplot), i], Xplot[nrow(Xplot):2, i]),
                          c(cov_lower[, i], rev(cov_upper[, i])),
                          border = NA,
                          col = "gray50")
                  lines(cov_mean[, i] ~ Xplot[2:nrow(Xplot), i])
                  lines(range(Xplot[2:nrow(Xplot), i]), c(0, 0), lty = 2)
                }
                
                dev.off()
                
              }
              
            }
            
          }
          
        }
        
      }
    }
    if (covar) {
      rownames(cov_res_all) <- c("mean_spwn_flow", "cov_ann_flow", "cov_spwn_flow")
      write.csv(cov_res_all, file = paste0("./outputs/covariate_inclusion_probs_", resp, ".csv"))
    } 
    write.csv(r2_all, file = paste0("./outputs/r2_values_", save.name, "_", resp, ".csv"))
  }
} 
