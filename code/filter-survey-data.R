# filter data down to a species and river subset

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
