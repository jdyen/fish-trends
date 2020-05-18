# helpers

# make BUGS file
make.model.file.hier <- function(filename,
                                 mod_type = c("int", "int_re", "covar", "trend", "covar_trend"),
                                 bugsdata,
                                 cont = TRUE) {
  
  sink(filename, append = FALSE)
  
  cat("
      
      Model {
      
      for (i in 1:N) {
      
      y[i] ~ dlnorm(mu[i], tau_main)
      res[i] <- log(y[i]) - mu[i]
      ",
      fill = TRUE)
  
  if(mod_type == "int")
    cat("    mu[i] <- alpha","\n")
  
  if(mod_type == "int_re")
    cat("    mu[i] <- alpha + sum(beta.comp[i, 1:Nbatch])","\n")
  
  if (mod_type == "trend")
    cat("    mu[i] <- alpha + logmean[system[i], yearf[i]] + sum(beta.comp[i, 2:Nbatch])","\n")
  
  if(mod_type == "covar")
    cat("    mu[i] <- alpha + mu.cov[i] + sum(beta.comp[i, 1:Nbatch])","\n")
  
  if(mod_type == "covar_trend")
    cat("    mu[i] <- alpha + logmean[system[i], yearf[i]] + mu.cov[i] + sum(beta.comp[i, 2:Nbatch])","\n")
  
  if (mod_type %in% c("int", "int_re", "covar")) {
    cat("
        for (b in 1:Nbatch) {
        beta.comp[i,b] <- beta[comp[i, b], b]
        }
  }
        ",
        fill = TRUE)
  } else {
    cat("
      for (b in 2:Nbatch) {
      beta.comp[i, b] <- beta[comp[i, b], b]
      }
      }
      
      for (j in 1:Nyear) {
      ",
        fill = TRUE)
    
    if (cont)
      cat("  logmean[Nsystem + 1, j] <- jump.mean[j] - jump.mean[1]", "\n")
    if (!cont)
      cat("  logmean[Nsystem + 1, j] <- jump.mean[j] - jump.mean[1] + jump.mean.step[j] - jump.mean.step[1]", "\n")
    
    
    for (r in seq_len(bugsdata$Nsystem)) {
      cat(paste0("    logmean[", r, ", j] <- beta[", r, ", 1] + logmean[Nsystem + 1, j] + systrend[", r, ", j] - mean(systrend[1:Nsystem, j]) - mean(beta[1:Nsystem, 1])"), "\n")
      if (cont)
        cat(paste0("    systrend[", r, ", j] <- jump.sys", r, "[j] - jump.sys", r, "[1]"), "\n")
      if (!cont)
        cat(paste("    systrend[", r, ", j] <- jump.sys", r, "[j] - jump.sys", r, "[1] + jump.sys.step", r, "[j] - jump.sys.step", r, "[1]"), "\n")
    }
    cat(" 
        }
              
              jump.mean[1:Nyear] <- jump.pw.poly.df.lin(yearx[1:Nyear], kmean[1], tau[2])
              
              ",
        fill = TRUE) 
    
    if (!cont)
      cat("  jump.mean.step[1:Nyear] <- jump.pw.poly.df.gen(yearx[1:Nyear], k.step[1], tau.step[2], 0, 0)", "\n", "\n")
    
    for (r in seq_len(bugsdata$Nsystem)) {
      cat(paste0("  jump.sys", r, "[1:Nyear] <- jump.pw.poly.df.lin(yearx[1:Nyear], kmean[", r + 1, "], tau[3])"), "\n")
      if (!cont)
        cat(paste0("  jump.sys.step", r, "[1:Nyear] <- jump.pw.poly.df.gen(yearx[1:Nyear], k.step[", r+1, "], tau.step[3], 0, 0)"), "\n")
    }
  }
  
  cat("
      
      alpha ~ dnorm(0,0.01)
      tau_main <- 1 / pow(sd_main, 2)
      sd_main ~ dunif(0, 1)
      
      for (b in 1:Nbatch) {
        for (l in 1:(Nlevels[b])) {
          beta[l, b] ~ dnorm(0, tau.rand[b])
        }
        tau.rand[b] <- 1 / pow(sd.rand[b], 2)
        sd.rand[b] ~ dunif(0, 1)
      }
      
      for (t in 1:3) {
      tau[t] <- 1 / pow(sd[t], 2)
      sd[t] ~ dunif(0, 1)
      }
      for (i in 1:(Nsystem + 1)) {
      kmean[i] ~ dbin(0.5, kmax)
      k.step[i] ~ dbin(0.5, kmax)\n",
      fill = TRUE)
  if (mod_type == "covar_trend")
    cat("for (j in 1:Nyear) {\n  fitted[i, j] <- exp(alpha + logmean[i, j] + mu.cov.mean[i, j])\n")
  if (mod_type == "trend")
    cat("for (j in 1:Nyear) {\n  fitted[i, j] <- exp(alpha + logmean[i, j])\n")
  if (mod_type %in% c("int", "int_re", "covar")) {
    if (mod_type == "covar") {
      cat("for (j in 1:Nyear) {\n  fitted[i, j] <- exp(alpha + mu.cov.mean[i, j])\n } \n } \n")
    } else {
      cat("      fitted[i] <- exp(alpha)
        }",
          fill = TRUE)
    }
  } else {
  cat("   pinc[i, j] <- step(trend[i, j] -  1)
          fitted.yr[i, j] <- exp(fitted.yr.log[i, j])\n",
      fill = TRUE)
    if (mod_type == "trend") {
    cat("     fitted.yr.log[i, j] <- alpha + logmean[i, j] + beta[j, 4] + beta[sysyr.ind[i, j], Nbatch]\n",
        fill = TRUE)
    } else {
      cat("     fitted.yr.log[i, j] <- alpha + logmean[i, j] + beta[j, 4] + mu.cov.mean[i, j] + beta[sysyr.ind[i, j], Nbatch]\n",
          fill = TRUE)
    }
    cat("      p.above[i, j] <- step(fitted.yr.log[i, j] - mean(fitted.yr.log[i, b.st:b.end]))
          p.ok[i, j] <- pinc[i,j] * p.above[i,j]
        }
        for (j in 2:Nyear) {
          log(trend[i, j]) <- logmean[i, j] - logmean[i, j - 1]
        }
        for (j in 2:Nyear) {
          pch[i, j] <- 1 - equals(trend[i, j], trend[i, j - 1])
        }
        trend[i, 1] <- trend[i, 2]
        pch[i, 1] ~ dnorm(0, 10000)
      
        for (tp in 1:Ntps) {
          mean.trend[i, tp] <- mean(trend[i, tpr[tp, 1]:tpr[tp, 2]])
          pmt[i, tp] <- step(mean.trend[i, tp] - 1)
        }
      
      }
      
      beta[Nlevels[Nbatch] + 1, Nbatch] <- 0
      ",
      fill = TRUE)
  }
  if (!cont) {
    cat("
      for (t in 1:3) {
        tau.step[t] <- 1 / pow(sd.step[t], 2)
        sd.step[t] ~ dunif(0, 1)
      }
        ",
        fill = TRUE)
  } 
  
  if (mod_type %in% c("covar", "covar_trend")) {
    cat("
      for (i in 1:N) {
        mu.cov[i] <- sum(cov.in[i, ])
        for (j in 1:Q) {
          cov.in[i, j] <- sum(cov.in.seg[, i, j])
          for (k in 1:nbreak) {
            cov.in.seg[k, i, j] <- step(Xcov[i, j] - xdum[(k + 1), j]) *
            beta.cov.all[k, j] * (Xcov[i, j] - xdum[(k + 1), j])
          }
          Xcov.res[i, j] <- Xcov[i, j] - mean(Xcov[, j])
          Xcov.res.sq[i, j] <- pow(Xcov.res[i, j], 2)
          slope.num[i, j] <- (cov.in[i, j] - mean(cov.in[, j])) * Xcov.res[i, j]
        }
      }
      for (i in 1:(Nsystem + 1)) {
        for (j in 1:Nyear) {
          mu.cov.mean[i, j] <- sum(cov.mean[i, j, ])
          for (q in 1:Q) {
            cov.mean[i, j, q] <- sum(cov.mean.seg[, i, j, q])
            for (k in 1:nbreak) {
              cov.mean.seg[k, i, j, q] <- step(flow.mean[i, j, q] - xdum[(k + 1), q]) *
                beta.cov.all[k, q] * (flow.mean[i, j, q] - xdum[(k + 1), q])
            }
          }
        }
      }
        
      for (i in 1:Nplot) {
        for (j in 1:Q) {
          cov.plot[i, j] <- sum(cov.plot.seg[, i, j])
          for (k in 1:nbreak) {
            cov.plot.seg[k, i, j] <- step(Xplot[i, j] - xdum[(k + 1), j]) *
            beta.cov.all[k, j] * (Xplot[i, j] - xdum[(k + 1), j])
          }
        }
      }
        ",
        fill = TRUE)
    
    cat(paste0("  for (i in 1:", bugsdata$nbreak + 1, ") {\n"), fill = TRUE)
    for (i in 1:bugsdata$Q) {
      cat(paste0("    x", i, "[i] <- xdum[i, ", i, "]\n"), fill = TRUE)
    } 
    cat("  }", fill = TRUE)
    
    cat(paste0("  for (i in 1:", bugsdata$nbreak, ") {\n"), fill = TRUE)
    for (i in 1:bugsdata$Q) {
      cat(paste0("    beta.cov.all[i, ", i, "] <- beta.cov", i, "[i + 1] - beta.cov", i, "[i]\n"))
    } 
    cat("  }\n", fill = TRUE)
    
    for (i in 1:bugsdata$Q) {
      cat(paste0("  beta.cov", i, "[1:", bugsdata$nbreak + 1,
                 "] <- jump.pw.poly.df.gen(x", i, "[1:",
                 bugsdata$nbreak + 1, "], kcov[",
                 i, "], tau.cov, 0, 0)\n"))
    } 
    
    cat("
        for (i in 1:Q) {
        kcov[i] ~ dshifted.cat(pc[], 0)
        effect[i] <- sum(slope.num[, i]) / sum(Xcov.res.sq[, i])
        inc.cov[i] <- 1 - equals(kcov[i], 0)
        }
        tau.cov <- 1 / pow(sdcov, 2)
        sdcov ~ dunif(0, 1)",
        fill = TRUE)
  } 
  cat("
      # ppp calcs
      for (i in 1:N) {
      y.sim[i] ~ dlnorm(mu[i], tau_main)
      res.sim[i] <- log(y.sim[i]) - mu[i]
      chi[i, 1] <- res[i] * res[i] / mu[i]
      chi[i, 2] <- res.sim[i] * res.sim[i] / mu[i]
      for (j in 1:Nsystem) {
      chi.r[i, 1, j] <- chi[i, 1] * equals(system[i], j)
      chi.r[i, 2, j] <- chi[i, 2] * equals(system[i], j)
      }
      }
      for (t in 1:2) {
      chi.sum[t] <- sum(chi[, t])
      for (j in 1:Nsystem) {
      chi.sum.r[j, t] <- sum(chi.r[, t, j])
      }
      }
      ppp[Nsystem + 1] <- step(chi.sum[2] - chi.sum[1])
      for (j in 1:Nsystem) {
      ppp[j] <- step(chi.sum.r[j, 2] - chi.sum.r[j, 1])
      }
      
      ",
      fill = TRUE)
  
  cat("}", "\n")
  
  sink()
} 

create_cov_matrix <- function(x, nbreak = 10) {
  
  xdum <- seq(min(x, na.rm = TRUE),
              max(x, na.rm = TRUE),
              length = (nbreak + 1))
  xdum <- c(min(x, na.rm = TRUE) - mean(diff(xdum)), xdum[1:nbreak])
  
  xdum
  
} 

na_rm_fun <- function(x) {
  
  if (any(is.na(x)))
    x[which(is.na(x))] <- median(x, na.rm = TRUE)
  
  x
  
} 

rebase_factor <- function(x) {
  as.integer(as.factor(x))
}

# prepare BUGS inputs (dummy matrix)
prepare_bugs_data <- function(surveys, flow, species, response, mod_type, nbreak, pc) {

  dat <- cbind(surveys, flow)  
  dat <- dat[dat$spp_formatted == species, ]
  dat$year <- year(dat$date_formatted)

  covar_names <- colnames(flow)
  
  dat_melt <- melt(dat,
                   id.vars = c("date_formatted", "year", "site", "system", "reach",
                               "spp_formatted", "ef_seconds_total")) 
  idx <- dat_melt$variable == "weight_g"
  idy <- dat_melt$variable == "abundance"
  dat_melt$value[idx] <- dat_melt$value[idx] / dat_melt$ef_seconds_total[idx]
  dat_melt$value[idy] <- dat_melt$value[idy] / dat_melt$ef_seconds_total[idy]
  
  dat <- dcast(dat_melt,
               date_formatted + year + site + system + reach + spp_formatted ~ variable,
               sum)
  dat2 <- dcast(dat_melt,
                date_formatted + year + site + system + reach + spp_formatted ~ variable,
                mean)
  dat[, covar_names] <- dat2[, covar_names]
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
  
  year <- dat$year
  
  yearf <- rebase_factor(year)
  site <- rebase_factor(dat$site)
  system <- rebase_factor(dat$system)
  site <- as.integer(ordered(as.factor(system):as.factor(site)))
  sysyear <- as.integer(ordered(as.factor(system):as.factor(yearf)))
  reach <- as.integer(dat$reach)
  reach <- as.integer(ordered(as.factor(system):as.factor(reach)))
  Nsystem <- max(system)
  
  components <- cbind(system, site, reach)
  if (mod_type %in% c("trend", "covar_trend"))
    components <- cbind(components, yearf, sysyear)

  Nlevels <- apply(components, 2, max)
  
  Nbatch <- length(Nlevels)
  N <- length(y)
  Nyear <- max(yearf)
  
  kmax <- max(floor(Nyear / 4), 1)
  
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
  
  if (mod_type %in% c("covar", "covar_trend")) {
    
    # extract plot variable from predictors
    Xcov <- dat[, covar_names]
    
    # interpolate NAs in covariate data
    Xcov <- apply(Xcov, 2, na_rm_fun)
    x_means <- apply(Xcov, 2, mean)
    x_sds <- apply(Xcov, 2, sd)
    Xcov <- sweep(Xcov, 2, x_means, "-")
    Xcov <- sweep(Xcov, 2, x_sds, "/")
    
    # calculate annual, system-specific mean values on standardized scale
    flow.mean <- array(NA, dim = c(Nsystem + 1, Nyear, ncol(Xcov)))
    for (q in seq_len(ncol(Xcov))) {
      flow.mean[seq_len(Nsystem), , q] <- tapply(Xcov[, q], list(system, yearf), mean)
      flow.mean[(Nsystem + 1), , q] <- tapply(Xcov[, q], list(yearf), mean)
      flow.mean[, , q] <- sweep(flow.mean[, , q], 1, apply(flow.mean[, , q], 1, mean, na.rm = TRUE),
                                function(x, y) ifelse(is.na(x), y, x))
    }

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
  
  inits <- function(){list(sd.rand = rep(1.0, Nbatch),
                           sd = rep(0.5, 3),
                           alpha = rnorm(1),
                           kmean = rep(1, Nsystem + 1))}
  params <- c("sd",
              "beta",
              "mu",
              "y.sim",
              "fitted",
              "ppp")
  if (mod_type %in% c("int", "int_re")) {
    bugsdata <- list(y = y,
                     N = N,
                     Nsystem = Nsystem,
                     system = system,
                     comp = components,
                     Nlevels = Nlevels, 
                     Nbatch = Nbatch,
                     kmax = kmax)
  } else {
    if (mod_type %in% c("trend", "covar_trend")) {
      bugsdata <- list(y = y,
                       N = N,
                       Nsystem = Nsystem,
                       system = system,
                       comp = components,
                       Nlevels = Nlevels, 
                       Nbatch = Nbatch,
                       kmax = kmax,
                       yearf = yearf,
                       yearx = yearx,
                       Nyear = Nyear,
                       sysyr.ind = sysyr.ind,
                       b.st = 1,
                       b.end = 3,
                       tpr = trend.period,
                       Ntps = Ntps)
    } else {
      bugsdata <- list(y = y,
                       N = N,
                       Nsystem = Nsystem,
                       system = system,
                       comp = components,
                       Nlevels = Nlevels, 
                       Nbatch = Nbatch,
                       kmax = kmax,
                       Nyear = Nyear)
    }
  } 
  if (mod_type %in% c("trend", "covar_trend")) {
    params <- c(params, "logmean", "trend", "pch", "pinc",
                "fitted.yr", "p.above", "p.ok", "mean.trend",
                "pmt")
  }
  covar_std <- NULL
  
  if (mod_type %in% c("covar", "covar_trend")) {
    bugsdata <- c(bugsdata,
                  Q = Q,
                  Xcov = list(Xcov),
                  xdum = list(xdum),
                  nbreak = nbreak,
                  pc = list(pc),
                  Nplot = Nplot,
                  Xplot = list(Xplot),
                  flow.mean = list(flow.mean))
    params <- c(params, "effect", "inc.cov", "beta.cov.all", "cov.plot")
    covar_std <- list(mean = x_means,
                      x = x_sds)
  }   
  
  file_tmp <- paste0("trend_model_", spp, ".txt")
  filename <- paste0(getwd(), "/code/temp_files/", file_tmp)
  
  # return outputs
  list(
    dat = dat,
    bugsdata = bugsdata,
    inits = inits,
    params = params,
    covar_std = covar_std,
    filename = filename,
    file_tmp = file_tmp
  )

}   

# extract model summary stats
summarise_fitted <- function(dat, fit, bugsdata, mod_type) {
  
  dat$yadj <- bugsdata$y
  dat$qyear <- dmy(paste0("01/01/", dat$year))
  yr <- sort(unique(dat$year))
  yr.date <- dmy(paste0("01/06/", yr))
  
  real_fitted <- (fit$summary[grep("mu", rownames(fit$summary)), ])
  r2_tmp <- cor(real_fitted[, "mean"], log(bugsdata$y)) ** 2
   
  if (mod_type %in% c("int", "int_re", "covar"))
    fitted <- data.frame(fit$summary[grep("fitted\\[", rownames(fit$summary)), ])
  else
    fitted <- data.frame(fit$summary[grep("fitted\\[", rownames(fit$summary)), ], yr, yr.date)
  if (mod_type %in% c("trend", "covar_trend")) {
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
    
    prior.ch <- (bugsdata$kmax * 0.5) / bugsdata$Nyear
    trends$ORs <- trends$p.change * (1 - prior.ch) / ((1 - trends$p.change) * prior.ch)
  }
  
  sysnames <- c(unique(as.character(dat$system))[order(unique(bugsdata$comp[, "system"]))],
                "All systems")
  
  ppps <- fit$mean$ppp
  names(ppps) <- sysnames
  
  if (mod_type == "trend" | mod_type == "covar_trend") {
    sysnos <- as.numeric(gsub(",.*", "", gsub("mean.trend\\[", "", rownames(mean.trends))))
    yearnos <- as.numeric(gsub("\\]", "", gsub(".*,", "", gsub("mean.trend\\[", "", rownames(mean.trends)))))
    start.year <- bugsdata$tpr[yearnos, 1] + max(dat$year) - bugsdata$Nyear
    end.year <- bugsdata$tpr[yearnos, 2] + max(dat$year) - bugsdata$Nyear
    mean.trends <- data.frame(sysnames[sysnos], start.year, end.year, mean.trends)
    names(mean.trends)[1] <- c("System")
    
    sysnos <- as.numeric(gsub(",.*", "", gsub("trend\\[", "", rownames(trends))))
    yearnos <- as.numeric(gsub("\\]", "", gsub(".*,", "", gsub("trend\\[", "", rownames(trends)))))
    yearvalue <- yearnos + max(dat$year) - bugsdata$Nyear
    trends <- data.frame(sysnames[sysnos], yearvalue, trends)
    names(trends)[1:2] <- c("System", "Year")
  } else {
    fitted.yr <- NULL
    trends <- NULL
    mean.trends <- NULL
  }
  
  if (mod_type != "trend") {
    cov_res <- fit$mean$inc.cov
    if (!is.null(colnames(bugsdata$Xcov)))
      names(cov_res) <- colnames(bugsdata$Xcov)
    cov_prior <- 1 - bugsdata$pc[1]
    or_cov <- cov_res * (1 - cov_prior) / ((1 - cov_res) * cov_prior)
    cov_inc <- cov_res
    cov_or <- or_cov
    cov_plot_vals <- fit$summary[grep("cov.plot\\[", rownames(fit$summary)), ]
  } else {  
    cov_inc <- rep(NA, 2)
    cov_or <- rep(NA, 2)
    cov_plot_vals <- NULL
  }  
  
  out <- list(r2 = r2_tmp,
              ppps = ppps,
              cov_inc = cov_inc,
              cov_or = cov_or,
              real_fitted = real_fitted,
              fitted = fitted,
              fitted.yr = fitted.yr,
              trends = trends,
              mean.trends = mean.trends,
              date = date,
              yr.date = yr.date,
              sysnames = sysnames,
              cov_plot_vals = cov_plot_vals,
              dat = dat)
  
  out
  
}  

# plot fitted models
plot_fitted <- function(mod_sum, bugsdata, sp_names, spp, resp, ylab = NULL) {
  
  sysnum <- bugsdata$Nsystem + 1
  fitted.r <- mod_sum$fitted[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$fitted)), ]
  
  xvals <- fitted.r$yr
  upper95 <- fitted.r$pc97.5
  lower95 <- fitted.r$pc2.5
  upper50 <- fitted.r$pc75
  lower50 <- fitted.r$pc25
  plot(mod_sum$dat$yadj ~ mod_sum$dat$year,
       ylim = range(c(mod_sum$dat$yadj, upper95, lower95, upper50, lower50)),
       type = "n",
       las = 1, bty = "l",
       xlab = "", ylab = "")
  polygon(c(xvals, rev(xvals)),
          c(upper95, rev(lower95)),
          col = ggplot2::alpha("grey75", 0.5),
          border = NA)
  polygon(c(xvals, rev(xvals)),
          c(upper50, rev(lower50)),
          col = ggplot2::alpha("grey50", 0.5),
          border = NA)
  lines(fitted.r$pc50 ~ xvals, lwd = 2,
        col = ggplot2::alpha("grey30", 0.5))
  points(mod_sum$dat$yadj ~ mod_sum$dat$year,
         pch = 16,
         col = ggplot2::alpha("grey30", 0.5))
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  if (is.null(ylab))
    ylab <- paste0(resp, " CPUE")
  mtext(ylab, side = 2, adj = 0.5, line = 3.0)
  
} 

plot_trend <- function(mod_sum, bugsdata, sp_names, spp, resp, system = NULL) {
  
  if (is.null(system)) {
    sysnum <- bugsdata$Nsystem + 1
  } else {
    sysnum <- system
  }
  
  trends.r <- mod_sum$trends[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$trends)), ]

  trend_upper95 <- trends.r$pc97.5
  trend_lower95 <- trends.r$pc2.5
  trend_upper50 <- trends.r$pc75
  trend_lower50 <- trends.r$pc25
  x_trend <- trends.r$Year
  ylim_set <- range(c(trend_lower95, trend_lower50, trend_upper50, trend_upper95))
  plot(trend_upper50 ~ x_trend,
       type = "n",
       ylim = ylim_set,
       bty = "l", las = 1,
       xlab = "", ylab = "")
  polygon(c(x_trend, rev(x_trend)),
          c(trend_upper95, rev(trend_lower95)),
          col = ggplot2::alpha("grey80", 0.5),
          border = NA)
  polygon(c(x_trend, rev(x_trend)),
          c(trend_upper50, rev(trend_lower50)),
          col = ggplot2::alpha("grey60", 0.5),
          border = NA)
  lines(trends.r$pc50 ~ x_trend,
        lwd = 2)
  lines(c(min(x_trend) - 1, max(x_trend) + 1), c(0, 0), lty = 2)
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  mtext("Trend", side = 2, adj = 0.5, line = 3.0)
  
}

# plot covariate effects
plot_covars <- function(mod_sum, bugsdata, resp, cov_names = NULL, covar_std) {
  
  resp <- ifelse(resp == "weight", "biomass", resp)
  
  if (is.null(cov_names))
    cov_names <- letters
  
  cov_mean <- exp(matrix(mod_sum$cov_plot_vals[, "mean"], ncol = bugsdata$Q, byrow = TRUE))
  cov_lower <- exp(matrix(mod_sum$cov_plot_vals[, "2.5%"], ncol = bugsdata$Q, byrow = TRUE))
  cov_upper <- exp(matrix(mod_sum$cov_plot_vals[, "97.5%"], ncol = bugsdata$Q, byrow = TRUE))
  
  if (is.null(cov_names))
    cov_names <- letters[seq_len(bugsdata$Q)]
  
  oldmfrow <- par()$mfrow
  round_val <- list(c(-2, rep(-3, 4)),
                    2, 2, 2, 2, 2, 2)
  par(mfrow = c(1, 2))
  for (i in seq_len(bugsdata$Q)) {
    x_set <- bugsdata$Xplot[2:nrow(bugsdata$Xplot), i]
    x_adj <- covar_std$mean[i] + (x_set * covar_std$x[i])
    plot(cov_mean[, i] ~ x_set,
         type = "l", bty = "l", las = 1,
         xaxt = "n",
         xlab = cov_names[i], ylab = paste0("Effect on ", resp),
         ylim = range(c(0, cov_mean[, i], cov_lower[, i], cov_upper[, i])))
    axis(1, at = seq(min(x_set), max(x_set), length = 5),
         labels = round(seq(min(x_adj), max(x_adj), length = 5), round_val[[i]]))
    polygon(c(bugsdata$Xplot[2:nrow(bugsdata$Xplot), i], bugsdata$Xplot[nrow(bugsdata$Xplot):2, i]),
            c(cov_lower[, i], rev(cov_upper[, i])),
            border = NA,
            col = ggplot2::alpha("gray50", 0.5))
    lines(cov_mean[, i] ~ bugsdata$Xplot[2:nrow(bugsdata$Xplot), i])
    lines(range(bugsdata$Xplot[2:nrow(bugsdata$Xplot), i]), c(1, 1), lty = 2)
  }  
  
}

# plot covariate effects
plot_covars_single <- function(mod_sum, bugsdata, resp, cov_names = NULL, covar_std,
                               subset = 1) {
  
  if (is.null(cov_names))
    cov_names <- letters
  
  cov_mean <- exp(matrix(mod_sum$cov_plot_vals[, "mean"], ncol = bugsdata$Q, byrow = TRUE))
  cov_lower <- exp(matrix(mod_sum$cov_plot_vals[, "2.5%"], ncol = bugsdata$Q, byrow = TRUE))
  cov_upper <- exp(matrix(mod_sum$cov_plot_vals[, "97.5%"], ncol = bugsdata$Q, byrow = TRUE))
  
  if (is.null(cov_names))
    cov_names <- letters[seq_len(bugsdata$Q)]
  
  oldmfrow <- par()$mfrow
  round_val <- list(c(-2, rep(-3, 4)),
                    2, 2, 2, 2, 2, 2)
  i <- subset
  x_set <- bugsdata$Xplot[2:nrow(bugsdata$Xplot), i]
  x_adj <- covar_std$mean[i] + (x_set * covar_std$x[i])
  plot(cov_mean[, i] ~ x_set,
       type = "l", bty = "l", las = 1,
       xaxt = "n",
       xlab = "", ylab = "",
       ylim = range(c(0, cov_mean[, i], cov_lower[, i], cov_upper[, i])))
  axis(1, at = seq(min(x_set), max(x_set), length = 5),
       labels = round(seq(min(x_adj), max(x_adj), length = 5), round_val[[i]]))
  polygon(c(bugsdata$Xplot[2:nrow(bugsdata$Xplot), i], bugsdata$Xplot[nrow(bugsdata$Xplot):2, i]),
          c(cov_lower[, i], rev(cov_upper[, i])),
          border = NA,
          col = ggplot2::alpha("gray50", 0.5))
  lines(cov_mean[, i] ~ bugsdata$Xplot[2:nrow(bugsdata$Xplot), i])
  lines(range(bugsdata$Xplot[2:nrow(bugsdata$Xplot), i]), c(1, 1), lty = 2)
  
  
  mtext(cov_names[i], side = 1, line = 2.1, adj = 0.5, cex = 0.7)
  
  if (k == 1 | k == 5) {
    mtext(paste0("Effect on ", ifelse(resp == "weight_g", "biomass", "abundance")),
        side = 2, line = 2.6, adj = 0.5, cex = 0.7)
  }

}

plot_fitted2 <- function(mod_sum, bugsdata, sp_names, spp, resp, system = 1, ylab = NULL) {
  
  sysnum <- system
  fitted.r <- mod_sum$fitted[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$fitted)), ]
  
  sysname <- mod$mod_sum$sysnames[sysnum]
  subset <- mod_sum$dat$system == sysname
  if (sysname == "All systems")
    subset <- seq_len(nrow(mod_sum$dat))
  
  xvals <- fitted.r$yr
  upper95 <- fitted.r$pc97.5
  lower95 <- fitted.r$pc2.5
  upper50 <- fitted.r$pc75
  lower50 <- fitted.r$pc25
  plot(mod_sum$dat$yadj[subset] ~ mod_sum$dat$year[subset],
       ylim = range(c(mod_sum$dat$yadj, upper95, lower95, upper50, lower50)),
       type = "n",
       log = "y",
       las = 1, bty = "l",
       xlab = "", ylab = "")
  polygon(c(xvals, rev(xvals)),
          c(upper95, rev(lower95)),
          col = ggplot2::alpha("grey75", 0.5),
          border = NA)
  polygon(c(xvals, rev(xvals)),
          c(upper50, rev(lower50)),
          col = ggplot2::alpha("grey50", 0.5),
          border = NA)
  lines(fitted.r$pc50 ~ xvals, lwd = 2,
        col = ggplot2::alpha("grey30", 0.5))
  points(mod_sum$dat$yadj[subset] ~ mod_sum$dat$year[subset],
         pch = 16,
         col = ggplot2::alpha("grey30", 0.5))
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  if (is.null(ylab))
    ylab <- paste0(ifelse(resp == "abundance", "Abundance ", "Biomass "), ' CPUE')
  mtext(ylab, side = 2, adj = 0.5, line = 4)
  
}

system_switch_fun <- function(x) {
  
  if (!is.character(x)) {
    x <- as.character(x)
  }
  
  out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    out[i] <- switch(substr(x[i], 1, 2),
                     "GO" = "GOULBURN",
                     "BR" = "BROKEN",
                     "PC" = "PYRAMID CK",
                     "LO" = "LODDON",
                     "CA" = "CAMPASPE")
  }  
  
  out
  
}
