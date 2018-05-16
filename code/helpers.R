# helpers

# make BUGS file
make.model.file.hier <- function(filename,
                                 covar = FALSE,
                                 bugsdata,
                                 cont = 1) {
  
  sink(filename, append = FALSE)
  
  cat("
      
      Model {
      
      for (i in 1:N) {
      
      y[i] ~ dlnorm(mu[i], tau_main)
      res[i] <- y[i] - mu[i]
      mu[i] ~ dnorm(mu0[i], tau[1])
      proc.res[i] <- mu[i] - mu00[i]
      ",
      fill = TRUE)
  
  if (!covar)
    cat("    mu00[i] <- alpha + logmean[system[i], yearf[i]] + sum(beta.comp[i, 2:Nbatch])","\n")
  
  if(covar)
    cat("    mu00[i] <- alpha + logmean[system[i], yearf[i]] + mu.cov[i] + sum(beta.comp[i, 2:Nbatch])","\n")
  
  cat("
      mu0[i] <- mu00[i]
      for (b in 2:Nbatch) {
      beta.comp[i,b] <- beta[comp[i, b], b]
      }
      }
      
      for (j in 1:Nyear) {
      ",
      fill = TRUE)
  
  if (cont)
    cat("  logmean[Nsystem + 1, j] <- jump.mean[j] - jump.mean[1]", "\n")
  if (!cont)
    cat("  logmean[Nsystem + 1, j] <- jump.mean[j] - jump.mean[1] + jump.mean.step[j] - jump.mean.step[1]", "\n")
  
  
  for (r in 1:bugsdata$Nsystem) {
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
  
  for (r in 1:bugsdata$Nsystem) {
    cat(paste0("  jump.sys", r, "[1:Nyear] <- jump.pw.poly.df.lin(yearx[1:Nyear], kmean[", r + 1, "], tau[3])"), "\n")
    if (!cont)
      cat(paste0("  jump.sys.step", r, "[1:Nyear] <- jump.pw.poly.df.gen(yearx[1:Nyear], k.step[", r+1, "], tau.step[3], 0, 0)"), "\n")
  }
  
  cat("
      
      rho ~ dnorm(0,2)
      alpha ~ dnorm(0,0.001)
      tau_main <- 1 / pow(sd_main, 2)
      sd_main ~ dunif(0, 10)
      
      for (b in 1:Nbatch) {
      for (l in 1:(Nlevels[b])) {
      beta[l, b] ~ dnorm(0, tau.rand[b])
      }
      tau.rand[b] <- 1 / pow(sd.rand[b], 2)
      sd.rand[b] ~ dunif(0, 5)
      }
      
      for (t in 1:3) {
      tau[t] <- 1 / pow(sd[t], 2)
      sd[t] ~ dunif(0, 1)
      }
      for (i in 1:(Nsystem + 1)) {
      kmean[i] ~ dbin(0.5, kmax)
      k.step[i] ~ dbin(0.5, kmax)
      for (j in 1:Nyear) {\n",
      fill = TRUE)
  if (covar) {
    cat("      fitted[i, j] <- exp(alpha + logmean[i, j] + mu.cov[i])\n")
  } else {
    cat("      fitted[i, j] <- exp(alpha + logmean[i, j])\n")
  }
  cat("      pinc[i, j] <- step(trend[i, j] - 1)
      fitted.yr[i, j] <- exp(fitted.yr.log[i, j])
      fitted.yr.log[i, j] <- alpha + logmean[i, j] + beta[j, 3] + beta[ryind[i, j], Nbatch]
      p.above[i, j] <- step(fitted.yr.log[i, j] - mean(fitted.yr.log[i, b.st:b.end]))
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
  if (!cont) {
    cat("
        for (t in 1:3) {
        tau.step[t] <- 1 / pow(sd.step[t], 2)
        sd.step[t] ~ dunif(0, 1)
        }
        ",
        fill = TRUE)
  }
  
  if (covar) {
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
                 i, "], tau.cov, 0, 0)"))
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
      res.sim[i] <- y.sim[i] - mu[i]
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
  xdum <- c(min(x, na.rm = TRUE) - 1, xdum[1:nbreak])
  
  xdum
  
}

na_rm_fun <- function(x) {
  
  if (any(is.na(x)))
    x[which(is.na(x))] <- median(x, na.rm = TRUE)
  
  x
  
}

# prepare BUGS inputs (dummy matrix)
prepare_bugs_data <- function(spp, resp, covar, nbreak, pc) {
  
  dat <- alldat[alldat$species == spp, ]
  
  dat_melt <- melt(dat, id.vars = c(colnames(dat)[c(1:6, 9)])) 
  dat_melt$value[dat_melt$variable == "biomass"] <- dat_melt$value[dat_melt$variable == "biomass"] / dat_melt$intensity[dat_melt$variable == "biomass"]
  dat_melt$value[dat_melt$variable == "abundance"] <- dat_melt$value[dat_melt$variable == "abundance"] / dat_melt$intensity[dat_melt$variable == "abundance"]
  
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
  bugsdata = list(y = y,
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
  
  out <- list(dat = dat,
              bugsdata = bugsdata,
              inits = inits,
              params = params,
              filename = filename,
              file_tmp = file_tmp)
  
  out
  
}  

# extract model summary stats
summarise_fitted <- function(dat, fit, bugsdata) {
  
  dat$yadj <- bugsdata$y
  dat$qyear <- dmy(paste0("01/01/", dat$year))
  yr <- sort(unique(dat$year))
  yr.date <- dmy(paste0("01/06/", yr))
  
  real_fitted <- fit$summary[grep("y.sim", rownames(fit$summary)), ]
  r2_tmp <- cor(real_fitted[, "mean"], bugsdata$y) ** 2
   
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
  
  prior.ch <- (bugsdata$kmax * 0.5) / bugsdata$Nyear
  trends$ORs <- trends$p.change * (1 - prior.ch) / ((1 - trends$p.change) * prior.ch)
  
  sysnames <- c(unique(as.character(dat$system))[order(unique(bugsdata$comp[, "system"]))],
                "All systems")
  
  ppps <- fit$mean$ppp
  names(ppps) <- sysnames
  
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
  
  if (covar) {
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
              fitted = fitted,
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
plot_fitted <- function(mod_sum, bugsdata, sp_names, spp, resp) {
  
  for (sysnum in seq_len(bugsdata$Nsystem + 1)) {
    
    fitted.r <- mod_sum$fitted[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$fitted)), ]
    trends.r <- mod_sum$trends[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$trends)), ]
    mean.trends.r <- mod_sum$mean.trends[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$mean.trends)), ]
    dat.r <- mod_sum$dat[mod_sum$dat$system == sysnum, ]
    if (length(unique(dat.r$system)) > 1)
      stop("check data, >1 system code in subset")
    if(sysnum == (bugsdata$Nsystem + 1))
      dat.r <- mod_sum$dat
    p <- ggplot() + geom_point(data = dat.r, aes(x = date, y = yadj))
    CIs <- ggplot() +
      geom_ribbon(data = fitted.r,
                  aes(x = mod_sum$yr.date, ymin = pc2.5, ymax = pc97.5),
                  alpha = 0.5,
                  fill = "grey80",
                  color = NA) + 
      geom_ribbon(data = fitted.r,
                  aes(x = mod_sum$yr.date, ymin = pc25, ymax = pc75),
                  fill = "grey",
                  color = NA) 
    count.plot <- CIs +
      geom_point(data = dat.r,
                 aes(x = date, y = yadj),
                 color = "grey30") + 
      geom_line(data = fitted.r,
                aes(x = mod_sum$yr.date, y = mean)) +
      ylab(paste0(sp_names$full[match(spp, sp_names$code)], " ", resp)) +
      xlab("Year") +
      ggtitle(paste0(mod_sum$sysnames[sysnum])) +
      theme_bw()
    p2 <- ggplot(data = trends.r) +
      geom_ribbon(aes(x = mod_sum$yr.date, ymin = pc2.5, ymax = pc97.5),
                  alpha = 0.5,
                  fill = "grey80",
                  color = NA) +
      geom_ribbon(aes(x = mod_sum$yr.date, ymin = pc25, ymax = pc75),
                  fill = "grey",
                  color = NA) +
      geom_line(aes(x = mod_sum$yr.date, y = mean), lwd = 2) 
    trend.plot <- p2 +
      xlab("Year") +
      ylab("Trend (proportional change per year)") +
      geom_abline(intercept = 0, slope = 0) +
      theme_bw()
    
    grid.arrange(count.plot, trend.plot, heights = c(1.5, 1))
    
  }
  
} 

# plot covariate effects
plot_covars <- function(mod_sum, bugsdata, resp, cov_names = NULL) {
  
  cov_mean <- matrix(mod_sum$cov_plot_vals[, "mean"], ncol = bugsdata$Q, byrow = TRUE)
  cov_lower <- matrix(mod_sum$cov_plot_vals[, "2.5%"], ncol = bugsdata$Q, byrow = TRUE)
  cov_upper <- matrix(mod_sum$cov_plot_vals[, "97.5%"], ncol = bugsdata$Q, byrow = TRUE)
  
  if (is.null(cov_names))
    cov_names <- letters[seq_len(bugsdata$Q)]
  
  oldmfrow <- par()$mfrow
  par(mfrow = c(2, 2))
  for (i in seq_len(bugsdata$Q)) {
    plot(cov_mean[, i] ~ bugsdata$Xplot[2:nrow(bugsdata$Xplot), i],
         type = "l", bty = "l", las = 1,
         xlab = cov_names[i], ylab = paste0("Effect on ", resp),
         ylim = range(c(cov_mean, cov_lower, cov_upper)))
    polygon(c(bugsdata$Xplot[2:nrow(bugsdata$Xplot), i], bugsdata$Xplot[nrow(bugsdata$Xplot):2, i]),
            c(cov_lower[, i], rev(cov_upper[, i])),
            border = NA,
            col = "gray50")
    lines(cov_mean[, i] ~ bugsdata$Xplot[2:nrow(bugsdata$Xplot), i])
    lines(range(bugsdata$Xplot[2:nrow(bugsdata$Xplot), i]), c(0, 0), lty = 2)
  } 
  
  par(mfrow = oldmfrow)
  
}
