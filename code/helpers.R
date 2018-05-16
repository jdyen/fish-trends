# helpers

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
  
  dat$yadj <- y
  dat$qyear <- dmy(paste0("01/01/", dat$year))
  yr <- sort(unique(dat$year))
  yr.date <- dmy(paste0("01/06/", yr))
  
  real_fitted <- fit$summary[grep("y.sim", rownames(fit$summary)), ]
  r2_tmp <- cor(real_fitted[, "mean"], y) ** 2
   
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
    cov_res <- fit$mean$inc.cov
    if (!is.null(colnames(bugsdata$Xcov)))
      names(cov_res) <- colnames(bugsdata$Xcov)
    cov_prior <- 1 - bugsdata$pc[1]
    or_cov <- cov_res * (1 - cov_prior) / ((1 - cov_res) * cov_prior)
    cov_inc <- cov_res
    cov_or <- or_cov
    cov_plot_vals <- fit$summary[grep("cov.plot\\[", rownames(fit$summary)), ]
  } else { 
    cov_res <- rep(NA, 2)
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
              cov_plot_vals = cov_plot_vals)
  
  out
  
} 

# plot fitted models
plot_fitted <- function(mod_sum, bugsdata, dat, sp_names, spp, resp) {
  
  for (sysnum in seq_len(bugsdata$Nsystem + 1)) {
    
    fitted.r <- mod_sum$fitted[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$fitted)), ]
    trends.r <- mod_sum$trends[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$trends)), ]
    mean.trends.r <- mod_sum$mean.trends[grep(paste0("\\[", sysnum, ","), rownames(mod_sum$mean.trends)), ]
    dat.r <- dat[system == sysnum, ]
    if (length(unique(dat.r$system)) > 1)
      stop("check data, >1 system code in subset")
    if(sysnum == (bugsdata$Nsystem + 1))
      dat.r <- dat
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
                 aes(x = mod_sum$date, y = yadj),
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
plot_covars <- function(mod_sum, bugsdata, cov_names) {
  
  cov_mean <- matrix(mod_sum$cov_plot_vals[, "mean"], ncol = bugsdata$Q, byrow = TRUE)
  cov_lower <- matrix(mod_sum$cov_plot_vals[, "2.5%"], ncol = bugsdata$Q, byrow = TRUE)
  cov_upper <- matrix(mod_sum$cov_plot_vals[, "97.5%"], ncol = bugsdata$Q, byrow = TRUE)
  
  oldmfrow <- par()$mfrow
  par(mfrow = c(2, 2))
  for (i in seq_len(bugsdata$Q)) {
    plot(cov_mean[, i] ~ bugsdata$Xplot[2:nrow(bugsdata$Xplot), i],
         type = "l", bty = "l", las = 1,
         xlab = cov_name[i], ylab = paste0("Effect on ", resp),
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
