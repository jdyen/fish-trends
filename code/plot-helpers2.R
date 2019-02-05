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