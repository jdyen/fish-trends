make.model.file.hier <- function(filename,
                                 covar = FALSE,
                                 npred,
                                 nbreak,
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


for (r in 1:Nsystem) {
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

for (r in 1:Nsystem) {
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
  
  cat(paste0("  for (i in 1:", nbreak + 1, ") {\n"), fill = TRUE)
  for (i in 1:npred) {
      cat(paste0("    x", i, "[i] <- xdum[i, ", i, "]\n"), fill = TRUE)
  }
  cat("  }", fill = TRUE)
  
  cat(paste0("  for (i in 1:", nbreak, ") {\n"), fill = TRUE)
  for (i in 1:npred) {
    cat(paste0("    beta.cov.all[i, ", i, "] <- beta.cov", i, "[i + 1] - beta.cov", i, "[i]\n"))
  }
  cat("  }\n", fill = TRUE)
  
  for (i in 1:npred) {
    cat(paste0("  beta.cov", i, "[1:", nbreak + 1,
               "] <- jump.pw.poly.df.gen(x", i, "[1:",
               nbreak + 1, "], kcov[",
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
