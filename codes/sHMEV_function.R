

data_nyears = function(data, nyears) {
  "---------------------------------------------------------------------
  from a dataframe data produce a new dataframe obtained taking the 
  first nyears years of observations for each station
  ----------------------------------------------------------------------"
  year_split = data %>% group_by(STATION) %>%
    dplyr::summarise(year_un = unique(YEAR)) %>%
    summarise(nend = nth(year_un, nyears))
  n_end = year_split$nend
  d = NULL
  n_station = length(unique(data$STATION))
  id_station = unique(data$STATION)
  for (i in 1:n_station) {
    d1 = data %>% filter(STATION == id_station[i] & YEAR <= n_end[i])
    d = rbind.data.frame(d, d1)
  }
  return(d)
}


gridPoints = function(name_state, cell_size = 0.1) {
  "---------------------------------------------------------------------
  cretate a grid of points with size of the cellse equal to cell_size 
  for the state name_state of USA.
  It also compute elevation for each point of the grid.
  ----------------------------------------------------------------------"
  us = getData('GADM', country = 'US', level = 1)
  nc = us[us$NAME_1 == name_state, ]
  
  # Create a grid of points within the bbox of the SpatialPolygonsDataFrame
  grid = makegrid(nc, cellsize = cell_size) # cellsize in map units
  grid = SpatialPoints(grid, proj4string = CRS(proj4string(nc)))
  grid = grid[nc,] # subset only point inside nc
  g = grid@coords
  g = data.frame(g)
  colnames(g) = c("x", "y") # x=long, y=lat
  
  # elevation for each point of the grid
  elev_grid = get_elev_point(g, prj = "EPSG:4326")
  elevation = elev_grid@data$elevation
  g = cbind.data.frame(g, elevation)
  g
}


shmev_ppd = function(model_fit, id_st, nyear_sim = 20) {
  "---------------------------------------------------------------------
   Compute ppd for yearly maxima for station with id equal to id_st

   arguments:
   model_fit -> stan object resulting from the fit of the sHMEV model
   id_st -> id of the station for which you want to compute ppd

   output: list with array of simulated annual maxima 
   and of daily data xij 
   (the nyear_sim is the number of years generated in the model)
  ----------------------------------------------------------------------"
  sim = rstan::extract(model_fit)
  ndraws = dim(sim$Ngen)[1]
  nyears = dim(sim$Ngen)[3]
  max = matrix(0, nrow = ndraws, ncol = nyears)
  daily_data = matrix(0, nrow = ndraws, ncol = nyears * 366)
  
  for (s in 1:ndraws) {
    daily_data_s = matrix(0, nrow = nyears, ncol = 366)
    for (j in 1:nyears) {
      if (sim$Ngen[s, id_st, j] > 0) {
        xij = rweibull(sim$Ngen[s, id_st, j], sim$ggen[s, id_st, j],
                       sim$dgen[s, id_st, j])
        # generate a number of values equal to the number of days with PRCP > 0
        daily_data_s[j, 1:sim$Ngen[s, id_st, j]] = xij
        max[s, j] <- max(xij)
      }
      daily_data[s,] = as.vector(daily_data_s)
    }
  }
  return(list(max = max, xij = daily_data))
}



plot_ppd = function(model_fit, id_st, data, ndrawsplot = 100) {
  "---------------------------------------------------------------------
  plot posterior predictive distributions for N, maxima and daily 
  data xij
  
  arguments:
  ndraws2plot -> number of draws to limit at for plotting only
  data -> observed data 
  id_st -> id of the station for which you want the plot
  ----------------------------------------------------------------------"
  
  id_station = unique(data$STATION)
  simdata = shmev_ppd(model_fit, id_st)
  sim = rstan::extract(model_fit)
  nyears = dim(sim$gamma)[3] # years used to fit the model
  Nt = 366
  d_train = data_nyears(data, nyears)
  d_train = d_train %>% filter(STATION == id_station[id_st])
  d_matrix = array(NA, dim = c(nyears, Nt))
  d_matrix = table_max(d_train, Nt = Nt)$data
  maxima = apply(d_matrix , 1, max)
  exceedances = pmax(d_matrix, 0)
  excesses = d_matrix[d_matrix > 0]
  N = rep(0, nyears)
  for (i in 1:nyears) {
    samplei = exceedances[i, ]
    N[i] = length(samplei[samplei > 0])
  }
  
  fig1 = ppc_dens_overlay(log(maxima), log(simdata$max[1:ndrawsplot,
                          1:nyears])) +
    labs(x = TeX(" $\\log{(max(x_{ij}))}$"))

  fig2 = ppc_dens_overlay(N, sim$Ngen[1:ndrawsplot, id_st, 1:nyears]) +
    labs(x = TeX("$n_{j}$"))
  
  fig3 <- ggplot()
  for (s in 1:ndrawsplot) {
    sample_s = simdata$xij[s,]
    wets_s = sample_s[sample_s > 0]
    fig3 = fig3 + 
      geom_line( aes_string(log(wets_s)), stat = 'density', adjust = 1 / 2,
      alpha = 0.3, col = 'steelblue2', lwd = 0.2)
  }
  fig3 = fig3 + 
    geom_line(aes(log(excesses)), stat = 'density', adjust = 1, col = 'black',
    lwd = 1, alpha = 0.8)
  
  fig3 = fig3 + scale_x_continuous(limits = c(-4, 6))
  fig3 = fig3 + labs(y = "",  x = TeX("$\\log{(x_{ij}) }$"))
  fig3 = fig3 + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.text.y = element_blank(),
    text = element_text(size = 21), axis.ticks = element_blank())
  
  fig2 = fig2 + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"), axis.text.y = element_blank(),
    text = element_text(size = 21), axis.ticks = element_blank())
  
  fig1 = fig1 + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.text.y = element_blank(),
    text = element_text(size = 21), axis.ticks = element_blank())
  
  fig <- grid.arrange(fig1, fig2, fig3, ncol = 3)
  # fig1 = fig1 + annotate("text", x=3.45, y= 1.8*0.85, label="a)",size = 10)
  # fig2 = fig2 + annotate("text", x=75, y= 0.1, label="b)", size = 10)
  # fig3 = fig3 + annotate("text", x=-3.6, y=0.87*0.4, label="c)", size = 10)
  # fig1 = fig1 + annotate("text", x=3.45, y= 1.8*0.95, label="a)",size = 10)
  # fig2 = fig2 + annotate("text", x=95, y= 0.072, label="b)", size = 10)
  # fig3 = fig3 + annotate("text", x=-3.6, y=0.83*0.4, label="c)", size = 10)
  return(fig)
}



ppd_map <- function(model_fit, grid) {
  "---------------------------------------------------------------------
  compute the posterior predictive distributions for the Weibull
  parameters (gamma and delta) and for the binomial success probability 
  (lambda) for each point of the grid grid
  
  arguments:
  model_fit -> stan object resulting from the estimation of the sHMEV
  model
  grid -> grid of points obtained from the function gridPoints
  
  output: mean and IQR of posterior distribution of lambda, delta and 
  gamma for each point of the grid
  ----------------------------------------------------------------------"
  sim = rstan::extract(model_fit)
  npoints = nrow(grid)
  
  # posterior estimation for regression coefficients
  alpha  = sim$ag
  beta = sim$bl
  nu = sim$nd
  sdl = sim$sdl
  sg = sim$sg
  B = dim(sim$lp__)
  
  mu.g = mu.d =  lambda = matrix(NA, ncol = npoints , nrow = B)
  # on each row there are the values for one MC iteration for grid's points
  int = rep(1, npoints)
  xx = cbind(int, grid$y, grid$x, grid$elevation, grid$dist_coast)
  xxs = xx
  xxs[, c(2:5)] = scale(xx[, 2:5])
  
  for (i in 1:B) {
    mu.g[i, ] = xxs %*% alpha[i, ]
    mu.d[i, ] = xxs %*% nu[i, ]
    pl = xxs %*% beta[i, ]
    lambda[i, ] = exp(pl) / (1 + exp(pl))
  }
  
  gm = dl = N = matrix(NA, nrow = B, ncol = npoints)
  for (b in 1:B) {
    for (i in 1:npoints) {
      # simulate Mgen times gamma, delta, N
      ds = extraDistr::rgumbel(n = 50, mu.d[b, i], sdl[b])
      gs = extraDistr::rgumbel(n = 50, mu.g[b, i], sg[b])
      Ns = rbinom(50, size = 366, prob = lambda[b, i])
      gm[b, i] = mean(gs)
      dl[b, i] = mean(ds)
      N[b, i] = mean(Ns)
    }
  }
  gmean = apply(gm, 2, mean)
  dmean = apply(dl, 2, mean)
  Nmean = apply(N, 2, mean)
  giqr = apply(gm, 2, IQR)
  diqr = apply(dl, 2, IQR)
  Niqr = apply(N, 2, IQR)
  mmu.g = apply(mu.g, 2, mean)
  iqrmu.g = apply(mu.g, 2, IQR)
  mmu.d = apply(mu.d, 2, mean)
  iqrmu.d = apply(mu.d, 2, IQR)
  mlambda = apply(lambda, 2, mean)
  iqrlambda = apply(lambda, 2, IQR)
  d_post = cbind.data.frame(grid, gmean, giqr, dmean,diqr, Nmean, Niqr, mmu.g,
                            iqrmu.g, mmu.d, iqrmu.d, mlambda, iqrlambda)
  return(d_post)
}


wei_cdf <- function(x, d = d, g = g, N = N) {
  'cdf of maxima of a Weibull'
  # only for scalar x
  cdf = mean(pweibull(x, g, d) ^ N)
  return(cdf)
}


quant_map <- function(model_fit, grid, Tr = 25) {
  "---------------------------------------------------------------------
  compute the predictive Tr-years return values for each point of 
  the grid grid
  
  arguments:
  model_fit -> stan object resulting from the estimation of the sHMEV 
  model
  grid -> grid of points obtained from the function gridPoints
  Tr -> return time (years)
  
  output: pointwise mean and pointwise credible intervals of return 
  values
  ----------------------------------------------------------------------"
  sim = rstan::extract(model_fit)
  npoints = nrow(grid)
  
  # posterior estimation for regression coefficients
  alpha  = sim$ag
  beta = sim$bl
  nu = sim$nd
  sdl = sim$sdl
  sg = sim$sg
  B = dim(sim$lp__)
  
  mu.g = mu.d =  lambda = matrix(NA, ncol = npoints , nrow = B)
  # on each row there are the values for one MC iteration for grid's points
  int = rep(1, npoints)
  xx = cbind(int, grid$y, grid$x, grid$elevation, grid$dist_coast)
  xxs = xx
  xxs[, c(2:5)] = scale(xx[, 2:5])
  for (i in 1:nrow(mu.g)) {
    mu.g[i, ] = xxs %*% alpha[i, ]
    mu.d[i, ] = xxs %*% nu[i, ]
    pl = xxs %*% beta[i, ]
    lambda[i, ] = exp(pl) / (1 + exp(pl))
  }
  
  nwarning_x0 = 0 # if true there are numerical problems with optimization
  Fival = 1 - 1 / Tr
  quants = matrix(NA, nrow = B, ncol = npoints) # return values
  
  for (b in 1:B) {
    for (i in 1:npoints) {
      # simulate Mgen times gamma, delta, N
      ds = extraDistr::rgumbel(n = 50, mu.d[b, i], sdl[b])
      gs = extraDistr::rgumbel(n = 50, mu.g[b, i], sg[b])
      Ns = rbinom(50, size = 366, prob = lambda[b, i])
      myfun = function(x)
        wei_cdf(x, d = ds, g = gs, N = Ns) - Fival
      F0 <- 0.1 # Tr 10 years to start 
      x0 <-  mean(ds) * (log(mean(Ns) / (1 - F0))) ^ (1 / mean(gs))
      optim = nleqslv(x0, myfun)
      # solve the equation F(y) = p_o
      quants[b, i] = optim$x
      termcd = optim$termcd
      if (termcd != 1) {
        print("comp_quant WARNING: nleqslv might not be converging")
        print("If 6 Jacobian is singular, new guess")
        print("termcd = ")
        print(termcd)
        nwarning_x0 = nwarning_x0 + 1
      }
    }
  }
  # compute uncertainty bands and average / mean values
  qmean =  apply(quants, 2, mean)
  qquant = apply(quants, 2, quantile, prob = c(0.05, 0.5, 0.95))

  return( list(qmean = qmean, qupper = qquant[3, ], qlower = qquant[1,], 
         qmedian = qquant[2,], quants = quants, Fival = Fival, Tr = Tr,
      nwarning_x0 = nwarning_x0))
}


quantile_validation = function(model_fit, maxval, covariate, id_st, trmin = 2) {
  "---------------------------------------------------------------------
  spatial validation for stations in the test set (not used to estimate
  the model)
  
  arguments:
  model_fit -> stan object resulting from the estimation of the sHMEV 
  model
  maxval -> maxima of all records
  covariate -> covariates of the test set
  id_st -> id of the station on the test set for which you want the
  validation
  trmin -> minimum return time for computing quantiles
  
  output: pointwise mean and pointwise credible intervals of return 
  values
  ----------------------------------------------------------------------"
    sim = rstan::extract(model_fit)
    Mgen = 50
    
    # posterior estimation for regression coefficients
    alpha  = sim$ag
    beta = sim$bl
    nu = sim$nd
    sdl = sim$sdl
    sg = sim$sg
    
    mu.g = mu.d =  lambda = rep(NA, nrow(beta))
    for (i in 1:length(mu.g)) {
      mu.g[i] = covariate[id_st, ] %*% alpha[i, ]
      mu.d[i] = covariate[id_st, ] %*% nu[i, ]
      pl = covariate[id_st, ] %*% beta[i, ]
      lambda[i] = exp(pl) / (1 + exp(pl))
    }
    
    nwarning_x0 = 0 # if true there are numerical problems with optimization
    Np = length(maxval) # number of years for the validation
    Fival = (1:Np) / (1 + Np) 
    # Weibull plotting position non exceedance frequency 
    Xival = sort(maxval)    # sorted maxima
    Trval = 1 / (1 - Fival)     # return time empirical values
    
    # compute quantiles only for Tr > trmin
    TrvalQ = Trval[Trval >= trmin]
    FivalQ = Fival[Trval >= trmin]
    XivalQ = Xival[Trval >= trmin]
    NpQ = length(XivalQ)
    
    B = dim(sim$lp__) # number of draws from the posterior
    quants = matrix(0, nrow = B, ncol =  NpQ) # return values
    
    for (b in 1:B) {
      # simulate Mgen times gamma, delta, N
      ds = extraDistr::rgumbel(n = 50, mu.d[b], sdl[b])
      gs = extraDistr::rgumbel(n = 50, mu.g[b], sg[b])
      Ns = rbinom(50, size = 366, prob = lambda[b])
      
      for (i in 1:NpQ) {
        # compute return values for Npq associated return times
        myfun = function(x)
          wei_cdf(x, d = ds, g = gs, N = Ns) - FivalQ[i]
        F0 <- 0.1 # Tr 10 years to start 
        x0 <-  mean(ds) * (log(mean(Ns) / (1 - F0))) ^ (1 / mean(gs))
        optim = nleqslv(x0, myfun)
        # solve the equation F(y) = p_o
        quants[b, i] = optim$x
        termcd = optim$termcd
        if (termcd != 1) {
          print("comp_quant WARNING: nleqslv might not be converging")
          print("If 6 Jacobian is singular, new guess")
          print("termcd = ")
          print(termcd)
          nwarning_x0 = nwarning_x0 + 1
        }
      }
    }
    # compute uncertainty bands and average / mean values
    qmean =  apply(quants, 2, mean)
    qquant = apply(quants, 2, quantile, prob = c(0.05, 0.5, 0.95))
    
    return(list(qmean = qmean, qupper = qquant[3, ], qlower = qquant[1,],
                qmedian = qquant[2,], quants = quants, Fival = Fival,
                Xival = Xival, Trval = Trval, FivalQ = FivalQ, XivalQ = XivalQ,
                TrvalQ = TrvalQ, Mgen = Mgen, nwarning_x0 = nwarning_x0))
}


quant_spatial <- function(model_fit, maxval, trmin = 2, id_st) {
  "---------------------------------------------------------------------
  Compute quantiles and goodness of fit measures for a fitted sHMEV model 
  
  arguments:
  model_fit -> stan object resulting from the fit of the sHMEV model
  maxval -> maxima of all records
  id_st -> id of the station for which you want the computation
  trmin -> minimum return time for computing quantiles
  
  output: pointwise mean and pointwise credible intervals of return 
  values, Fractional square error (FSE), Mean bias (mbias) and 
  mean credibility interval width (mwidth)
  ----------------------------------------------------------------------"
  model_sim = rstan::extract(model_fit)
  Mgen = dim(model_sim$Ngen)[3]
  Nt = 366
  
  nwarning_x0 = 0 # if true there are numerical problems with optimization
  Np = length(maxval) # number of total  years 
  Fival = (1:Np) / (1 + Np) 
  # Weibull plotting position non exceedance frequency 
  Xival = sort(maxval) # sorted maxima
  Trval = 1 / (1 - Fival) # return time empirical values
  
  # for quantiles only, compute only for Tr > trmin
  TrvalQ = Trval[Trval >= trmin]
  FivalQ = Fival[Trval >= trmin]
  XivalQ = Xival[Trval >= trmin]
  NpQ = length(XivalQ)
  
  B = dim(model_sim$lp__) # number of draws from the posterior
  quants = matrix(0, nrow = B, ncol =  NpQ) # return values
  epsi   = matrix(0, nrow = B, ncol =  NpQ) # needed to compute FSE 
  
  for (b in 1:B) {
    # print(b)
    ds = model_sim$dgen[b, id_st, ] # delta_b of posterior
    gs = model_sim$ggen[b, id_st, ] # gamma_b of posterior
    Ns = model_sim$Ngen[b, id_st, ] # N_b of posterior
    for (i in 1:NpQ) {
      # compute return values for Npq associated return times
      myfun = function(x)
        wei_cdf(x, d = ds, g = gs, N = Ns) - FivalQ[i]
      F0 <- 0.1 # Tr 10 years to start 
      x0 <-  mean(ds) * (log(mean(Ns) / (1 - F0))) ^ (1 / mean(gs))
      optim = nleqslv(x0, myfun)
      # solve the equation F(y) = p_o
      quants[b, i] = optim$x
      termcd = optim$termcd
      if (termcd != 1) {
        print("comp_quant WARNING: nleqslv might not be converging")
        print("If 6 Jacobian is singular, new guess")
        print("termcd = ")
        print(termcd)
        nwarning_x0 = nwarning_x0 + 1
      }
      epsi[b, i]   = (quants[b, i] - XivalQ[i]) / XivalQ[i]
    }
  }
  
  # compute uncertainty bands and average / mean values
  qmean =  apply(quants, 2, mean) 
  qquant = apply(quants, 2, quantile, prob = c(0.05, 0.5, 0.95))
  
  fse_Tr = rep(NA, NpQ) # FSE for a given non exceedance frequency
  bias_Tr = rep(NA, NpQ) # bias
  width_Tr = rep(NA, NpQ) # width of credibility interval
  for (i in 1:NpQ) {
    # averaged over the S draws from the posterior
    fse_Tr[i]  = sqrt(mean(epsi[, i] ^ 2))
    bias_Tr[i] = mean(epsi[, i])
    width_Tr[i] = qquant[3, i] - qquant[1, i]
  }
  fse = mean(fse_Tr)
  mbias = mean(bias_Tr)
  mwidth = mean(width_Tr)
  
  return(list(qmean = qmean, qupper = qquant[3, ], qlower = qquant[1,], 
              qmedian = qquant[2,], quants = quants, Fival = Fival, 
              Xival = Xival, Trval = Trval, FivalQ = FivalQ, XivalQ = XivalQ, 
              TrvalQ = TrvalQ, fse = fse, fse_Tr = fse_Tr, mbias = mbias, 
              bias_Tr = bias_Tr, mwidth = mwidth, width_Tr = width_Tr,
              Mgen = Mgen, nwarning_x0 = nwarning_x0))
}




simul_data = function(J, S, ptrue, ntrue, dist, Nt = 366, seed = "123") {
  "---------------------------------------------------------------------
  Generate syntetic data according to a given specification. 
  
  arguments:
  S -> number of stations in the data to be generated
  J -> length of the time series to be generated (in years / blocks)
  ptrue -> list of parameters describing the distribution of event
  magnitudes
  ntrue -> vector of dimensione S describing the success probability 
  for the binomial distribution for the number of event 
  dist -> distribution of the event magnitudes, can be 'wei', 'gam'
  or 'wei_gp'
  Nt -> the length of each block

  output: generates an array xij of synthetic data with dimension S*J*Nt 
  and non-zero events are stacked at the 'beginning' of each block 
  ----------------------------------------------------------------------"
  
  set.seed(seed)
  nj = matrix(0, nrow = S, ncol = J)
  for (s in 1:S) {
    nj[s, ] = rbinom(J, Nt, ntrue[s])
  }
  dailyp = array(0, dim = c(S, J, Nt))
  
  if (dist == "wei") {
    for (s in 1:S) {
      for (j in 1:J) {
        Cj = extraDistr::rgumbel(1, mu = ptrue$md[s], sigma = ptrue$sd)
        wj = extraDistr::rgumbel(1, mu = ptrue$mg[s], sigma = ptrue$sg)
        dailyp[s, j, 1:nj[s, j]] = rweibull(nj[s, j], wj, Cj)
      }
    }
  }
  else if (dist == "gam") {
    for (s in 1:S) {
      for (j in 1:J) {
        dailyp[s, j, 1:nj[s, j]] = rgamma(nj[s, j], ptrue$a[s], ptrue$b[s])
      }
    }
  }
  else if (dist == "wei_gp") {
    cov_g = exp_cov_matrix(ptrue$dst, sill = ptrue$sill_g, 
                           range = ptrue$range_g, nugget = ptrue$nugget_g)
    cov_d = exp_cov_matrix(ptrue$dst, sill = ptrue$sill_d, 
                           range = ptrue$range_d, nugget = ptrue$nugget_d)
    require(mvtnorm)
    eps_g = mvtnorm::rmvnorm(1, sigma = cov_g)
    eps_d = mvtnorm::rmvnorm(1, sigma = cov_d)
    
    for (s in 1:S) {
      for (j in 1:J) {
        Cj = extraDistr::rgumbel(1, mu = ptrue$md[s] + eps_d[s],
                                 sigma = ptrue$sd)
        wj = extraDistr::rgumbel(1, mu = ptrue$mg[s] + eps_g[s],
                                 sigma = ptrue$sg)
        dailyp[s, j, 1:nj[s, j]] = rweibull(nj[s, j], wj, Cj)
      }
    }
  }
  else {
    print('ERROR: invalid dist')
  }
  
  maxima = matrix(0, nrow = S, ncol = J)
  for (s in 1:S) {
    for (j in 1:J) {
      maxima[s, j] = max(dailyp[s, j, ], na.rm = TRUE)
    }
  }
  Fi = (1:J) / (1 + J)
  Tr = 1 / (1 - Fi)
  Xi = matrix(0, nrow = S, ncol = J)
  for (s in 1:S) {
    Xi[s, ] = sort(maxima[s, ])
  }
  return(list(
    dailyp = dailyp, N = nj, maxima = maxima, S = S, nyears = J, Tr = Tr, 
    Fi = Fi, Xi = Xi))
}


exp_cov_matrix = function (x, sill, range, nugget) {
  'exponential covariance matrix of a gaussian process
  with x matrix of distance'
  S = nrow(x) # number of sites
  neg_inv_range = -1.0 / range
  cov = matrix(NA, nrow = S, ncol = S)
  for (i in 1:(S - 1)) {
    for (j in (i + 1):S) {
      cov[i, j] = sill * exp(neg_inv_range * x[i, j])
      cov[j, i] = cov[i, j]
    }
  }
  for (k in 1:S) {
    cov[k, k] = sill + nugget
  }
  return (cov)
}


qshmev = function(Tr, ptrue, ntrue, id_st = i, Nt = 366, nsamples = 50) {
  'compute theoretical shmev quantiles for the values in Fi'
   # values in Tr must be  return times
  Fi = 1 - 1 / Tr
  Ngen = rbinom(nsamples, Nt, ntrue[i])
  dgen = extraDistr::rgumbel(nsamples, mu = ptrue$md[i], sigma = ptrue$sd)
  ggen = extraDistr::rgumbel(nsamples, mu = ptrue$mg[i], sigma = ptrue$sg)
  Nx = length(Fi)
  res = rep(0, Nx)
  for (i in 1:Nx) {
    res[i] <- hbev_wei_quant(Fi[i], C = dgen, W = ggen, N = Ngen)
  }
    return(res)
}
