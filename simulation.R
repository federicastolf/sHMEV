
rm(list = ls())

library(tidyverse)
library(hmevr)
library(extraDistr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(nleqslv)
library(gridExtra)

source("sHMEV_function.R")


#----------------------# Simulate synthetic datasets #-------------------------#

## analysis parameter
S = 27 # number of sites
ny_train = 20 # number of years in the training set
ny_test = 100 # number of years in the test set
Mgen = 50 # number of sample points for sHMEV latent parameters
Nt = 366 # number of events / block = number of days/year

lat = runif(S)
long = runif(S)
covariate = as.matrix(cbind(rep(1,S), lat, long))

# parameters for generate n_j (number of events)
coeff_true0 = list(l0 = -0.7, l1 = 0.15, l2 = 0.25)
# linear predictor for lambda
pll = coeff_true0$l0 + coeff_true0$l1*lat + coeff_true0$l2*long 
# success probability for the binomial distribution
lambda = exp(pll)/(1+exp(pll)) 


#--- scenario 1 (WEI) ---#

coeff_true1 = list(d0 = 9, d1 = -0.4, d2 = 0.8, g0 = 0.85, g1 = 0.15, g2 = 0.05)
coeff_true = append(coeff_true0, coeff_true1)
ptrue = list(md = coeff_true$d0 + coeff_true$d1*lat + coeff_true$d2*long,
             # gumbel location parameter for delta
             mg = coeff_true$g0 + coeff_true$g1*lat + coeff_true$g2*long, 
             # gumbel location parameter for gamma
             sd = 1, # gumbel scale parameter for delta
             sg = 0.1) # gumbel scale parameter for gamma

dtsim = simul_data(ny_train, S, ptrue, lambda, dist = "wei")
# independent dataset for validation purposes
dvsim = simul_data(ny_test, S, ptrue, lambda, dist = "wei")


#--- scenario 2 (GAM) ---#

ptrue_gam = list(a = 0.9 + 0.15*lat + 0.5*long,
                 b = 0.15 - 0.1*lat + 0.05*long)

dtsim_gam = simul_data(ny_train, S, ptrue_gam, lambda, dist = "gam")
# independent dataset for validation purposes
dvsim_gam = simul_data(ny_test, S, ptrue_gam, lambda, dist = "gam")


#--- scenario 3 (WEI_gp) ---#

ptrue_gp = list(dst = as.matrix(dist(cbind(lat,long))), # distance matrix
             sill_g = 0.005, # sill for exponential covariance for gamma
             range_g = 22, # range for exponential covariance for gamma
             nugget_g = 0.001, # nugget for exponential covariance for gamma
             sill_d = 0.4, # sill for exponential covariance  for delta
             range_d = 22, # range for exponential covariance for delta
             nugget_d = 0.001, # nugget for exponential covariance for delta
             md =  9  - 0.4*lat + 0.8*long,
             # gumbel location parameter for delta
             mg = 0.85 + 0.15*lat + 0.05*long, 
             # gumbel location parameter for gamma
             sd = 1, # gumbel scale parameter for delta
             sg = 0.1) # gumbel scale parameter for gamma

dtsim_gp = simul_data(ny_train, S, ptrue_gp, lambda, dist = "wei_gp")
# independent dataset for validation purposes
dvsim_gp = simul_data(ny_test, S, ptrue_gp, lambda, dist="wei_gp")


#----------------------------# Fit sHMEV model #-------------------------------#

# model data for the three scenarios
model_data = list(M = dtsim$nyears, Nt = Nt, y = dtsim$dailyp, N = dtsim$N,
                  S = dtsim$S, X = covariate, Mgen = Mgen)
model_data.gam0 = list(M = dtsim_gam$nyears, Nt = Nt, y = dtsim_gam$dailyp, 
                       N = dtsim_gam$N, S = dtsim_gam$S, X = covariate, 
                       Mgen = Mgen)
model_data.gp = list(M = dtsim_gp$nyears, Nt = Nt, y = dtsim_gp$dailyp, 
                     N = dtsim_gp$N, S = dtsim_gp$S, X = covariate, Mgen = Mgen)

# hyperparameters 
mod_ippars = list(sg0prior = c(10, 0.35), # a, b of inverse gamma for gamma
                  sd0prior = c(10, 19.75), #  a, b of inverse gamma for delta
                  ma0prior = c(0, 0, 0), # mean of the normal for gamma
                  sa0prior = c(0.5, 0.5, 0.5), # sd of the normal for gamma
                  mn0prior = c(0, 0,0), # mean of the normal for delta
                  sn0prior = c(0.5, 0.5, 0.5), # sd of the normal for delta
                  mb0prior = c(-0.5, 0, 0), # mean of the normal for lambda
                  sb0prior = c(1, 1, 1)) # sd of the normal for lambda

model_datas = append(model_data, mod_ippars)
model_data.gam = append(model_data.gam0, mod_ippars)
model_datas.gp = append(model_data.gp, mod_ippars)

spat.sim = rstan::stan(file="sHMEV_sim.stan", data = model_datas,
                       iter = 2000, chains = 4, 
                       control = list(adapt_delta = 0.8, max_treedepth = 10))

spat.sim.gam = rstan::stan(file="sHMEV_sim.stan", data = model_data.gam, 
                           iter = 2000, chains = 4, 
                           control = list(adapt_delta = 0.8, max_treedepth = 10))

spat.sim.gp = rstan::stan(file="sHMEV_sim.stan", data = model_datas.gp, 
                          iter = 2000, chains = 4, 
                          control = list(adapt_delta = 0.8, max_treedepth = 10))


#--- check convergence of the chains ---#

# sigma_gamma and sigma_delta 
parss = c("sg", "sdl")
stan_trace(spat.sim, pars = parss) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gam, pars = parss) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gp, pars = parss) +
  theme(plot.title = element_text(hjust = 0.5)) 

# regression coefficients for lambda
parsb = c("bl[1]", "bl[2]", "bl[3]") 
stan_trace(spat.sim, pars = parsb) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gam, pars = parsb) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gp, pars = parsb) +
  theme(plot.title = element_text(hjust = 0.5)) 

# regression coefficients for mu_delta
parsn = c("nd[1]", "nd[2]", "nd[3]") 
stan_trace(spat.sim, pars = parsn) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gam, pars = parsn) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gp, pars = parsn) +
  theme(plot.title = element_text(hjust = 0.5)) 

# regression coefficients for mu_gamma
parsa = c("ag[1]", "ag[2]","ag[3]") 
stan_trace(spat.sim, pars = parsa) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gam, pars = parsa) +
  theme(plot.title = element_text(hjust = 0.5)) 
stan_trace(spat.sim.gp, pars = parsa) +
  theme(plot.title = element_text(hjust = 0.5)) 


#--------------------------# results scenario WEI #----------------------------#

d = as.data.frame(spat.sim)
par_names = c(parss, parsb, parsn, parsa)
dd = d %>% dplyr::select(par_names)
colnames(dd)[3:11] = c("bl0", "bl1", "bl2", "bd0", "bd1", "bd2",
                       "bg0", "bg1", "bg2") 

# credible intervals of level 95%
ic = t(apply(dd, 2, quantile, prob = c(0.025, 0.975)))
colnames(ic) = c("lower", "upper")

# true values of the parameters used to simulate the data
par_true = c(ptrue$sg, ptrue$sd, coeff_true$l0, coeff_true$l1, coeff_true$l2,
             coeff_true$d0, coeff_true$d1, coeff_true$d2, coeff_true$g0, 
             coeff_true$g1, coeff_true$g2)
cbind(par_true,ic)

vlist = names(dd)
lab = c(expression(sigma[gamma]), expression(sigma[delta]), 
        expression(beta[lambda][","][0]), expression(beta[lambda][","][1]),
        expression(beta[lambda][","][2]), expression(beta[delta][","][0]),
        expression(beta[delta][","][1]), expression(beta[delta][","][2]),
        expression(beta[delta][","][2]), expression(beta[gamma][","][1]),
        expression(beta[gamma][","][2])) 

# plot of posterior distributions
list_of_ggplots = lapply(1:length(vlist), function(i) {
  p = ggplot(dd, aes_string(x = vlist[i])) + 
    geom_histogram(color = "black", fill = "gray") +
    geom_vline(aes(xintercept = par_true[i]), color = "red", size = 0.9) +
    geom_vline(aes(xintercept = ic[i,1]), color = "red", linetype = "dashed", 
               size = 0.7) +
    geom_vline(aes(xintercept = ic[i,2]), color = "red", linetype = "dashed", 
               size = 0.7) +
    xlab(lab[i]) + ylab("") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_bw() + 
    theme(text = element_text(size = 15), axis.text.x = element_text(size = 10)) 
  p
})

hist_plot = grid.arrange(grobs = list_of_ggplots, ncol = 4)
# ggsave("hist_sim.png", plot = hist_plot, device = 'png', width = 8,
# height = 7)


#-------------------# comparison with competing methods #----------------------#


#--- scenario 1 (WEI) ---#

fse_sim1 = mbias_sim1 = mwidth_sim1 = NULL

for(i in 1:S){
  
  print(paste("station", i))
  
  # HMEV 
  hmev = hmevr::fit_ev_model(dtsim$dailyp[i,,], model = 'wei_dgu_bin',  
                             iter = 2000, chains = 4, Mgen = 50)
  hmevq = hmevr::comp_quant(hmev, dvsim$maxima[i,], trmin = 2)
  
  #GEV
  gev = hmevr::fit_ev_model(dtsim$dailyp[i,,], model = 'gev', iter = 2000, 
                            chains = 4, Mgen = 50)
  gevq = hmevr::comp_quant(gev, dvsim$maxima[i,], trmin = 2)
  
  # sHMEV 
  shmev = quant_spatial(spat.sim, maxval = dvsim$maxima[i,], id_st = i)
  
  adt_f = cbind.data.frame(i,round(gevq$fse,3), round(hmevq$fse,3), 
                           round(shmev$fse,3))
  adt_w = cbind.data.frame(i,round(gevq$mwidth,3), round(hmevq$mwidth,3), 
                           round(shmev$mwidth,3))
  adt_b = cbind.data.frame(i,round(gevq$mbias,3), round(hmevq$mbias,3), 
                           round(shmev$mbias,3))
  colnames(adt_f) = colnames(adt_w) = colnames(adt_b) = c("station", "gev",
                      "hmev", "shmev")
  fse_sim1 = rbind(fse_sim1,adt_f)
  mwidth_sim1 = rbind(mwidth_sim1, adt_w)
  mbias_sim1 = rbind(mbias_sim1, adt_b)
}

#--- scenario 2 (GAM) ---#

fse_sim_gam = mbias_sim_gam = mwidth_sim_gam = NULL

for(i in 1:27){
  
  print(paste("station", i))
  
  # HMEV 
  dyn1 = hmevr::fit_ev_model(dtsim_gam$dailyp[i,,], model = 'wei_dgu_bin',  
                             iter = 2000, chains = 4, Mgen = 50)
  hmq = hmevr::comp_quant(dyn1, dvsim_gam$maxima[i,], trmin = 2)
  
  #GEV
  gev1 = hmevr::fit_ev_model(dtsim_gam$dailyp[i,,], model = 'gev', iter = 2000,
                             chains = 4, Mgen = 50)
  gevq = hmevr::comp_quant(gev1, dvsim_gam$maxima[i,], trmin = 2)
  
  # sHMEV 
  shmq = quant_spatial(spat.sim.gam, maxval = dvsim_gam$maxima[i,], id_st = i)
  
  adt_f = cbind.data.frame(i,round(gevq$fse,3), round(hmq$fse,3),
                           round(shmq$fse,3))
  adt_w = cbind.data.frame(i,round(gevq$mwidth,3), round(hmq$mwidth,3),
                           round(shmq$mwidth,3))
  adt_b = cbind.data.frame(i,round(gevq$mbias,3), round(hmq$mbias,3),
                           round(shmq$mbias,3))
  colnames(adt_f) = colnames(adt_w) = colnames(adt_b) =  c("station", "gev", 
                    "hmev", "shmev")
  fse_sim_gam = rbind(fse_sim_gam, adt_f)
  mwidth_sim_gam = rbind(mwidth_sim_gam, adt_w)
  mbias_sim_gam = rbind(mbias_sim_gam, adt_b)
}

#--- scenario 3 (WEI_gp) ---#

fse_sim_wgp = mbias_sim_wgp = mwidth_sim_wgp = NULL

for(i in 1:27){
  
  print(paste("station", i))
  
  # HMEV 
  dyn1 = hmevr::fit_ev_model(dtsim_gp$dailyp[i,,], model = 'wei_dgu_bin',  
                             iter = 2000, chains = 4, Mgen = 50)
  hmq = hmevr::comp_quant(dyn1, dvsim_gp$maxima[i,], trmin = 2)
  
  #GEV
  gev1 = hmevr::fit_ev_model(dtsim_gp$dailyp[i,,], model = 'gev', 
                             iter = 2000, chains = 4, Mgen = 50)
  gevq = hmevr::comp_quant(gev1, dvsim_gp$maxima[i,], trmin = 2)
  
  # sHMEV
  shmq = quant_spatial(spat.sim.gp, maxval = dvsim_gp$maxima[i,], id_st = i)
  
  adt_f = cbind.data.frame(i,round(gevq$fse,3), round(hmq$fse,3), 
                           round(shmq$fse,3))
  adt_w = cbind.data.frame(i,round(gevq$mwidth,3), round(hmq$mwidth,3), 
                           round(shmq$mwidth,3))
  adt_b = cbind.data.frame(i,round(gevq$mbias,3), round(hmq$mbias,3), 
                           round(shmq$mbias,3))
  colnames(adt_f) = colnames(adt_w) = colnames(adt_b) = c("station", "gev", 
                      "hmev", "shmev")
  fse_sim_wgp = rbind(fse_sim_wgp,adt_f)
  mwidth_sim_wgp = rbind(mwidth_sim_wgp, adt_w)
  mbias_sim_wgp = rbind(mbias_sim_wgp, adt_b)
}


#------ plot FSE, mwidth, mbias ------#

# FSE
fse_sim1$spec = rep("WEI", nrow(fse_sim1))
fse_sim_gam$spec = rep("GAM", nrow(fse_sim_gam))
fse_sim_wgp$spec = rep("WEIgp", nrow(fse_sim_wgp))
fse_sim = rbind(fse_sim1, fse_sim_gam, fse_sim_wgp)
fse_sim$spec = as.factor(fse_sim$spec)
colnames(fse_sim)[2:4] = c("GEV", "HMEV", "sHMEV")

fse_long = gather(fse_sim, "model", "FSE", 2:4)
fse_long$model = as.factor(fse_long$model)

fse_plot = ggplot(fse_long, aes(x = model, y = FSE, fill = model))+
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values=c("green", "blue", "red")) +
  facet_grid(.~ spec, scales = 'free') +
  theme_bw()  +
  theme(legend.position = 'none', text = element_text(size = 15),
         strip.background = element_rect(fill = "white"))

fse_plot
# ggsave('fse_sim.png', plot = fse_plot, device = 'png', width = 8, height = 5)

# mbias
mbias_sim1$spec = rep("WEI", nrow(mbias_sim1))
mbias_sim_gam$spec = rep("GAM", nrow(mbias_sim_gam))
mbias_sim_wgp$spec = rep("WEIgp", nrow(mbias_sim_wgp))
mbias_sim1 = rbind(mbias_sim1, mbias_sim_gam, mbias_sim_wgp)
mbias_sim1$spec = as.factor(mbias_sim1$spec)
colnames(mbias_sim1)[2:4] = c("GEV", "HMEV", "sHMEV")

mbias_long = gather(mbias_sim1, "model", "mbias", 2:4)
mbias_long$model = as.factor(mbias_long$model)

mbias_plot = ggplot(mbias_long, aes(x = model, y = mbias, fill = model))+
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values=c("green", "blue", "red")) +
  facet_grid(.~ spec, scales = 'free') +
  theme_bw()  +
  labs(y = expression(b[q]) )+
  theme(legend.position = 'none', text = element_text(size = 15),
         strip.background = element_rect(fill = "white"))

mbias_plot
# ggsave('mbias_sim.png', plot = mbias_plot, device = 'png', width = 8, 
# height = 5)

# mwidth
mwidth_sim1$spec = rep("WEI", nrow(mwidth_sim1))
mwidth_sim_gam$spec = rep("GAM", nrow(mwidth_sim_gam))
mwidth_sim_wgp$spec = rep("WEIgp", nrow(mwidth_sim_wgp))
mwidth_sim1 = rbind(mwidth_sim1, mwidth_sim_gam, mwidth_sim_wgp)
mwidth_sim1$spec = as.factor(mwidth_sim1$spec)
colnames(mwidth_sim1)[2:4] = c("GEV", "HMEV", "sHMEV")

mwidth_long = gather(mwidth_sim1, "model", "mwidth", 2:4)
mwidth_long$model = as.factor(mwidth_long$model)

mwidth_plot = ggplot(mwidth_long, aes(x = model, y = mwidth, fill = model))+
  geom_boxplot(  alpha = 0.6) +
  scale_fill_manual(values=c("green", "blue", "red")) +
  facet_grid(.~ spec, scales = 'free') +
  labs(y=expression(Delta[q[90]]))+
  theme_bw()  +
  theme(legend.position = 'none', text = element_text(size = 15),
        strip.background =element_rect(fill = "white"))

mwidth_plot
# ggsave('mwidth_sim.png', plot = mwidth_plot, device = 'png', width = 8,
# height = 5)


#------ quantile versus return time plots ------#


#-- station 1 example --#
i = 1

# HMEV 
dyn1 = hmevr::fit_ev_model(dtsim$dailyp[i,,], model = 'wei_dgu_bin',  
                           iter = 2000, chains = 4, Mgen = 50)
hmq = hmevr::comp_quant(dyn1, dvsim$maxima[i,], trmin = 2)

#GEV
gev1 = hmevr::fit_ev_model(dtsim$dailyp[i,,], model = 'gev', iter = 2000,
                           chains = 4, Mgen = 50)
gevq = hmevr::comp_quant(gev1, dvsim$maxima[i,], trmin = 2)

# sHMEV 
shmq = quant_spatial(spat.sim, maxval = dvsim$maxima[i,], id_st = i)

trtrue = c(2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 30, 40, 50, 
           75, 100)
qutrue = qshmev(trtrue, ptrue = ptrue, ntrue = lambda, id_st = i, Nt = Nt, 
                nsamples = 50)

shadealpha = 0.4
linesize = 2
qplot1 = ggplot() +
  geom_ribbon(aes(x = gevq$TrvalQ, ymax = gevq$qupper, ymin = gevq$qlower),
              alpha = shadealpha,  fill = "green") +
  geom_line(aes(gevq$TrvalQ, gevq$qupper), linetype = "dotted", 
            color = "green") +
  geom_line(aes(gevq$TrvalQ, gevq$qlower), linetype = "dotted", 
            color = "green") +
  geom_ribbon(aes(x = hmq$TrvalQ, ymax = hmq$qupper, ymin = hmq$qlower),
              alpha = shadealpha,  fill = "blue") +
  geom_line(aes(hmq$TrvalQ, hmq$qupper), linetype= "dotted", color = "blue") +
  geom_line(aes(hmq$TrvalQ, hmq$qlower), linetype = "dotted", color = "blue") +
  geom_ribbon(aes(x = shmq$TrvalQ, ymax = shmq$qupper, ymin = shmq$qlower),
              alpha = shadealpha,  fill = "red") +
  geom_line(aes(shmq$TrvalQ, shmq$qupper), linetype = "dotted", color = "red") + 
  geom_line(aes(shmq$TrvalQ, shmq$qlower), linetype = "dotted", color = "red") + 
  geom_line(aes(hmq$TrvalQ, hmq$qmean,color = "blue"), size = linesize) +
  geom_line(aes(shmq$TrvalQ, shmq$qmean,color = "red"), size = linesize) +
  geom_line(aes(gevq$TrvalQ, gevq$qmean,color = "green"), size = linesize) +
  geom_point(aes( dtsim$Tr, dtsim$Xi[i,]), color = 'black', shape = 1,
             size = 2.5, stroke = 1.2) +
  geom_line(aes(trtrue, qutrue),  linetype="solid", color = 'black', 
            size = 1.2*linesize) +
  coord_trans(x = "log10") +
  labs(y = "Quantile [mm/day]", x= "Return Time [years]") +
  scale_x_continuous(breaks=c(2, 4,7, 10, 20, 50, 100), limits = c(2, 110)) + 
  scale_color_identity(name = "Model fit",
                       breaks = c("red", "blue", "green"),
                       labels = c("sHMEV", "HMEV", "GEV"),
                       guide = "legend") +
  theme_bw() + 
  theme(legend.position = c(0.35, 0.8),
        legend.background = element_blank(),
        text = element_text(size = 14),
        legend.box.background = element_rect(colour = "black"))

qplot1

#-- station 2 example --#
i4 = 11

# HMEV 
dyn14 = hmevr::fit_ev_model(dtsim$dailyp[i4,,], model = 'wei_dgu_bin',  
                            iter = 2000, chains = 4, Mgen = 50)
hmq4 = hmevr::comp_quant(dyn14, dvsim$maxima[i4,], trmin = 2)

#GEV
gev14 = hmevr::fit_ev_model(dtsim$dailyp[i4,,], model = 'gev', iter = 2000, 
                            chains = 4, Mgen = 50)
gevq4 = hmevr::comp_quant(gev14, dvsim$maxima[i4,], trmin = 2)

# sHMEV
shmq4 = quant_spatial(spat.sim, maxval = dvsim$maxima[i4,], id_st = i4)

qutrue4 = qhbev1(trtrue, ptrue = ptrue, ntrue = lambda, id_st = i4,
                 Nt = 366, nsamples = 50)

qplot2 = ggplot() +
  geom_ribbon(aes(x = gevq4$TrvalQ, ymax = gevq4$qupper, ymin = gevq4$qlower),
              alpha = shadealpha,  fill = "green") +
  geom_line(aes(gevq4$TrvalQ, gevq4$qupper), linetype = "dotted",
            color = "green") +
  geom_line(aes(gevq4$TrvalQ, gevq4$qlower), linetype = "dotted",
            color = "green") +
  geom_ribbon(aes(x = hmq4$TrvalQ, ymax = hmq4$qupper, ymin = hmq4$qlower),
              alpha = shadealpha,  fill = "blue") +
  geom_line(aes(hmq4$TrvalQ, hmq4$qupper), linetype = "dotted", 
            color = "blue") +
  geom_line(aes(hmq4$TrvalQ, hmq4$qlower), linetype = "dotted", 
            color = "blue") +
  geom_ribbon(aes(x = shmq4$TrvalQ, ymax = shmq4$qupper, ymin = shmq4$qlower),
              alpha = shadealpha,  fill = "red") +
  geom_line(aes(shmq4$TrvalQ, shmq4$qupper), linetype = "dotted",
            color = "red") + 
  geom_line(aes(shmq4$TrvalQ, shmq4$qlower), linetype = "dotted", 
            color = "red") + 
  geom_line(aes(hmq4$TrvalQ, hmq4$qmean,color = "blue"), size = linesize) +
  geom_line(aes(shmq4$TrvalQ, shmq4$qmean,color = "red"), size = linesize) +
  geom_line(aes(gevq4$TrvalQ, gevq4$qmean,color = "green"), size = linesize) +
  geom_point(aes( dtsim$Tr, dtsim$Xi[i4,]), color = 'black', shape = 1,
             size = 2.5, stroke = 1.2) +
  geom_line(aes( trtrue, qutrue4),  linetype="solid", color = 'black', 
            size = 1.2*linesize) +
  coord_trans(x = "log10") +
  labs(y = "Quantile [mm/day]", x = "Return Time [years]") +
  scale_x_continuous(breaks=c(2, 4,7, 10, 20, 50, 100), limits = c(2, 110)) + 
  scale_color_identity(name = "Model fit",
                       breaks = c("red", "blue", "green"),
                       labels = c("sHMEV", "HMEV", "GEV"),
                       guide = "legend") +
  theme_bw() + 
  theme(legend.position = c(0.35, 0.8),
        legend.background = element_blank(),
        text = element_text(size = 14),
        legend.box.background = element_rect(colour = "black"))

qplot2

qplot = grid.arrange(qplot1, qplot2, ncol = 2)
# ggsave("quant_sim.png", plot = qplot, device = 'png', width = 8, height = 4.2)







