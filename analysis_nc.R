
rm(list = ls())

library(tidyverse)
library(gridExtra)
library(data.table)
library(elevatr)
library(raster)
library(maps) 
library(sp)
library(PBSmapping)
library(hmevr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(extraDistr)
library(bayesplot)
library(viridis)
library(latex2exp)
library(nleqslv)

load("datafit.RData")
source("sHMEV_function.R")


#------------------------------------------------------------------------------#

id_station_fit = unique(datafit$STATION)
id_station_test = unique(datatest$STATION)

# matrix with covariates 
covariate = datafit %>%  
  group_by(STATION) %>% 
  summarise(LATITUDE = max(LATITUDE),
            LONGITUDE = max(LONGITUDE),
            ELEVATION = max(ELEVATION),
            dist_coast = max(dist_coast)) %>%
  mutate(intercept = rep(1, length(id_station_fit)), .after = STATION) %>%
  dplyr::select(-STATION)

cov_scale = covariate = as.matrix(covariate) 
# standardize covariates
cov_scale[,c(2:5)] = scale(covariate[,2:5])

# analysis parameters
nyears = 20 # number of years for the training set
Mgen=50 # number of sample points for sHMEV latent parameters
Nt = 366 # number of events / block = number of days/year
S = nrow(cov_scale) # number of sites


#----------------------------# Define model data #-----------------------------#

# training set
d_train = data_nyears(datafit, nyears)

d_matrix = array(NA, dim=c(S,nyears,Nt)) # array with PRCP 
for(s in 1:S){
  dt = d_train %>% filter(STATION == id_station_fit[s])
  d_matrix[s,,] = table_max(dt, Nt = Nt)$data
}

# number of event with PRCP>0 for each station
N = matrix(NA, nrow = S, ncol = nyears)
for(s in 1:S){
  for(i in 1:nyears){
    samplei = d_matrix[s,i,]
    N[s,i] = length(samplei[samplei > 0])
  }
}

model_datas = list(M = nyears, Nt = Nt, y = d_matrix, N = N, S = S, 
                    X = cov_scale, Mgen = Mgen)

# hyperparameters 
mod_ippars = list(sg0prior = c(10, 0.35), # a, b of inverse gamma for gamma
                  sd0prior = c(10, 19.75), #  a, b of inverse gamma for delta
                  ma0prior = c(0.75, 0, 0, 0, 0), # mean of the normal for gamma
                  sa0prior = c(0.05, 0.05, 0.05, 0.05, 0.05),
                  # sd of the normal for gamma
                  mn0prior = c(9.8, 0, 0, 0, 0),# mean of the normal for delta
                  sn0prior = c(1.5, 0.5, 0.5, 0.5, 0.5), 
                  # sd of the normal for delta
                  mb0prior = c(-0.5, 2, -8, 10, 8), 
                  # mean of the normal for lambda
                  sb0prior = c(1, 3, 4, 2.5, 4)) # sd of the normal for lambda

model_data = append(model_datas, mod_ippars)

#----------------------------# Fit sHMEV model #-------------------------------#
sHMEV.nc = rstan::stan(file="sHMEV_nc.stan", data = model_data, 
                       iter = 2000, chains = 4, 
                       control = list(adapt_delta = 0.8,max_treedepth = 10))

#--- check convergence of the chains ---#

# sigma_gamma and sigma_delta 
parss = c("sg", "sdl")
stan_trace(sHMEV.nc, pars = parss) +
  theme(plot.title = element_text(hjust = 0.5)) 

# regression coefficients for lambda
parsb = c("bl[1]", "bl[2]", "bl[3]", "bl[4]", "bl[5]") 
stan_trace(sHMEV.nc, pars = parsb) +
  theme(plot.title = element_text(hjust = 0.5))

# regression coefficients for mu_delta
parsn = c("nd[1]", "nd[2]", "nd[3]", "nd[4]", "nd[5]") 
stan_trace(sHMEV.nc, pars = parsn) +
  theme(plot.title = element_text(hjust = 0.5)) 

# regression coefficients for mu_gamma
parsa = c("ag[1]", "ag[2]","ag[3]","ag[4]","ag[5]") 
stan_trace(sHMEV.nc, pars = parsa) +
  theme(plot.title = element_text(hjust = 0.5)) 


#--- mean, sd and credible intervals  ---#
param = c(parss, parsb, parsn, parsa)
d = as.data.frame(sHMEV.nc)
dd = d %>% dplyr::select(param)
ic = t(apply(dd, 2, quantile, probs=c(0.025, 0.975)))
colnames(ic) = c("lower", "upper")
meanp = apply(dd, 2, mean)
se = apply(dd, 2, sd)
cbind(meanp,se, ic)



#----------------------# Posterior predictive checks #-------------------------#

ppd5 = plot_ppd(model_fit = sHMEV.nc, id_st = sexample[2], data = datafit)
# ggsave('ppdc52.png', plot = ppd5, device = 'png', width = 12, height = 6.5)

ppd7 = plot_ppd(model_fit=sHMEV.nc, id_st = sexample[3], data = datafit)
#  ggsave('ppdc7f.png', plot = ppd7, device = 'png', width = 12, height = 6.5)


#-----------------------# Posterior predictive maps #--------------------------#

# create grid of points of North Carolina
g = gridPoints(name_state = "North Carolina")

# compute for each point the distance from the coast
nc_borders = read.table("nc_borders.txt")
colnames(nc_borders) = c("lat", "lon")

# borders on the coast
nc_borders = nc_borders %>% filter(lon > -79)

dist_coast = rep(NA, nrow(g))
for (i in 1: nrow(g)){
  dist = pointDistance(p1 = nc_borders[,2:1], p2 = g[i,1:2], lonlat = T)
  # distance (in meters) from the coast for each grid point
  dist_coast[i] = min(dist) /1000 # convert to km
}

g = cbind.data.frame(g, dist_coast)
statemap = map_data('state')
setnames(statemap, c('X', 'Y', 'PID', 'POS', 'region', 'subregion'))
statemap = statemap %>% filter(region == "north carolina")
statemap = clipPolys(statemap, xlim = c(-84.5, -75), ylim = c(33.5, 37))

# ppd of the parameters lambda, gamma and delta for each grid point  
ppd_mapnc = ppd_map(model_fit = sHMEV.nc, grid = g)

# plot the maps
vfill_ppd = c("gmean", "giqr", "dmean", "diqr", "mlambda", "iqrlambda")
ltitle_ppd = c(TeX("$\\gamma$  (mean)"), TeX("$\\gamma$  (IQR)"), 
           TeX("$\\delta$  (mean)"), TeX("$\\delta$  (IQR)"),
           TeX("$\\lambda$  (mean)"), TeX("$\\lambda$  (IQR)"))

plist_ppd = lapply(1:length(vfill_ppd), function(i) {
  p = ggplot(ppd_mapnc) +
    geom_tile(aes_string(x = "x", y = "y", fill = vfill_ppd[i])) +
    scale_fill_viridis(option = "inferno", name = "") +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle(ltitle_ppd[i]) +
    coord_map(xlim = c(-84,-75), ylim = c(33.9, 36.6)) +
    geom_polygon( data = statemap, aes(X, Y, group = PID), color = 'grey50',
      fill = NA) +
    theme_minimal() +
    theme(legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(size = 13))
  p
})

# plot for gamma 
gamma_map = grid.arrange(plist_ppd[[1]], plist_ppd[[2]], nrow = 1)
# ggsave('gamma_map.png', plot = gamma_map, device = 'png', width = 9,
# height = 8)

# plot for delta 
delta_map = grid.arrange(plist_ppd[[3]], plist_ppd[[4]], nrow = 1)
# ggsave('delta_map.png', plot = delta_map, device = 'png', width = 9, 
# height = 8)

# plot for lambda 
lambda_map = grid.arrange(plist_ppd[[5]], plist_ppd[[6]], nrow = 1)
# ggsave('lambda_map.png', plot = lambda_map, device = 'png', width = 9,
# height = 8)


#---------------------------# return level maps #------------------------------#

# return values for Tr=25 and Tr=50 years for each grid point 
qmap25 = quant_map(model_fit = sHMEV.nc, grid = g, Tr = 25)
qmap50 = quant_map(model_fit = sHMEV.nc, grid = g, Tr = 50)

data_qmap = cbind.data.frame(g, qm25 = qmap25$qmean, ql25 = qmap25$qlower,
  qu25 = qmap25$qupper, qm50 = qmap50$qmean, ql50 = qmap50$qlower,
  qu50 = qmap50$qupper)

# plot the maps
my_breaks = c(100,  130, 160, 200, 250,  300,  350)
vfill_q = c("qm25", "ql25", "qu25", "qm50", "ql50", "qu50")

plist_q = lapply(1:length(vfill_q), function(i) {
  p = ggplot(data_qmap) +
    geom_tile(aes_string(x = "x", y = "y", fill = vfill_q[i])) +
    scale_fill_viridis(option = "turbo", name = "", breaks = my_breaks,
                      labels = my_breaks, trans = "log", limits = c(100, 350)) +
    xlab("Longitude") + ylab("Latitude") +
    coord_map(xlim = c(-84,-75), ylim = c(33.9, 36.6)) +
    geom_polygon(data = statemap, aes(X, Y, group = PID), color = 'grey50',
                 fill = NA) +
    theme_minimal()
  
  if (i == 2 | i == 5)
    p = p + theme(legend.key.width  = unit(1.8, "cm"), 
                  plot.title = element_text(hjust = 0.5), 
                  legend.position = "top", text = element_text(size = 13),
                  legend.key.height = unit(0.55, 'cm'))
  else
    p = p + theme(legend.position = "none")
  
  return(p)
})

qmap_mean = grid.arrange(plist_q[[1]], plist_q[[4]], nrow = 1)
# ggsave('qmap_mean.png', plot = qmap_mean, device = 'png', width = 9, 
# height = 8)

qmap_low = grid.arrange(plist_q[[2]], plist_q[[5]], nrow = 1)
# ggsave('qmap_low1.png', plot = qmap_low, device = 'png', width = 9, 
# height = 8)

qmap_up = grid.arrange(plist_q[[3]], plist_q[[6]], nrow = 1)
# ggsave('qmap_up1.png', plot = qmap_up, device = 'png', width = 9, height = 8)



#--------------------------# spatial validation #------------------------------#

# covariates for the test set
cov_test = datatest %>%  
  group_by(STATION) %>% 
  summarise(LATITUDE = max(LATITUDE),
            LONGITUDE = max(LONGITUDE),
            ELEVATION = max(ELEVATION),
            dist_coast = max(dist_coast)) %>%
  mutate(intercept = rep(1, length(id_station_test)), .after=STATION) %>%
  dplyr::select(-STATION)

cov_test = covs_test = as.matrix(cov_test) 
# standardize covariates
covs_test[,c(2:5)] = scale(cov_test[,2:5])

# first station of the test set
d22 = datatest %>% filter(STATION == id_station_test[1])
dts22 = split_obs_data(d22, M_cal = 20, Nt = 366, cross_val = FALSE, 
                       reshuffle =  FALSE, flip_time = FALSE, 
                       reshuffle_days = FALSE, decluster = FALSE)
max_all22 = dts22$dataval$max # maxima of all series
res22 = table_max(d22, Nt = Nt)

# second station of the test set
d26 = datatest %>% filter(STATION == id_station_test[2])
dts26 = split_obs_data(d26, M_cal = 20, Nt = 366, cross_val = FALSE,
                       reshuffle =  FALSE, flip_time = FALSE, 
                       reshuffle_days = FALSE, decluster = FALSE)
max_all26 = dts26$dataval$max # maxima of all series
res26 = table_max(d26, Nt = Nt)


# compute the quantile on the test set
qval22 = quantile_validation(sHMEV.nc, maxval = max_all22, 
                             covariate = covs_test, id_st = 1)
qval26 = quantile_validation(sHMEV.nc, maxval = max_all26,
                             covariate = covs_test, id_st = 2) 


# quantile versus return time plots
shadealpha = 0.3
qv1 = ggplot() +
  geom_line(aes(qval22$TrvalQ, qval22$qupper), linetype = "dotted", 
            color = "red") +
  geom_line(aes(qval22$TrvalQ, qval22$qlower), linetype = "dotted",
            color = "red") +
  geom_ribbon(aes(x = qval22$TrvalQ, ymax = qval22$qupper, 
                  ymin = qval22$qlower), alpha = shadealpha, fill = "red") +
  geom_line(aes(qval22$TrvalQ, qval22$qmean, color = "red"), show.legend = F,
    size = 1) +
  geom_point(aes(res22$Tr, res22$Xi), size = 2.5, color = "black", shape = 1) +
  coord_trans(x = "log10") + # log scale
  ggtitle("Statesville") +
  labs(y = "Quantile [mm/day]", x = "Return Time [years]") +
  scale_x_continuous(breaks = c(2, 5, 10, 20, 50, 100, 150), 
                     limits = c(2, 131)) +
  theme_bw() +   
  theme(text = element_text(size = 14))

qv2 = ggplot() +
  geom_line(aes(qval26$TrvalQ, qval26$qupper), linetype = "dotted",
            color = "red") +
  geom_line(aes(qval26$TrvalQ, qval26$qlower), linetype = "dotted",
            color = "red") +
  geom_ribbon(aes(x = qval26$TrvalQ, ymax = qval26$qupper, 
                  ymin = qval26$qlower),
    alpha = shadealpha, fill = "red") +
  geom_line(aes(qval26$TrvalQ, qval26$qmean, color = "red"), show.legend = F,
    size = 1) +
  geom_point(aes(res26$Tr, res26$Xi), size = 2.5, color = "black", shape = 1) +
  coord_trans(x = "log10") +  # log scale
  ggtitle("Wilson") +
  labs(y = "Quantile [mm/day]", x = "Return Time [years]") +
  scale_x_continuous(breaks = c(2, 5, 10, 20, 50, 100, 150), 
                     limits = c(2, 131)) +
  theme_bw() + 
  theme(text = element_text(size = 14))

qv = grid.arrange(qv1, qv2, nrow = 1)
# ggsave('qvalidation.png', plot=qv, device = 'png', width = 8, height = 4.2)



#-------------------# comparison with competing methods #----------------------#

## compute fse, mbias, mwidth
nc_fse = nc_mbias = nc_mwidth = NULL

for(i in 1:S){
  
  d = datafit %>% filter(STATION == id_station_fit[i])
  res = table_max(d, Nt=366)
  dts = split_obs_data(d, M_cal = 20, Nt = 366, cross_val = FALSE, 
                       reshuffle =  FALSE, flip_time = FALSE, 
                       reshuffle_days = FALSE, decluster = FALSE)
  dat1   = dts$datacal$data # data for the estimation
  max_all = dts$dataval$max # maxima of all series
  
  print(paste("station", i))
  
  # HMEV 
  hmev = hmevr::fit_ev_model(dat1, model ='wei_dgu_bin', iter = 2000, 
                             chains = 4, Mgen = 50, adapt_delta = 0.8, 
                             thresh_hbev = 0)
  hmq = hmevr::comp_quant(hmev, max_all, trmin = 2)
  
  #GEV
  gev = hmevr::fit_ev_model(dat1, model = 'gev', iter =2000, chains=4, 
                            Mgen = 50, adapt_delta = 0.8)
  gevq = hmevr::comp_quant(gev, max_all, trmin = 2)
  
  # sHMEV 
  shmq = quant_spatial(sHMEV.nc, maxval = max_all, i=i)
  
  adt_f = cbind.data.frame(i,round(gevq$fse,3), round(hmq$fse,3),
                           round(shmq$fse,3))
  colnames(adt_f) = c("station", "gev", "hmev", "shmev")
  nc_fse = rbind(nc_fse,adt_f)
  
  adt_w = cbind.data.frame(i,round(gevq$mwidth,3), 
                           round(hmq$mwidth,3), round(shmq$mwidth,3))
  colnames(adt_w) = c("station", "gev", "hmev", "shmev")
  nc_mwidth = rbind(nc_mwidth, adt_w)
  
  adt_b = cbind.data.frame(i,round(gevq$mbias,3), 
                           round(hmq$mbias,3), round(shmq$mbias,3))
  colnames(adt_b) = c("station", "gev", "hmev", "shmev")
  nc_mbias = rbind(nc_mbias, adt_b)
}     


colnames(nc_fse)[2:4] = colnames(nc_mbias)[2:4] = colnames(nc_mwidth)[2:4] = 
  c("GEV", "HMEV", "sHMEV")

fse_long = gather(nc_fse, "model", "FSE", 2:4)
fse_long$model = as.factor(fse_long$model)
mbias_long = gather(nc_mbias, "model", "mbias", 2:4)
mwidth_long = gather(nc_mwidth, "model", "mwidth", 2:4)
data_fse = cbind(fse_long, mbias_long$mbias, mwidth_long$mwidth)
colnames(data_fse)[4:5] = c("mbias", "mwidth")

#plot the results
yvar = c("FSE", "mbias", "mwidth")

plist_fse = lapply(1:length(yvar), function(i) {
  p = ggplot(data_fse, aes_string(x = "model", y = yvar[i], fill = "model"))+
    geom_boxplot(alpha = 0.6) +
    scale_fill_manual(values = c("green", "blue", "red")) +
    theme_bw()  +
    theme( legend.position = 'none',
           text = element_text(size = 15),
           strip.background = element_rect(fill = "white")) 
  return(p)
})

fse_nc = grid.arrange(grobs = plist_fse, nrow = 1)
# ggsave('fse_nc.png', plot = fse_nc ,device = 'png', width = 8, height = 5)










