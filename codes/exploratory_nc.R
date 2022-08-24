
rm(list=ls()) 
library(tidyverse)
library(devtools)
library(gridExtra)
library(data.table)
library(latex2exp)
library(elevatr)
library(raster)
library(maps) 
library(sp)
library(PBSmapping)
library(hmevr)

source("sHMEV_function.R")
load("nc.RData")


#------------------------------------------------------------------------------#

# splitting data:  training, test and station for illustrative purpose data
id_station = unique(data$STATION)
n_station = as.integer(length(id_station))
stest = sample(1:n_station,2)
sexample = c(3,5,7)

datafit = data %>% filter(STATION %in% id_station[-stest])
datatest = setdiff(data, datafit)
# save(datafit, datatest, stest, sexample, file= "datafit.RData")


#------------------------------------------------------------------------------#
# Map of North Carolina

data_station = data %>% 
  group_by(STATION, LATITUDE, LONGITUDE) %>% 
  summarise(id = 1) 
data_station$id = c(1:nrow(data_station))

station_fit = data_station %>% filter(STATION %in% id_station[-stest])
station_test = setdiff(data_station, station_fit)
station_fit1 = station_fit %>% filter(STATION %in% id_station[-sexample])
station_fit2 = setdiff(station_fit, station_fit1)

# create grid of points of North Carolina
g = gridPoints(name_state = "North Carolina")

statemap = map_data('state')
setnames(statemap, c('X', 'Y', 'PID', 'POS', 'region', 'subregion'))
statemap = statemap %>% filter(region == "north carolina")
statemap = clipPolys(statemap, xlim = c(-84.5, -75), ylim = c(33.5, 37))

map_nc = ggplot(g) + 
  geom_tile(aes(x = x, y = y, fill = elevation)) + 
  scale_fill_gradientn(colours = terrain.colors(100)) +
  xlab("Longitude") + ylab("Latitude") +
  coord_map(xlim = c(-84, -75), ylim = c(33.9, 36.6))+
  geom_polygon(data = statemap, aes(X ,Y, group = PID), 
               color = 'grey50', fill = NA)+
  geom_point(data = station_fit1, aes(x = LONGITUDE, y = LATITUDE), size = 3,
             shape = 21, fill = "grey17") +
  geom_point(data = station_fit2, aes(x = LONGITUDE, y = LATITUDE), shape = 24,
             fill = "grey17", size = 3) +
  geom_point(data = station_test, aes(x = LONGITUDE, y = LATITUDE), 
             col = "blue", size = 3) +
  theme_minimal() 

map_nc
# ggsave('map_nc.png', plot=map_nc,device = 'png', width = 5, height = 2)


#------------------------------------------------------------------------------#
#--------------------------# Exploratory analysis #----------------------------#

# stations taken as an example
data_ex = data %>% filter(STATION %in% id_station[sexample])
id_ex = unique(data_ex$STATION)

#--- ACF ---#
par(mfrow=c(1,3))
data_ex %>% filter(STATION == id_ex[1]) %>%
  dplyr::select(PRCP) %>% acf(lag=25, main="Edenton")
data_ex %>% filter(STATION == id_ex[2]) %>%
  dplyr::select(PRCP) %>% acf(lag=25, main="Fayetteville")
data_ex %>% filter(STATION == id_ex[3]) %>%
  dplyr::select(PRCP) %>% acf(lag=25, main="Hendersonville")
par(mfrow=c(1,1))


#--- annual maxima ---#
data_year = data %>%
  group_by(STATION, LATITUDE, LONGITUDE, dist_coast, YEAR) %>%
  summarise(max_year = max(PRCP),
            mean_year = mean(PRCP),
            mean0_year = mean(PRCP[PRCP>0]),
            sd_year = sd(PRCP),
            iqr_year = IQR(PRCP),
            dayp = length(PRCP[PRCP>0]))

spat_avg = data_year %>% group_by(STATION, LATITUDE, LONGITUDE, dist_coast) %>% 
  summarise(media_max = mean(max_year),
            sd_max = sd(max_year),
            media_PRCP = mean(mean_year), 
            media_PRCP0 = mean(mean0_year), 
            media_sd = mean(sd_year), 
            media_day = mean(dayp))
el.nc = data %>% group_by(STATION) %>%
  dplyr::summarise(ELEVATION = max(ELEVATION))
spat_avg = cbind.data.frame(spat_avg, el.nc$ELEVATION)
colnames(spat_avg)[11] = "ELEVATION"


plot_max = ggplot(spat_avg) +
  geom_point(aes(x = LONGITUDE, y = LATITUDE, size = sd_max, 
                 colour = media_max), show.legend=F, alpha=1) +
  geom_point(aes(x = LONGITUDE, y = LATITUDE, size = sd_max, 
                 colour = media_max), show.legend=T, alpha=1) +
  labs(color = "mean maxima") +
  scale_colour_gradient() +
  xlab("Longitude") + ylab("Latitude") +
  coord_map(xlim = c(-84, -75), ylim = c(33.9, 36.6))+
  geom_polygon(data = statemap, aes(X, Y, group = PID), color = 'grey50', 
               fill = NA) +
  theme_minimal() + 
  theme(text = element_text(size = 14))

plot_max
# ggsave('annual_maxima.png', plot=plot_max, device = 'png', 
# width = 9, height = 8)


#--- scatter plot ---#

gamma.hat = beta.hat = rep(0, n_station)
for(i in 1:n_station){
  data_we = data %>% filter(STATION == id_station[i] & PRCP >0 )
  m1_emp = mean(data_we$PRCP)
  m2_emp = 1/nrow(data_we)*sum(data_we$PRCP^2)
  var_emp = var(data_we$PRCP)
  f = function(x) gamma(1 + 2/x)/gamma(1 + 1/x)^2 - m2_emp/m1_emp^2
  st = (m1_emp/sqrt(var_emp))^(1.086)
  res = uniroot(f, c(st, 1000), lowe = 0.1)
  gamma.hat[i] = res$root
  beta.hat[i] = m1_emp/gamma(1 + 1/gamma.hat[i])
}
spat_avg = cbind.data.frame(spat_avg, gamma.hat, beta.hat)

xlabel = c("Latitude", "Longitude", "Distance from coast [km]", "Elevation")
vlist = c("LATITUDE", "LONGITUDE", "dist_coast", "ELEVATION")

# plot for delta
plot_delta = lapply(1:length(xlabel), function(i) {
  p = ggplot(spat_avg) +
    geom_point(aes_string(vlist[i],  "beta.hat")) +
    xlab(xlabel[i]) + ylab(TeX("$\\tilde{\\delta}$")) +
    geom_smooth(aes_string(vlist[i],  "beta.hat"), method = lm, se = FALSE) + 
    theme_bw() + theme(text = element_text(size = 15))
})
delta_scatter = grid.arrange(grobs = plot_delta, ncol = 2)

#plot for gamma
plot_gamma = lapply(1:length(xlabel), function(i) {
  p = ggplot(spat_avg) +
    geom_point(aes_string(vlist[i],  "gamma.hat")) +
    xlab(xlabel[i]) + ylab(TeX("$\\tilde{\\gamma}$")) +
    geom_smooth(aes_string(vlist[i],  "gamma.hat"), method = lm, se = FALSE) + 
    theme_bw() + theme(text = element_text(size = 15))
})
gamma_scatter = grid.arrange(grobs = plot_gamma, ncol = 2)


# plot for n_j
plot_n = lapply(1:length(xlabel), function(i) {
  p = ggplot(spat_avg) +
    geom_point(aes_string(vlist[i],  "media_day")) +
    xlab(xlabel[i]) + ylab("Annual day of rain (mean)")+ 
    geom_smooth(aes_string(vlist[i],  "media_day"), method = lm, se = FALSE) + 
    theme_bw() + theme(text = element_text(size = 15))
})
n_scatter = grid.arrange(grobs = plot_n, ncol = 2)
