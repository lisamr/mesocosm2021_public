# understanding individual-level disease risk in the four different treatments
# models run seperately for each treatment


rm(list = ls())
library(cowplot)
library(tidyverse)
library(rethinking)
library(bayesplot)
library(loo)
zscore <- function(x) (x - mean(x)) / (2*sd(x))
anti_zscore <- function(x, z) 2*sd(x)*z + mean(x)

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5), 
                   legend.position = 'none')
sp_pal <- c('#F16745', '#FFC65D', '#7BC8A4', '#4CC3D9', '#93648D', '#404040')

#===============================================================================
#load data
#===============================================================================

plantdf <- read_csv('output/plant_level_data.csv')
traydf <- read_csv('output/tray_level_data.csv')
trts <- readRDS('output/treatments_list.RDS')

#For sub/stoch: add in the single species trays for species other than radish too. It would only help.
extras <- traydf %>% filter(trayID >=200 & trayID < 300) %>% pull(trayID)
trts$sub_stoch <- trts$sub_stoch %>% add_case(trayID = extras)

get_dataframes <- function(treatment, spIDs = 1:6){
  
  plantdata <- plantdf %>% 
    filter(trayID %in% treatment$trayID,
           trayID != 132) %>% 
    filter(spID %in% spIDs) %>% 
    mutate(trayID_consec = as.integer(as.factor(trayID)),
           spID = as.integer(spID)) %>% 
    group_by(trayID,trayID_consec, species, spID) %>% 
    summarise(I = sum(state_final == 'I', na.rm = T ),
              n = sum(!is.na(state_final), na.rm = T )) %>% 
    ungroup() %>% 
    arrange(trayID_consec)
  
  traydata <- traydf %>% 
    filter(trayID %in% treatment$trayID,
           trayID != 132) %>% 
    mutate(trayID_consec = as.integer(as.factor(trayID)),
           day = as.integer(day)) %>% 
    mutate_at(vars(temp_PM, richness, n_radish, n_allP.na), list(z=zscore)) %>% 
    arrange(trayID_consec)
  
  return(list(plantdata = plantdata, traydata = traydata))
}

d_AD <- get_dataframes(trts$add_det, 1:4)
d_AS <- get_dataframes(trts$add_stoch)
d_SD <- get_dataframes(trts$sub_det, 1:4)
d_SS <- get_dataframes(trts$sub_stoch)



#===============================================================================
#prep data for stan
#===============================================================================

f_datlist <- function(dataframes, variables, spIDs = 1:6, poly = F){
  traydata <- dataframes$traydata
  plantdata <- dataframes$plantdata
  traydatamat <- as.matrix(select(traydata, variables))
  densmat <- as.matrix(select(traydata, n_radish:n_butter))[,spIDs]
  
  #new data to simulate predictions
  minr <- min(traydata$richness_z)
  maxr <- max(traydata$richness_z)
  newdat <- expand_grid(
    spID = spIDs,
    intercept = 1,
    richness_z = seq(minr, maxr, length.out = 25)
  ) %>% 
    mutate(richness_z2 = richness_z^2)
  rich_sim <- anti_zscore(x = traydata$richness, newdat$richness_z) 
  richness <- traydata$richness_z[plantdata$trayID_consec]
  
  #form data list
  datlist <- list(
    N = nrow(plantdata),
    S = max(plantdata$spID),
    M = max(plantdata$trayID_consec),
    K = ncol(traydatamat),
    X = traydatamat,
    y = plantdata$I,
    n = plantdata$n,
    spID = plantdata$spID,
    mID = plantdata$trayID_consec,
    Day = traydata$day[plantdata$trayID_consec],
    XS = cbind(1, richness),
    
    #things for predicting to new data
    N_sim = nrow(newdat),
    spID_sim = newdat$spID,
    XS_sim = as.matrix(newdat[,2:3]),
    rich_sim = rich_sim,#not used in model, but useful for plotting
    
    ##things for the CC models (not used, but keeping anyways)
    Ks = ncol(densmat),
    X_dens = densmat
  )
  
  if(poly == T){
    datlist$XS <- cbind(1, richness, richness^2)
    datlist$XS_sim <-  as.matrix(newdat[,2:4])
  }
  
  return(datlist)
}



#===============================================================================
# run models
#===============================================================================

#stan model
stan_X <- stan_model('Stan_models/no_temporal_component/specieslevel_betabinom_notime_X.stan')
#stan_X <- stan_model('Stan_models/no_temporal_component/specieslevel_binom_notime_X.stan')
#stan_X2 <- stan_model('Stan_models/no_temporal_component/specieslevel_binom_notime_X2.stan')
# DISEASE RISK ~ A_SP + BETA_SP*RICHNESS + A_TRAY, controlling for temperature and day planted


#data lists
datlist_SD <- f_datlist(d_SD, c('temp_PM_z'), 1:4)
#datlist_SD2 <- f_datlist(d_SD2, c('temp_PM_z'), 1:4, poly = T)
datlist_SS <- f_datlist(d_SS, c('temp_PM_z')) 
datlist_AD <- f_datlist(d_AD, c('temp_PM_z'), 1:4) 
datlist_AS <- f_datlist(d_AS, c('temp_PM_z')) 

#fit models
fit_SD <- sampling(stan_X, data = datlist_SD, chains = 4, iter = 2000, cores = 4)
#fit_SD2 <- sampling(stan_X2, data = datlist_SD2, chains = 1, iter = 2000, cores =4 ) #includes polynomial term
fit_SS <- sampling(stan_X, data = datlist_SS, chains = 4, iter = 2000, cores = 4)
fit_AD <- sampling(stan_X, data = datlist_AD, chains = 4, iter = 2000, cores = 4)
fit_AS <- sampling(stan_X, data = datlist_AS, chains = 4, iter = 2000, cores = 4)

#check out posteriors
precis(fit_SD, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
precis(fit_SS, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
precis(fit_AD, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
precis(fit_AS, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)




#===============================================================================
#check out model fit
#===============================================================================
#looks perfect for all!
post <- extract(fit_SS)
ppc_dens_overlay(datlist_SS$y, post$y_rep[1:50,]) 
ppc_scatter_avg(datlist_SS$y, post$y_rep)

#looks good :)
pairs(fit_SD, pars = c('a0[1]', 'a0[2]', 'betaS', 'Beta', 'zS[2,1]')) 
pairs(fit_SD, pars = c('zS[2,1]', 'zS[2,2]')) 


#===============================================================================
#Plot predictions against observed data
#===============================================================================

#function for summarizing posterior predictions
f_psim <- function(posterior, datalist){
  #posterior = extract(fit_SS)
  #datalist = datlist_SS
  psim = posterior$p_sim
  median = apply(psim, 2, median)
  lower = apply(psim, 2, HPDI, .9)[1,]
  upper = apply(psim, 2, HPDI, .9)[2,]
  df <- data.frame(spID = datalist$spID_sim,
             richness = datalist$rich_sim, 
             median, lower, upper)
  
  #does the richness effect cross zero?
  tmp <- apply(posterior$sp_slope1, 2, HPDI, .9) %>% 
    apply(., 2, prod) 
  df <- df %>% 
    left_join(data.frame(spID = 1:length(tmp), crosseszero = tmp < 0)) %>% 
    mutate(lty = ifelse(crosseszero == T, 2, 1) )

  # remove predictions where data doesn't exist
  if(length(tmp) == 4){
    df <- df %>% 
      filter(!(spID >= 3 & richness < 4),
             !(spID == 2 & richness < 2))
  }
  return(df)
  
}

#plotting prevalence
f_plot <- function(fit, data, datalist, ribbon=F, Title){
  #fit = fit_SS
  #data = d_SS
  #datalist = datlist_SS
  psim <- f_psim(extract(fit), datalist)
  spIDs <- 1:datalist$S
  d <- data$plantdata %>%
    left_join(data$traydata) %>% 
    mutate(richnessplus = richness - .2 + spID*.1) %>% 
    filter(spID %in% spIDs)
  
  p <- ggplot(d, aes(group = as.factor(spID))) +
    geom_line(data = psim, aes(richness, median, color = as.factor(spID)), linetype = psim$lty, lwd = 1) +
    geom_point( aes(richnessplus, I/n, color = as.factor(spID)), size = 2, alpha = 1) +
    scale_color_manual(values = sp_pal[spIDs]) +
    scale_fill_manual(values = sp_pal[spIDs]) +
    ggtitle(Title)
  p
  if(ribbon == T){
    p <- p+  
      geom_ribbon(data = psim, alpha = .4,
                  aes(richness, ymin = lower, ymax = upper, fill = as.factor(spID)))
  }
    
  return(p)
}

#===============================================================================
# plot species level Prevalence
#===============================================================================

f_plot(fit_SD, d_SD, datlist_SD, T, 'Substitutive, Deterministic') + facet_wrap(~as.factor(spID))
f_plot(fit_SS, d_SS, datlist_SS, T, 'Substitutive, Stochastic') + facet_wrap(~as.factor(spID))
f_plot(fit_AS, d_AS, datlist_AS, T, 'Additive, Stochastic') + facet_wrap(~as.factor(spID))
f_plot(fit_AD, d_AD, datlist_AD, T, 'Additive, Deterministic') + facet_wrap(~as.factor(spID))

grid1 <- plot_grid(
  f_plot(fit_AD, d_AD, datlist_AD, F, 'Additive, Deterministic') + noaxislab,  
  f_plot(fit_SD, d_SD, datlist_SD, F, 'Substitutive, Deterministic')+ noaxislab,  
  f_plot(fit_AS, d_AS, datlist_AS, F, 'Additive, Stochastic') + noaxislab,  
  f_plot(fit_SS, d_SS, datlist_SS, F, 'Substitutive, Stochastic') + noaxislab,  
  scale = .9, labels = "auto", hjust = -2
)+
  draw_label('Richness', x = .5, y = 0, vjust = -.5, size = 11) + 
  draw_label('Species-level disease prevalence', x = .01, y = .5, angle = 90, size = 11)

#make legend for species
sp <- c('Radish', 'Arugula', 'Basil', 'Red romaine', 
        'Green lettuce', 'Butter lettuce')
plegend <- data.frame(spID = 1:6, 
                      Species = factor(sp, levels = sp)) %>% 
  ggplot(., aes(spID, Species, color = Species)) +
  geom_point(size = 2) +
  scale_color_manual(values = sp_pal) +
  theme(legend.position = 'bottom')

final_plot <- plot_grid(grid1, get_legend(plegend), nrow = 2, rel_heights = c(.9, .1), scale = c(.95, .85))
final_plot

pdf('figures/four_trts/species_DDR.pdf', width = 7.12598, height = 4)
final_plot
dev.off()
#save_plot(filename = 'figures/four_trts/species_DDR', plot = final_plot, device = 'pdf', units = 'mm', base_width = 181, base_height = 100)


#===============================================================================
# Export model output
#===============================================================================

precis(fit_SD, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'z', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
precis(fit_SS, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
precis(fit_AD, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
precis(fit_AS, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
# 
# > precis(fit_SD, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'z', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
# mean   sd    5%   95% n_eff Rhat4
# a0[1]          1.16 0.57  0.20  2.06  1043     1
# a0[2]          0.12 0.57 -0.88  1.02  1180     1
# a0[3]          0.50 0.57 -0.46  1.44  1116     1
# sp_int[1]      1.47 0.24  1.09  1.86  1148     1
# sp_int[2]      1.17 0.25  0.77  1.58  1476     1
# sp_int[3]      0.00 0.33 -0.51  0.58  1884     1
# sp_int[4]      0.48 0.31 -0.03  1.00  1739     1
# sp_slope1[1]  -0.96 0.44 -1.69 -0.23  1521     1
# sp_slope1[2]  -1.02 0.47 -1.82 -0.27  1619     1
# sp_slope1[3]  -1.33 0.63 -2.45 -0.37  1746     1
# sp_slope1[4]  -0.99 0.55 -1.90 -0.10  1844     1
# betaS         -1.00 0.50 -1.82 -0.17  1499     1
# z[1,1]         0.88 0.47  0.20  1.74  1129     1
# z[1,2]         0.57 0.47 -0.11  1.43  1181     1
# z[1,3]        -0.59 0.47 -1.36  0.21  1204     1
# z[1,4]        -0.11 0.47 -0.84  0.72  1175     1
# z[2,1]         0.04 0.34 -0.52  0.57  2116     1
# z[2,2]        -0.03 0.35 -0.63  0.48  2018     1
# z[2,3]        -0.33 0.49 -1.33  0.18  1831     1
# z[2,4]         0.01 0.37 -0.60  0.58  2689     1
# Beta[1]        0.00 0.41 -0.67  0.68  1430     1
# corr_mat[1,2]  0.10 0.44 -0.65  0.79  3555     1
# > precis(fit_SS, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
# mean   sd    5%   95% n_eff Rhat4
# a0[1]         -0.59 0.56 -1.49  0.38  1286     1
# a0[2]         -0.62 0.55 -1.50  0.31  1260     1
# a0[3]         -1.33 0.58 -2.24 -0.36   956     1
# sp_int[1]      1.17 0.36  0.57  1.75  1368     1
# sp_int[2]     -0.56 0.28 -1.02 -0.10  1894     1
# sp_int[3]     -2.17 0.35 -2.75 -1.61  2405     1
# sp_int[4]     -1.63 0.32 -2.15 -1.12  1723     1
# sp_int[5]     -1.87 0.33 -2.42 -1.33  1497     1
# sp_int[6]     -2.06 0.30 -2.53 -1.56  2000     1
# sp_slope1[1]  -1.63 0.51 -2.44 -0.76   945     1
# sp_slope1[2]   0.16 0.44 -0.56  0.88  1656     1
# sp_slope1[3]   0.05 0.50 -0.77  0.86  1877     1
# sp_slope1[4]   0.52 0.46 -0.26  1.27  1909     1
# sp_slope1[5]   0.61 0.48 -0.18  1.41  1562     1
# sp_slope1[6]   0.53 0.45 -0.21  1.27  1854     1
# betaS         -0.07 0.47 -0.83  0.71  1120     1
# Beta[1]       -0.43 0.37 -1.03  0.14  1838     1
# corr_mat[1,2] -0.54 0.30 -0.90  0.04  2653     1
# > precis(fit_AD, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
# mean   sd    5%   95% n_eff Rhat4
# a0[1]          0.34 0.51 -0.49  1.17   842     1
# a0[2]         -0.27 0.50 -1.08  0.57  1027     1
# a0[3]         -0.36 0.51 -1.21  0.46   769     1
# sp_int[1]      0.50 0.20  0.18  0.84  1226     1
# sp_int[2]      0.27 0.22 -0.08  0.62  1361     1
# sp_int[3]     -0.72 0.33 -1.24 -0.15  1677     1
# sp_int[4]     -0.38 0.32 -0.88  0.18  1601     1
# sp_slope1[1]   1.63 0.41  0.96  2.33  1350     1
# sp_slope1[2]   1.48 0.43  0.77  2.20  1580     1
# sp_slope1[3]   0.96 0.62 -0.08  1.92  1696     1
# sp_slope1[4]   1.02 0.59 -0.01  1.94  1759     1
# betaS          1.14 0.55  0.16  1.98  1108     1
# Beta[1]       -0.60 0.37 -1.22 -0.02  1277     1
# corr_mat[1,2]  0.23 0.42 -0.51  0.84  2327     1
# > precis(fit_AS, pars = c('a0', 'sp_int', 'sp_slope1', 'betaS', 'Beta', 'corr_mat[1,2]'), depth = 3, prob = .9)
# mean   sd    5%   95% n_eff Rhat4
# a0[1]         -0.68 0.37 -1.27 -0.03   585  1.01
# a0[2]         -0.72 0.38 -1.32 -0.09   738  1.01
# a0[3]         -1.25 0.38 -1.85 -0.60   735  1.00
# sp_int[1]     -0.07 0.16 -0.33  0.19  1752  1.00
# sp_int[2]     -0.68 0.20 -1.02 -0.36  1698  1.00
# sp_int[3]     -1.46 0.23 -1.83 -1.08  1792  1.00
# sp_int[4]     -1.06 0.20 -1.39 -0.73  2204  1.00
# sp_int[5]     -1.30 0.23 -1.68 -0.94  1885  1.00
# sp_int[6]     -1.39 0.23 -1.78 -1.01  2599  1.00
# sp_slope1[1]   1.32 0.29  0.85  1.81  1613  1.00
# sp_slope1[2]   1.56 0.39  0.97  2.23  1403  1.00
# sp_slope1[3]   0.97 0.39  0.31  1.57  1543  1.00
# sp_slope1[4]   1.17 0.33  0.60  1.71  1628  1.00
# sp_slope1[5]   1.42 0.37  0.85  2.05  1569  1.00
# sp_slope1[6]   1.20 0.35  0.61  1.77  1769  1.00
# betaS          1.23 0.33  0.71  1.74  1314  1.00
# Beta[1]       -0.48 0.27 -0.92 -0.03  1298  1.00
# corr_mat[1,2]  0.09 0.41 -0.61  0.74  2662  1.00