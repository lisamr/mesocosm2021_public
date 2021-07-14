# understanding community level changes in disease risk across the four treatments

rm(list = ls())
library(cowplot)
library(tidyverse)
library(rethinking)
library(loo)
library(bayesplot)

#===============================================================================
#functions
#===============================================================================

zscore <- function(x) (x - mean(x)) / (2*sd(x))
anti_zscore <- function(x, z) 2*sd(x)*z + mean(x)
#for prepping data list for stan
f_dat <- function(df, Xmat, X_phimat, Xsim, rich_sim){
  dat_prev <- list(
    N = nrow(Xmat),
    K = ncol(Xmat),
    Kphi = ncol(X_phimat),
    X = Xmat,
    day = as.integer(df$day),
    X_phi = X_phimat,
    I = as.integer(df$I_all),
    n = as.integer(df$n_allP.na),
    Nsim = nrow(Xsim),
    Xsim = Xsim,
    trayID = df$trayID,
    rich_sim = rich_sim
  )
  return(dat_prev)
}
f_dat2 <- function(data){
  
  #for simulating to new data
  richnessz_sim <- seq(min(data$richnessz), max(data$richnessz), length.out = 10)
  richness_sim <- anti_zscore(x = data$richness, richnessz_sim) 
  
  datlist <- list(
    #polynomial richness, phi doesn't vary by anything
    f_dat(df = data,
          Xmat = cbind(data$tempz, data$richnessz, I(data$richnessz^2)),
          X_phimat = cbind(data$richnessz),
          Xsim = cbind(0, richnessz_sim, richnessz_sim^2),
          rich_sim = richness_sim),
    #no polynomial richness, phi doesn't vary by anything
    f_dat(df = data,
          Xmat = cbind(data$tempz, data$richnessz),
          X_phimat = cbind(data$richnessz),
          Xsim = cbind(0, richnessz_sim),
          rich_sim = richness_sim)
  )
  return(datlist)
}

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5))
pal <- rev(wesanderson::wes_palette("Darjeeling1")[c(1,4,2,5)])



#===============================================================================
#load data
#===============================================================================

plantdf <- read.csv('output/plant_level_data.csv')
traydf <- read.csv('output/tray_level_data.csv')
trts <- readRDS('output/treatments_list.RDS')


#===============================================================================
# Statistical models
# test out linear and polynomial relationships between disease and richness, 
# run seperately for each treatment. allow variance terms to vary by richness.
#===============================================================================


#===============================================================================
# prep data. functions above can be used for both infection prevalence and density.
#===============================================================================
set.seed(123)
df_AD <- traydf %>% 
  filter(trayID %in% trts$add_det$trayID) %>% 
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM))
df_SD <- traydf %>% 
  filter(trayID %in% trts$sub_det$trayID) %>% 
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM))
df_AS <- traydf %>% 
  filter(trayID %in% trts$add_stoch$trayID) %>% 
  filter(!trayID %in% sample(1:10, 8)) %>% #need to remove 8 radish trays
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) 
df_SS <- traydf %>% 
  filter(trayID %in% trts$sub_stoch$trayID) %>% 
  filter(!trayID %in% sample(21:30, 9)) %>% #need to remove 9 radish trays
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) 
  
datlist_AD <- f_dat2(df_AD)
datlist_SD <- f_dat2(df_SD)
datlist_AS <- f_dat2(df_AS)
datlist_SS <- f_dat2(df_SS)


#===============================================================================
#prevalence--binomial likelihood
#===============================================================================

#stanBBphi <- stan_model('Stan_models/four_trts/betabinom_phivary.stan')
stanBB <- stan_model('Stan_models/four_trts/betabinom.stan')


runstan_p <- function(datlist){
  fit_p <- lapply(datlist[1:2], function(x){
    sampling(stanBB, data = x, 
             chains = 4, iter = 2000, cores = 4, control=list(adapt_delta=0.90))
  })
  precisoutput <- lapply(fit_p, function(x) 
    precis(x, pars = c('a0', 'beta', 'phi'), depth = 2, prob = .9))
  return(list(
    fits = fit_p,
    outputs = precisoutput
  ))
}

#fits
fits_pAD <- runstan_p(datlist_AD)
fits_pSD <- runstan_p(datlist_SD)
fits_pAS <- runstan_p(datlist_AS)
fits_pSS <- runstan_p(datlist_SS)

#lOOIC
L_pAD <- lapply(fits_pAD$fits, loo, moment_match = T)
L_pSD <- lapply(fits_pSD$fits, loo, moment_match = T)
L_pAS <- lapply(fits_pAS$fits, loo, moment_match = T)
L_pSS <- lapply(fits_pSS$fits, loo, moment_match = T)


#check out model fit. All look good.
post <- extract.samples(fits_pSD$fits[[1]])
ppc_dens_overlay(datlist_SD[[1]]$I, post$y_rep[1:50,])

#1 = polynomial richness, phi doesn't vary by anything
#2 = no polynomial richness, phi doesn't vary by anything

loo_compare(L_pAD)
fits_pAD$outputs
# elpd_diff se_diff
# model2  0.0       0.0   
# model1 -0.2       0.4
# [[2]]
# mean   sd    5%   95% n_eff Rhat4
# a0[1]    0.42 0.27 -0.03  0.88  4942     1
# a0[2]   -0.21 0.28 -0.67  0.24  5126     1
# a0[3]   -0.08 0.26 -0.50  0.34  5514     1
# beta[1] -0.55 0.32 -1.07 -0.03  5011     1
# beta[2]  0.64 0.30  0.15  1.13  5640     1
# phi      1.55 0.73  0.42  2.85  3546     1


loo_compare(L_pSD)
fits_pSD$outputs
# elpd_diff se_diff
# model1  0.0       0.0   
# model2 -8.5       1.4  
# [[1]]
# mean   sd    5%   95% n_eff Rhat4
# a0[1]    0.90 0.33  0.37  1.45  3035     1
# a0[2]    0.02 0.33 -0.49  0.57  2893     1
# a0[3]    0.32 0.31 -0.18  0.83  3446     1
# beta[1]  0.07 0.32 -0.46  0.58  4422     1
# beta[2] -1.64 0.32 -2.17 -1.11  4094     1
# beta[3]  2.23 0.69  1.07  3.37  2210     1
# phi      2.04 1.03  0.54  3.86  2664     1


loo_compare(L_pAS)
fits_pAS$outputs
# elpd_diff se_diff
# model1  0.0       0.0   
# model2 -0.3       1.1 
# [[2]]
# mean   sd    5%   95% n_eff Rhat4
# a0[1]   -0.63 0.22 -0.98 -0.27  6147     1
# a0[2]   -0.69 0.23 -1.08 -0.31  6378     1
# a0[3]   -0.97 0.24 -1.36 -0.57  6234     1
# beta[1] -0.18 0.27 -0.61  0.25  5485     1
# beta[2]  1.18 0.26  0.75  1.61  5243     1
# phi      4.23 1.40  2.17  6.70  4626     1


loo_compare(L_pSS)
fits_pSS$outputs
# elpd_diff se_diff
# model2  0.0       0.0   
# model1 -0.3       0.1   
#[[2]]
#          mean   sd    5%   95% n_eff Rhat4
# a0[1]   -1.16 0.33 -1.69 -0.62  3252     1
# a0[2]   -0.86 0.30 -1.35 -0.37  3365     1
# a0[3]   -1.20 0.29 -1.70 -0.73  3817     1
# beta[1] -0.18 0.34 -0.75  0.39  3624     1
# beta[2]  0.61 0.34  0.06  1.17  3413     1
# phi      0.69 0.47  0.08  1.54  3794     1


#===============================================================================
#density--neg. binom 
#===============================================================================

#wrote the models with brms. see 'four_trts_comlevel_brms.R'

#===============================================================================
# Plot predictions against observations
#===============================================================================

#function for summarizing posterior predictions
f_musim <- function(posterior, datlist){
  median = apply(posterior$mu_sim, 2, median)
  lower = apply(posterior$mu_sim, 2, HPDI, .9)[1,]
  upper = apply(posterior$mu_sim, 2, HPDI, .9)[2,]
  data.frame(richness = datlist$rich_sim, median, lower, upper)
}

#plotting prevalence
f_plot_p <- function(fit, datlist, Color, Title, Linetype = 1){
  tmp = f_musim(extract.samples(fit), datlist)
  traydf %>% 
    filter(trayID %in% datlist$trayID) %>% 
    ggplot(.) +
    geom_ribbon(data = tmp, aes(richness, ymin = lower, ymax = upper), alpha = .2) +
    geom_line(data = tmp, aes(x = richness, y = median),lty = Linetype, lwd = 1.5, color = grey(.6)) +
    geom_jitter(aes(richness, I_all/n_allP.na),
                height = 0, width = .05, size = 2, alpha = .7, color = Color) +
    labs(y = 'Community infection prevalence', title = Title) +
    scale_y_continuous(limits = c(0, 1))
}


#===============================================================================
# Community level Prevalence
#===============================================================================

plot_grid(
  f_plot_p(fits_pAD$fits[[2]], datlist_AD[[2]], pal[1], 'Additive, Deterministic')+ noaxislab,  
  f_plot_p(fits_pSD$fits[[2]], datlist_SD[[2]], pal[2], 'Substitutive, Deterministic')+ noaxislab, 
  f_plot_p(fits_pAS$fits[[2]], datlist_AS[[2]], pal[3], 'Additive, Stochastic')+ noaxislab, 
  f_plot_p(fits_pSS$fits[[2]], datlist_SS[[2]], pal[4], 'Substitutive, Stochastic', 2)+ noaxislab, 
  scale = .9, labels = "auto", hjust = -2
) +
  draw_label('Richness', x = .5, y = 0, vjust = -.5, size = 11) + 
  draw_label('Community disease prevalence', x = .01, y = .5, angle = 90, size = 11)
ggsave('figures/four_trts/comm_prev.pdf', device = 'pdf', units = 'mm', width = 181, height = 100)
