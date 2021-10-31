# understanding species-level disease risk using data from all trays


rm(list = ls())
library(cowplot)
library(tidyverse)
library(brms)

zscore <- function(x) (x - mean(x)) / (2*sd(x))
anti_zscore <- function(x, z) 2*sd(x)*z + mean(x)
sumtozero <- function(df){
  contrasts(df$day) <- contr.sum( length(unique( df$day )))
  return(df)
}
source('Analysis_scripts/betabinom_custom_brms.R') #for using beta-binomial likelihood with brms

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'none')

#===============================================================================
#load data
#===============================================================================

plantdf <- read_csv('Data/plant_level_data.csv')
traydf <- read_csv('Data/tray_level_data.csv')


#merge data sets
plantdata <- plantdf %>% 
  filter(trayID < 300) %>% #remove uninoculated trays
  filter(trayID != 132, tray_type != 'augmented') %>%  #and remove the unwatered tray
  group_by(trayID, species, spID) %>% 
  summarise(I = sum(state0_v2=='S' & state_final == 'I', na.rm = T ),
            n = sum(state0_v2=='S' & !is.na(state_final), na.rm = T )) %>% 
  ungroup() 

traydata <- traydf %>% 
  filter(trayID %in% unique(plantdata$trayID)) %>% 
  mutate(n_others = n_allP.na - n_radish, 
         n_others2 = n_allP.na - n_radish - n_arugula) %>% 
  mutate_at(vars(temp_PM, richness, n_others, n_others2, n_allP.na), list(z=zscore)) %>% 
  mutate_at(vars(n_radish:n_butter), list(zsq = function(x) zscore(sqrt(x)))) %>% 
  mutate_at(vars(n_radish:n_butter), list(p = function(x) x/.$n_allP.na))

merged <- left_join(plantdata, traydata) %>% 
  mutate_at(vars(trayID, day, spID), as.factor) 

merged <- sumtozero(merged)



#===============================================================================
#run brms model
#===============================================================================
# model   tray cov              
# 0:      --                     
# 1:      nall                   
# 2:      nradish
# 3:      nothers, nradish
# 4:      nothers2, nradish, narugula
# 5:      nradish, narugula...nbutter
# 6:      nradish, narugula

bf0 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + (1|trayID) + (1 + richness_z | spID))
bf1 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_allP.na_z + (1|trayID) + (1 + richness_z | spID))
bf2 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + (1|trayID) + (1 + richness_z | spID))
bf3 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_others_z + n_radish_zsq + (1|trayID) + (1 + richness_z | spID))
bf4 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_others2_z + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))
bf5 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + n_arugula_zsq + n_basil_zsq + n_red_rom_zsq + n_green_rom_zsq + n_butter_zsq + (1|trayID) + (1 + richness_z | spID))
#bf6 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))

fit0 <- brm(bf0, data = merged,
            prior = c(prior(normal(0, 1.5), class = Intercept),
                      prior(normal(0, 1), class = b),
                      prior(lkj(2), class = cor),
                      prior(normal(0,2), class = sd), 
                      prior(exponential(1), class = phi)),
            family = beta_binomial2, stanvars = stanvars,
            iter = 2000, chains = 4, cores = 4, 
            control = list(adapt_delta = 0.95, max_treedepth = 15), 
            save_pars = save_pars(all = TRUE))
fit1 <- update(fit0, formula. = bf1, newdata = merged, save_pars = save_pars(all = TRUE))
fit2 <- update(fit0, formula. = bf2, newdata = merged, save_pars = save_pars(all = TRUE))
fit3 <- update(fit0, formula. = bf3, newdata = merged, save_pars = save_pars(all = TRUE))
fit4 <- update(fit0, formula. = bf4, newdata = merged, save_pars = save_pars(all = TRUE))
fit5 <- update(fit0, formula. = bf5, newdata = merged, save_pars = save_pars(all = TRUE))


#===============================================================================
#post-process--checking and comparing models
#===============================================================================
#stan functions need to be exposed for each fit. produces lots of messages, they're ok, ignore them.
expose_functions(fit4, vectorize = T)

pp_check(fit4, ndraws = 50)

K0 <- kfold(fit0, K = 10)
K1 <- kfold(fit1, K = 10)
K2 <- kfold(fit2, K = 10)
K3 <- kfold(fit3, K = 10)
K4 <- kfold(fit4, K = 10)
K5 <- kfold(fit5, K = 10)
loo_compare(K0, K1, K2, K3, K4, K5)
# model   tray cov              
# 0:      --                     
# 1:      nall                   
# 2:      nradish
# 3:      nothers, nradish
# 4:      nothers2, nradish, narugula
# 5:      nradish, narugula...nbutter
# 6:      nradish, narugula

# elpd_diff se_diff
# fit4   0.0       0.0  
# fit3  -4.8       3.6  
# fit5  -5.5       4.1  
# fit2  -9.8       5.7  
# fit0 -43.8       7.5  
# fit1 -53.2       8.2  

loo_compare(K1, K3, K4, K5)
# elpd_diff se_diff
# fit4   0.0       0.0  
# fit3  -4.8       3.6  
# fit5  -5.5       4.1  
# fit1 -53.2       8.2  
#Export sutff for figure making-------------------------------------------------
write_rds(fit1, 'Outputs/splevel_alltrays_fit1_noaug.rds')
write_rds(fit3, 'Outputs/splevel_alltrays_fit3_noaug.rds')
write_rds(fit4, 'Outputs/splevel_alltrays_fit4_noaug.rds')
write_rds(fit5, 'Outputs/splevel_alltrays_fit5_noaug.rds')

