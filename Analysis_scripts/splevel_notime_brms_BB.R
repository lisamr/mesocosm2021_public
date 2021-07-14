# understanding individual-level disease risk in the four different treatments
# models run seperately for each treatment


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
source('Stan_models/betabinom_custom_brms.R') #for using beta-binomial likelihood with brms

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'none')

#===============================================================================
#load data
#===============================================================================
#M0: zscored(sqrt(radish density))
#M1: zscored(density of others) + zscored(sqrt(radish density))
#M2: zscored(density of others) + zscored(sqrt(radish density)) + zscored(sqrt(arugula density))*radish

plantdf <- read_csv('output/plant_level_data.csv')
traydf <- read_csv('output/tray_level_data.csv')
trts <- readRDS('output/treatments_list.RDS')

#merge data sets
plantdata <- plantdf %>% 
  filter(trayID < 300) %>% #remove uninoculated trays
  filter(trayID != 132) %>%  #and remove the unwatered tray
  group_by(trayID, species, spID) %>% 
  summarise(I = sum(state0_v2=='S' & state_final == 'I', na.rm = T ),
            n = sum(state0_v2=='S' & !is.na(state_final), na.rm = T )) %>% 
  ungroup() 

traydata <- traydf %>% 
  filter(trayID < 300) %>% #remove uninoculated trays
  filter(trayID != 132) %>%  #and remove the unwatered tray
  mutate(n_others = n_allP.na - n_radish, 
         n_others2 = n_allP.na - n_radish - n_arugula, 
         radish_10 = as.numeric(n_radish >= 10)) %>% 
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
# 6:      nradish*(narugula...nbutter)
# 7:      nradish, narugula

bf0 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + (1|trayID) + (1 + richness_z | spID))
bf1 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_allP.na_z + (1|trayID) + (1 + richness_z | spID))
bf2 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + (1|trayID) + (1 + richness_z | spID))
bf3 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_others_z + n_radish_zsq + (1|trayID) + (1 + richness_z | spID))
bf4 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_others2_z + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))
bf5 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + n_arugula_zsq + n_basil_zsq + n_red_rom_zsq + n_green_rom_zsq + n_butter_zsq + (1|trayID) + (1 + richness_z | spID))
bf6 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + n_arugula_zsq + n_basil_zsq + n_red_rom_zsq + n_green_rom_zsq + n_butter_zsq + n_green_rom_zsq:radish_10 + n_butter_zsq:radish_10 + (1|trayID) + (1 + richness_z | spID))
bf7 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))

fit0 <- brm(bf0, data = merged,
            prior = c(prior(normal(0, 1.5), class = Intercept),
                      prior(normal(0, 1), class = b),
                      prior(lkj(2), class = cor),
                      prior(normal(0,1), class = sd)),
            family = beta_binomial2, stanvars = stanvars,
            iter = 2000, chains = 1, cores = 1, 
            control = list(adapt_delta = 0.95, max_treedepth = 15), 
            save_pars = save_pars(all = TRUE))
fit1 <- update(fit0, formula. = bf1, newdata = merged, save_pars = save_pars(all = TRUE))
fit2 <- update(fit0, formula. = bf2, newdata = merged, save_pars = save_pars(all = TRUE))
fit3 <- update(fit0, formula. = bf3, newdata = merged, save_pars = save_pars(all = TRUE))
fit4 <- update(fit0, formula. = bf4, newdata = merged, save_pars = save_pars(all = TRUE))
fit5 <- update(fit0, formula. = bf5, newdata = merged, save_pars = save_pars(all = TRUE))
fit6 <- update(fit0, formula. = bf6, newdata = merged, save_pars = save_pars(all = TRUE))
fit7 <- update(fit0, formula. = bf7, newdata = merged, save_pars = save_pars(all = TRUE),
               iter = 2000, chains = 4, cores = 4)



#===============================================================================
#post-process--checking and comparing models
#===============================================================================
#stan functions need to be exposed for each fit. produces lots of warnings, supress them.
expose_functions(fit0, vectorize = T)

pp_check(fit7, nsamples = 50)

# a0 + 
fit5 %>% 
  spread_draws(b_Intercept) %>% 
  median_hdci()

K0 <- kfold(fit0, K = 10)
K1 <- kfold(fit1, K = 10)
K2 <- kfold(fit2, K = 10)
K3 <- kfold(fit3, K = 10)
K4 <- kfold(fit4, K = 10)
K5 <- kfold(fit5, K = 10)
K6 <- kfold(fit6, K = 10)
K7 <- kfold(fit7, K = 10)
loo_compare(K0, K1, K2, K3, K4, K5, K6, K7)
loo_compare(K0, K1, K2, K3, K4, K5, K7)
# elpd_diff se_diff
# fit7   0.0       0.0  
# fit3  -5.7       5.5  
# fit5  -9.1       5.0  
# fit2 -15.7       7.3  
# fit4 -16.4      10.1  
# fit0 -32.9       9.4  
# fit1 -57.9      10.1  
loo_compare(K0, K1, K2, K3, K4, K5)
# elpd_diff se_diff
# fit3   0.0       0.0  
# fit5  -3.4       5.9  
# fit2 -10.0       7.1  
# fit4 -10.7      10.4  
# fit0 -27.2      10.1  
# fit1 -52.2      10.6 
loo_compare(K0, K1, K3, K4, K5)
# elpd_diff se_diff
# fit3   0.0       0.0  
# fit5  -3.4       5.9  
# fit4 -10.7      10.4  
# fit0 -27.2      10.1  
# fit1 -52.2      10.6  
#===============================================================================
#posterior predictions, plotting for only 1 species
#===============================================================================
fit <- fit3
get_variables(fit)

merged_filtered <- merged %>% filter(spID == 1)
#get posteriors in a matrix
B <- fit %>% 
  spread_draws(b_Intercept, `r_spID[1,Intercept]`, b_richness_z, b_n_others_z, b_n_radish_zsq) %>%
  mutate(b_Intercept = b_Intercept + `r_spID[1,Intercept]`) %>% 
  select(!c(.chain:.draw, `r_spID[1,Intercept]`)) %>% as.matrix
B

#multiply by new dataframe
newdat <- expand_grid(
  a0 = 1, 
  richness = 0,
  n_others = unname(quantile(merged_filtered$n_others_z, probs = c(0,.5, 1))),
  #n_arugula = unname(quantile(merged$n_radish_zsq, probs = c(0,.5, 1))),
  n_radish = seq(min(merged_filtered$n_radish_zsq), max(merged_filtered$n_radish_zsq), length.out = 30)
)


predictions <- matrix(NA, ncol = nrow(newdat), nrow = nrow(B))
for(i in 1:nrow(B)){
  predictions[i,] <- t(as.matrix(newdat) %*% B[i,])
}
predictions <- inv_logit(predictions)

#summarize predictions
predictionsdf <- cbind.data.frame(newdat, 
                 data.frame(median = apply(predictions, 2, median), 
                            lower = apply(predictions, 2, HPDI)[1,],
                            upper = apply(predictions, 2, HPDI)[2,] ))


#plot predictions against observations. Do for a single species.
pdf('figures/radish_prevalence.pdf', width = 3, height = 2.5)
ggplot(predictionsdf) +
  geom_ribbon(aes(n_radish, ymin = lower, ymax = upper, group = n_others, fill = n_others), alpha = .2) +
  geom_line(aes(n_radish, median, group = n_others, color = n_others)) +
  geom_point(data = merged_filtered, aes(n_radish_zsq, I/n, color = n_others_z), size = 1) +
  scale_color_viridis_c(option = 'D', direction = -1) +
  scale_fill_viridis_c(option = 'D', direction = -1) +
  labs(color = 'density of \nother species', fill = 'density of \nother species',
       y = 'Disease prevalence of radish', x = 'Density of radish') + 
  theme(legend.position = 'bottom', legend.key.size = unit(.1, 'in'), 
        text = element_text(size = 7), legend.text = element_text(size = 4), 
        legend.title = element_text(size = 5))
dev.off()


#===============================================================================
#posterior predictions, plotting for only 1 species
# Redo again with fit7.
#===============================================================================
fit <- fit7
merged_filtered <- merged %>% filter(spID == 1) 
 B <- fit %>% 
   spread_draws(b_Intercept, `r_spID[1,Intercept]`, b_richness_z, b_n_arugula_zsq, b_n_radish_zsq) %>%
   mutate(b_Intercept = b_Intercept + `r_spID[1,Intercept]`) %>% 
   select(!c(.chain:.draw, `r_spID[1,Intercept]`)) %>% as.matrix
 
 #multiply by new dataframe
 newdat <- expand_grid(
   a0 = 1, 
   richness = 0,
   n_arugula = unname(quantile(merged_filtered$n_arugula_zsq, probs = c(0,.5, 1))),
   n_radish = seq(min(merged_filtered$n_radish_zsq), max(merged_filtered$n_radish_zsq), length.out = 30)
 )
 
 
 predictions <- matrix(NA, ncol = nrow(newdat), nrow = nrow(B))
 for(i in 1:nrow(B)){
   predictions[i,] <- t(as.matrix(newdat) %*% B[i,])
 }
 predictions <- inv_logit(predictions)
 
 #summarize predictions
 predictionsdf <- cbind.data.frame(newdat, 
                                   data.frame(median = apply(predictions, 2, median), 
                                              lower = apply(predictions, 2, HPDI)[1,],
                                              upper = apply(predictions, 2, HPDI)[2,] ))

ggplot(predictionsdf) +
  geom_ribbon(aes(n_radish, ymin = lower, ymax = upper, group = n_arugula, fill = n_arugula), alpha = .2) +
  geom_line(aes(n_radish, median, group = n_arugula, color = n_arugula)) +
  geom_point(data = merged_filtered, aes(n_radish_zsq, I/n, color = n_arugula_zsq), size = 2) +
  scale_color_viridis_c(option = 'D', direction = -1) +
  scale_fill_viridis_c(option = 'D', direction = -1)
