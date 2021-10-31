# understanding individual-level disease risk using data from all trays
# will include a spatial covariate (distance to nearest challenged neighbor) to control for randomly located inoculum

rm(list = ls())
library(cowplot)
library(tidyverse)
library(brms)

# functions
zscore <- function(x) (x - mean(x)) / (2*sd(x))
anti_zscore <- function(x, z) 2*sd(x)*z + mean(x)
sumtozero <- function(df){
  contrasts(df$day) <- contr.sum( length(unique( df$day )))
  return(df)
}
dist_NN <- function(DF){
  # function to find distance to nearest challenged neighbor
  D <- fields::rdist(cbind(DF$x, DF$y)) # distance matrix
  D <- D[,DF$state0_v2 == 'C'] # only keep columns of challenged individuals
  apply(D, 1, min) # find minimum distance. 
}
source('Analysis_scripts/betabinom_custom_brms.R') #for using beta-binomial likelihood with brms

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'none')

# load data
plantdf <- read_csv('Data/plant_level_data.csv')
traydf <- read_csv('Data/tray_level_data.csv')

#===============================================================================
#wrangle data
#===============================================================================

# get distances to nearest challenged neighbor. 
# even if a plant didn't emerge, inoculum was still added, so keep those in

# split dataframe by tray, put into list.
plantdf_list <- plantdf %>%  
  filter(trayID < 300) %>% #remove uninoculated trays
  split(., .$trayID)

# get distances to nearest challenged neighbor
for(i in 1:length(plantdf_list)){
  plantdf_list[[i]]$nearestC <- dist_NN(plantdf_list[[i]]) 
} 

# filter plant data
plantdf2 <- plantdf_list %>% 
  bind_rows() %>% # bind the list back together
  # do some filtering
  filter(!is.na(state_final), #unemerged or herbivore-eaten plants 
         trayID != 132, # unwatered tray
         state0_v2=='S') # only assessing disease on non-challenged plants

# merge with tray-level varialbes
traydata <- traydf %>% 
  filter(trayID < 300) %>% #remove uninoculated trays
  filter(trayID != 132) %>%  #and remove the unwatered tray
  mutate(n_others = n_allP.na - n_radish, 
         n_others2 = n_allP.na - n_radish - n_arugula) %>% 
  mutate_at(vars(temp_PM, richness, n_others, n_others2, n_allP.na), list(z=zscore)) %>% 
  mutate_at(vars(n_radish:n_butter), list(zsq = function(x) zscore(sqrt(x)))) %>% 
  mutate_at(vars(n_radish:n_butter), list(p = function(x) x/.$n_allP.na))

# final dataset
merged <- left_join(plantdf2, 
                    traydata %>% select(trayID, day, temp_PM_z:n_butter_zsq), 
                    by = 'trayID') %>% 
  mutate(I = as.numeric(state_final == 'I')) %>% 
  mutate_at(vars(trayID, day, spID), as.factor) 
merged <- sumtozero(merged)


#===============================================================================
# explore distance functions
#===============================================================================

# seems reasonable to have the effect of nearest neighbor decay with 
# inverse distance relationship: ex/ 1/(x^y)
f <- function(x) 1/(x^.5) #could vary up power parameter to see if it affects results
x <- seq(0, 10, length = 100)
plot(x, f(x), type = 'l')

merged <- merged %>% 
  mutate(nearestC_transz = zscore(nearestC^-.5) ) #transformed (d^-.5) and zscored



#===============================================================================
#run individual-level brms model
#===============================================================================

# hopes: results don't change when I include or exclude the 'extra' trays. 
# ideally, i'd just remove them from the anlaysis since they're confusing. 
# richness should not have a positive effect in models including radish.
# also need to check species-specific effects of richness, not just average



# model   tray cov              
# 0:      --                     
# 1:      nall                   
# 2:      nradish
# 3:      nothers, nradish
# 4:      nothers2, nradish, narugula
# 5:      nradish, narugula...nbutter
# 6:      nradish, narugula
bf0 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz +  (1|trayID) + (1 + richness_z | spID))
bf1 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz + n_allP.na_z + (1|trayID) + (1 + richness_z | spID))
bf2 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz + n_radish_zsq + (1|trayID) + (1 + richness_z | spID))
bf3 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz + n_others_z + n_radish_zsq + (1|trayID) + (1 + richness_z | spID))
bf4 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz + n_others2_z + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))
bf5 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz + n_radish_zsq + n_arugula_zsq + n_basil_zsq + n_red_rom_zsq + n_green_rom_zsq + n_butter_zsq + (1|trayID) + (1 + richness_z | spID))
bf6 <- bf(formula = I ~ day + temp_PM_z + richness_z + nearestC_transz + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))

get_prior(bf0, data = merged)

set.seed(2021)
randomtrays <- sample(unique(merged$trayID), 30)
merged_filtered <- merged %>% 
  filter(trayID %in% randomtrays)
merged_noaug <- merged %>% 
  filter(tray_type != 'augmented')

fit1 <- brm(bf1, data = merged_noaug,
            prior = c(prior(normal(0, 1.5), class = Intercept),
                      prior(normal(0, 1), class = b),
                      prior(lkj(2), class = cor),
                      prior(normal(0,2), class = sd)),
            family = bernoulli(),
            iter = 2000, chains = 3, cores = 3, 
            control = list(adapt_delta = 0.95), 
            save_pars = save_pars(all = TRUE))
fit2 <- update(fit1, formula. = bf2, newdata = merged_noaug)
fit5 <- update(fit1, formula. = bf5, newdata = merged_noaug)
fit0 <- update(fit1, formula. = bf0, newdata = merged_noaug)


#===============================================================================
#compare models
#===============================================================================
l1 <- loo(fit1)
l2 <- loo(fit2)
l5 <- loo(fit5)
l0 <- loo(fit0)
loo_compare(l0, l1, l2, l5)

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
bf6 <- bf(formula = I|vint(n) ~ day + temp_PM_z + richness_z + n_radish_zsq + n_arugula_zsq + (1|trayID) + (1 + richness_z | spID))

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
fit6 <- update(fit0, formula. = bf6, newdata = merged, save_pars = save_pars(all = TRUE))

# write_rds(fit3, 'Outputs/splevel_alltrays_fit3.rds')
# write_rds(fit4, 'Outputs/splevel_alltrays_fit4.rds')
# write_rds(fit5, 'Outputs/splevel_alltrays_fit5.rds')
# write_csv(merged, 'Outputs/merged_splevelanalysis.csv') 
