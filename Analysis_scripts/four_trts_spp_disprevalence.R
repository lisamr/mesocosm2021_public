# SPECIES-LEVEL DISEASE PREVALENCE. vs. richness across the four treatments

rm(list = ls())
library(cowplot)
library(tidyverse)
library(tidybayes)
library(brms)
library(loo)
library(bayesplot)
library(PNWColors)

#===============================================================================
#functions
#===============================================================================

zscore <- function(x) (x - mean(x)) / (2*sd(x))
anti_zscore <- function(x, z) 2*sd(x)*z + mean(x)
sumtozero <- function(df){
  df$day <- as.factor(df$day)
  contrasts(df$day) <- contr.sum( length(unique( df$day )))
  return(df)
}
source('Analysis_scripts/betabinom_custom_brms.R') #for using beta-binomial likelihood with brms

#===============================================================================
#load data
#===============================================================================

plantdf <- read_csv('Data/plant_level_data.csv')
traydf <- read_csv('Data/tray_level_data.csv')
trts <- readRDS('Data/treatments_list.RDS')

get_datasets <- function(trays, spIDs){
  plantdata <- plantdf %>% 
    filter(trayID %in% trays,
           trayID != 132) %>% 
    filter(spID %in% spIDs) %>% 
    group_by(trayID, species, spID) %>% 
    summarise(I = sum(state0_v2 != 'C' & state_final == 'I', na.rm = T ),
              n = sum(state0_v2 != 'C' & !is.na(state_final), na.rm = T )) %>% 
    ungroup() 
  chosen_trays <- unique(plantdata$trayID)
  traydata <- traydf %>% 
    filter(trayID %in% chosen_trays) %>%
    mutate(richnessz = zscore(richness),
           tempz = zscore(temp_PM)) %>% sumtozero()
  left_join(plantdata, traydata)
}

df_AD <- get_datasets(trts$add_det$trayID, 1:4)
df_SD <- get_datasets(trts$sub_det$trayID, 1:4) 
df_AS <- get_datasets(trts$add_stoch$trayID, 1:6)
df_SS <- get_datasets(c(trts$sub_stoch$trayID, trts$single_species$trayID), 1:6)




#===============================================================================
#run brms model
#===============================================================================

# Beta-binomial model
bf <- bf(formula = I|vint(n) ~ as.factor(day) + tempz + richnessz + 
           (1 | trayID) + (1 + richnessz | spID),
         family = beta_binomial2)
get_prior(bf, data = df_AD)
fit_AD <- brm(bf, data = df_AD,
              prior = c(prior(normal(0, 1.5), class = Intercept),
                        prior(normal(0, 1), class = b),
                        prior(lkj(2), class = cor),
                        prior(exponential(1), class = phi),
                        prior(cauchy(0,2), class = sd)),
              stanvars = stanvars,
              iter = 2000, chains = 4, cores = 4,
              control = list(adapt_delta = 0.95)) 
fit_SD <- update(fit_AD, newdata = df_SD) 
fit_AS <- update(fit_AD, newdata = df_AS) 
fit_SS <- update(fit_AD, newdata = df_SS) 

#stan functions need to be exposed to do predictions. 
#produces lots of warnings/messages; ignore them.
expose_functions(fit_SS, vectorize = T)

#check out fit--they're fine.
pp_check(fit_AD, nsamples = 100)
pp_check(fit_SD, nsamples = 100)
pp_check(fit_AS, nsamples = 100)
pp_check(fit_SS, nsamples = 100)


#===============================================================================
#show effects of richness + intercepts
#===============================================================================
get_variables(fit_AD)

get_effects <- function(fit, modname){
  fit %>% 
    gather_draws(b_Intercept, b_richnessz) %>% 
    rename(term = .variable, main_ef = .value) %>% 
    mutate(term = gsub('^b_', '', term)) %>% 
    left_join(spread_draws(fit, r_spID[spID,term])) %>% 
    mutate(sp_ef = r_spID + main_ef) %>% 
    group_by(term, spID) %>% 
    median_hdci(sp_ef, .width = .9) %>% 
    mutate(nonzero = ifelse(.lower * .upper > 0, 1, 2),
           model = modname)
}

effects_list <- list(
  get_effects(fit_AD, 'AD'), 
  get_effects(fit_SD, 'SD'), 
  get_effects(fit_AS, 'AS'),
  get_effects(fit_SS, 'SS')
)
post_ALL <- bind_rows(effects_list) %>% 
  mutate_at(vars(sp_ef:.width), function(x) round(x, 2))
print(post_ALL, n = Inf)

write_csv(post_ALL, file = 'Outputs/four_trts_spdisprev.csv')



#===============================================================================
#get predictions
#===============================================================================

#average over trayID, temp, day. isolate average and species-specific effects

ppred <- function(fit=fit_AD, d=df_AD){
  post <- posterior_samples(fit)
  R1 <- min(d$richnessz)
  R2 <- max(d$richnessz)
  rseq <- seq(R1, R2, length.out = 25)
  
  #sorry for the repetitive code. I usually predict in the generated quantities block, 
  # but I don't have a streamlined workflow with brms/tidybayes atm and don't feel like
  # figuring out a cleaner way.
  predictions <- list(NULL)
  predictions[[1]] <- sapply(rseq, function(x){
    with(post, inv_logit_scaled(b_Intercept + `r_spID[1,Intercept]` + 
                                       x*(b_richnessz + `r_spID[1,richnessz]`) ))
  }) 
  predictions[[2]] <- sapply(rseq, function(x){
    with(post, inv_logit_scaled(b_Intercept + `r_spID[2,Intercept]` + 
                                       x*(b_richnessz + `r_spID[2,richnessz]`) ))
  }) 
  predictions[[3]] <- sapply(rseq, function(x){
    with(post, inv_logit_scaled(b_Intercept + `r_spID[3,Intercept]` + 
                                       x*(b_richnessz + `r_spID[3,richnessz]`) ))
  }) 
  predictions[[4]] <- sapply(rseq, function(x){
    with(post, inv_logit_scaled(b_Intercept + `r_spID[4,Intercept]` + 
                                       x*(b_richnessz + `r_spID[4,richnessz]`) ))
  }) 
  if(max(d$spID) > 4){
    predictions[[5]] <- sapply(rseq, function(x){
      with(post, inv_logit_scaled(b_Intercept + `r_spID[5,Intercept]` + 
                                         x*(b_richnessz + `r_spID[5,richnessz]`) ))
    }) 
    predictions[[6]] <- sapply(rseq, function(x){
      with(post, inv_logit_scaled(b_Intercept + `r_spID[6,Intercept]` + 
                                         x*(b_richnessz + `r_spID[6,richnessz]`) ))
    }) 
  }
  
  #summarize the posterior predictions
  tmp <- d %>% distinct(trayID, richness)
  richness <- anti_zscore(x = tmp$richness, z = rseq)
  
  #species specific effects
  for(i in 1:length(predictions)){
    median = apply(predictions[[i]], 2, median)
    lower = apply(predictions[[i]], 2, rethinking::HPDI, .9)[1,]
    upper = apply(predictions[[i]], 2, rethinking::HPDI, .9)[2,]
    predictions[[i]] <- cbind(spID = i, richness, median, lower, upper)
  }
  predictions <- as.data.frame(do.call(rbind, predictions))
  
  if(max(d$spID) == 4){
    predictions <- predictions %>% 
      #make predictions span range of data only
      filter(!(spID>=3 & richness<3.9), !(spID ==2 & richness < 2)) 
  }
  
  return(predictions)
}



#get predictions for each model
datalists <- list(df_AD, df_SD, df_AS, df_SS)
fitlist <- list(fit_AD, fit_SD, fit_AS, fit_SS)
postpreds <- map2(fitlist, datalists, ppred)

#add in whether the slope includes zero in the 90% CI
postpreds <- map2(effects_list, postpreds,
  function(x, y){
    tmp <- x %>% 
      filter(term == 'richnessz') %>% 
      select(spID, nonzero)
    y %>% 
      left_join(tmp) 
}) 


#===============================================================================
# plot predictions
#===============================================================================
# ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'none')
sp_pal <- rev(pnw_palette(name = 'Bay', n=6))

#function for plotting
f_plot <- function(i, predictions, d, Title){
  spmax = max(predictions$spID)
  d %>% 
    mutate(richnessplus = richness - .2 + spID*.1) %>% 
    ggplot(.) +
    #geom_ribbon(data=predictions, aes(richness, ymin=lower, ymax=upper), fill='grey', alpha=.2) +
    geom_line(data = predictions, lwd = 1,
              aes(richness, median, color = as.factor(spID)), lty = predictions$nonzero) +
    labs(title = Title) +
    geom_point(aes(richnessplus, I/n, color = as.factor(spID)), size = 1.5, alpha = .8) +
    scale_linetype_manual(values = c(1, 2)) + 
    scale_color_manual(values = sp_pal[1:spmax]) 
}


#make legend for species
sp <- c('Radish', 'Arugula', 'Basil', 'Red romaine', 
        'Green lettuce', 'Butter lettuce')
plegend <- data.frame(spID = 1:6, 
                      Species = factor(sp, levels = sp)) %>% 
  ggplot(., aes(spID, Species, color = Species)) +
  geom_point(size = 2) +
  scale_color_manual(values = sp_pal) +
  theme(legend.position = 'bottom')


#remove predictions for richness 1 and 2 for species 3 and 4 in deterministic trts
p1 <- plot_grid(
  f_plot(1, postpreds[[1]], datalists[[1]], 'Additive, Non-random') + noaxislab,
  f_plot(2, postpreds[[2]], datalists[[2]], 'Substitutive, Non-random') + noaxislab, 
  f_plot(3, postpreds[[3]], datalists[[3]], 'Additive, Random')+ noaxislab, 
  f_plot(4, postpreds[[4]], datalists[[4]], 'Substitutive, Random')+ noaxislab, 
  scale = .9, labels = "auto", hjust = -2
) +
  draw_label('Richness', x = .5, y = 0, vjust = -.5, size = 11) + 
  draw_label('Species-specific disease prevalence', x = .01, y = .5, angle = 90, size = 11)
plot_grid(p1, get_legend(plegend), nrow = 2, rel_heights = c(.9, .1), scale = c(.95, .85))

ggsave('Figures/four_trts_spp_prev.pdf', device = 'pdf', units = 'in', width = 6, height = 3.9)
