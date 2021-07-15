# COMMUNITY-LEVEL DISEASE PREVALENCE vs. richness across the four treatments


rm(list = ls())
library(cowplot)
library(tidyverse)
library(tidybayes)
library(brms)
library(loo)
library(bayesplot)

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

#get datasets for each of the four treatments.
df_AD <- traydf %>% 
  filter(trayID %in% trts$add_det$trayID,
         trayID != 132) %>% #remove tray132. missed a water day. 
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) %>% sumtozero()
df_SD <- traydf %>% 
  filter(trayID %in% trts$sub_det$trayID) %>% 
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) %>% sumtozero()

# for the 'stochastic' (random) species loss, extra steps are involved.
# We only need 1 or 2 monospecific radish trays, but we have data for 10.
# Resample radish trays until all permutations exhausted, run models on each 
# of those datasets, combine posterior.
# It's a brute-force way to run a weighted regression, but the models are fast.
perms_AS <- combn(1:10, 8) #need to remove 8 radish trays, leaving 2 for each dataset
df_AS <- list(NULL)
for(i in 1:ncol(perms_AS)){
  df_AS[[i]] <- traydf %>% 
    filter(trayID %in% trts$add_stoch$trayID) %>% 
    filter(!trayID %in% perms_AS[,i]) %>%
    mutate(richnessz = zscore(richness),
           tempz = zscore(temp_PM)) %>% 
    sumtozero() 
} 
perms_SS <- combn(21:30, 9) #need to remove 9 radish trays, leaving 1 for each dataset
df_SS <- list(NULL)
for(i in 1:ncol(perms_SS)){
  df_SS[[i]] <- traydf %>% 
    filter(trayID %in% trts$sub_stoch$trayID) %>% 
    filter(!trayID %in% perms_SS[,i]) %>% #need to remove 9 radish trays
    mutate(richnessz = zscore(richness),
           tempz = zscore(temp_PM)) %>% 
    sumtozero() 
} 


#===============================================================================
#run brms model
#===============================================================================

# Beta-binomial model
bf <- bf(formula = I_all|vint(n_allP.na) ~ as.factor(day) + tempz + richnessz,
         family = beta_binomial2)
get_prior(bf, data = df_AD)
fit_AD <- brm(bf, data = df_AD,
              prior = c(prior(normal(0, 1.5), class = Intercept),
                        prior(normal(0, 1), class = b)),
              stanvars = stanvars,
              iter = 2000, chains = 4, cores = 4) 
fit_SD <- update(fit_AD, newdata = df_SD) 
fit_AS <- lapply(df_AS, function(x) update(fit_AD, newdata = x)) 
fit_SS <- lapply(df_SS, function(x) update(fit_AD, newdata = x)) 

#stan functions need to be exposed to do predictions. 
#produces lots of warnings/messages; ignore them.
expose_functions(fit_AD, vectorize = T)

#check out fit--they're fine.
pp_check(fit_AD, nsamples = 100)
pp_check(fit_SD, nsamples = 100)
pp_check(fit_AS[[1]], nsamples = 100)
pp_check(fit_SS[[1]], nsamples = 100)


#===============================================================================
#show effects of richness
#===============================================================================

post_AD <- fit_AD %>% 
  spread_draws(b_richnessz) %>% 
  median_hdci(.width = .9) %>% 
  mutate(model = 'AD')
post_SD <- fit_SD %>% 
  spread_draws(b_richnessz) %>% 
  median_hdci(.width = .9)%>% 
  mutate(model = 'SD')
post_AS <- lapply(fit_AS, function(x) spread_draws(x, b_richnessz)) %>% 
  bind_rows() %>% 
  median_hdci(.width = .9)%>% 
  mutate(model = 'AS')
post_SS <- lapply(fit_SS, function(x) spread_draws(x, b_richnessz)) %>% 
  bind_rows() %>% 
  median_hdci(.width = .9)%>% 
  mutate(model = 'SS')
post_ALL <- bind_rows(post_AD, post_SD, post_AS, post_SS) %>% 
  mutate_at(vars(b_richnessz:.width), function(x) round(x, 2))

write_csv(post_ALL, file = 'Outputs/four_trts_comdisprev.csv')



#===============================================================================
#plot predictions
#===============================================================================

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5))
pal <- rev(wesanderson::wes_palette("Darjeeling1")[c(1,4,2,5)])


#function for summarizing posterior predictions
ppreds <- function(post, df){
  #get predictions
  R1 <- min(df$richnessz)
  R2 <- max(df$richnessz)
  rseq <- seq(R1, R2, length.out = 10)
  postpred <- sapply(rseq, function(x) inv_logit_scaled(post$b_Intercept + x*post$b_richnessz)) 
  
  #summarize into df
  median = apply(postpred, 2, median)
  lower = apply(postpred, 2, rethinking::HPDI, .9)[1,]
  upper = apply(postpred, 2, rethinking::HPDI, .9)[2,]
  richness = anti_zscore(df$richness, rseq)
  return(data.frame(richness, median, lower, upper))
}

f_plot_d <- function(predictions, df, Color, Title, Linetype=1){
  ggplot(df) +
    geom_line(data=predictions, aes(richness, median), 
              lty = Linetype, lwd = 1, color = grey(.6)) +
    geom_ribbon(data=predictions, aes(richness, ymin = lower, ymax = upper), alpha=.2) +
    geom_jitter(aes(richness, I_all / n_allP.na),
                height = 0, width = .1, size = 2, alpha = .7, color = Color) +
    labs(y = 'Community infection prevalence', title = Title) +
    scale_y_continuous(limits = c(0,1))
}


#plot them
pred_AD <- ppreds(posterior_samples(fit_AD), df_AD)
pred_SD <- ppreds(posterior_samples(fit_SD), df_SD)
pred_AS <- lapply(fit_AS, function(x) posterior_samples(x)) %>% 
  bind_rows() %>% ppreds(., df_AS[[1]])
pred_SS <- lapply(fit_SS, function(x) posterior_samples(x)) %>% 
  bind_rows() %>% ppreds(., df_SS[[1]])
plot_grid(
  f_plot_d(pred_AD, df_AD, pal[1], 'Additive, Non-random')+ noaxislab,  
  f_plot_d(pred_SD, df_SD, pal[2], 'Substitutive, Non-random')+ noaxislab, 
  f_plot_d(pred_AS, df_AS[[1]], pal[3], 'Additive, Random')+ noaxislab, 
  f_plot_d(pred_SS, df_SS[[1]], pal[4], 'Substitutive, Random', 2)+ noaxislab, 
  scale = .9, labels = "auto", hjust = -2
) +
  draw_label('Richness', x = .5, y = 0, vjust = -.5, size = 11) + 
  draw_label('Community disease prevalence', x = .01, y = .5, angle = 90, size = 11)
ggsave('Figures/four_trts_comm_prev.pdf', device = 'pdf', units = 'mm', width = 181, height = 100)











