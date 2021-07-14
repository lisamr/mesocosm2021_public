# understanding community level changes in disease risk across the four treatments
# use brms to model effects, compare polynomial term or not. don't deal with varying phi.

rm(list = ls())
library(cowplot)
library(tidyverse)
library(brms)
library(loo)
library(bayesplot)

#===============================================================================
#functions
#===============================================================================

zscore <- function(x) (x - mean(x)) / (2*sd(x))

sumtozero <- function(df){
  df$day <- as.factor(df$day)
  contrasts(df$day) <- contr.sum( length(unique( df$day )))
  return(df)
}

#===============================================================================
#load data
#===============================================================================

plantdf <- read_csv('output/plant_level_data.csv')
traydf <- read_csv('output/tray_level_data.csv')
trts <- readRDS('output/treatments_list.RDS')


set.seed(123)
df_AD <- traydf %>% 
  filter(trayID %in% trts$add_det$trayID) %>% 
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) %>% sumtozero()
df_SD <- traydf %>% 
  filter(trayID %in% trts$sub_det$trayID) %>% 
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) %>% sumtozero()
df_AS <- traydf %>% 
  filter(trayID %in% trts$add_stoch$trayID) %>% 
  filter(!trayID %in% sample(1:10, 8)) %>%#need to remove 8 radish trays
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) %>% 
   sumtozero() 
df_SS <- traydf %>% 
  filter(trayID %in% trts$sub_stoch$trayID) %>% 
  filter(!trayID %in% sample(21:30, 9)) %>% #need to remove 9 radish trays
  mutate(richnessz = zscore(richness),
         tempz = zscore(temp_PM)) %>% 
  sumtozero() 



#===============================================================================
#run brms model
#===============================================================================

# binomial likelihood: produces awful fits
# beta binomial model is too difficult to figure out in brms, so using stan.
# n=4
# bf1 <- bf(formula = I_all|trials(n_allP.na) ~  as.factor(day) + tempz + richnessz,
#           family = binomial())
# get_prior(bf1, data = df_AD)
# fit_AD <- brm(bf1, data = df_AD,
#               prior = c(prior(normal(0, 1.5), class = Intercept),
#                         prior(normal(0, 1), class = b)),
#               iter = 2000, chains = 4, cores = 4)
# fit_SD <- update(fit_AD, newdata = sumtozero(df_SD))
# fit_AS <- update(fit_AD, newdata = df_AS)
# fit_SS <- update(fit_AD, newdata = df_SS)
# 
# 
# pp_check(fit_SS) #oof, pretty bad. 



# negative binomial model
bf2 <- bf(formula = I_all ~ as.factor(day) + tempz + richnessz,
          family = negbinomial())
get_prior(bf2, data = df_AD)
fit_AD <- brm(bf2, data = df_AD,
              prior = c( prior(normal(0, 1), class = b)),
              iter = 2000, chains = 4, cores = 4) #richness = pos
fit_SD <- update(fit_AD, newdata = df_SD) #richness = Neg
fit_AS <- update(fit_AD, newdata = df_AS) #richness = pos
fit_SS <- update(fit_AD, newdata = df_SS) #richness = effect crosses zero

#check out fit--they're fine.
pp_check(fit_AD, nsamples = 100)
pp_check(fit_SD, nsamples = 100)
pp_check(fit_AS, nsamples = 100)
pp_check(fit_SS, nsamples = 100)



#===============================================================================
#plot predictions
#===============================================================================

#ggplot theme stuff
theme_set(theme_classic())
noaxislab <- theme(axis.title = element_blank(), plot.title = element_text(hjust = .5))
pal <- rev(wesanderson::wes_palette("Darjeeling1")[c(1,4,2,5)])


#function for summarizing posterior predictions
ppreds <- function(fit, df){
  #get predictions
  post <- posterior_samples(fit)
  R1 <- min(df$richnessz)
  R2 <- max(df$richnessz)
  rseq <- seq(R1, R2, length.out = 10)
  postpred <- sapply(rseq, function(x) exp(post$b_Intercept + x*post$b_richnessz)) 

  #summarize into df
  median = apply(postpred, 2, median)
  lower = apply(postpred, 2, HPDI, .9)[1,]
  upper = apply(postpred, 2, HPDI, .9)[2,]
  richness = anti_zscore(df_AD$richness, rseq)
  return(data.frame(richness, median, lower, upper))
}
f_plot_d <- function(fit, df, Color, Title, Linetype=1){
  tmp <- ppreds(fit, df)
  ggplot(df) +
    geom_line(data=tmp, aes(richness, median), lty = Linetype, lwd = 1.5, color = grey(.6)) +
    geom_ribbon(data=tmp, aes(richness, ymin = lower, ymax = upper), alpha=.2) +
    geom_jitter(aes(richness, I_all),
                height = 0, width = .05, size = 2, alpha = .7, color = Color) +
    labs(y = 'Community infection prevalence', title = Title) 
}


#plot them

plot_grid(
  f_plot_d(fit_AD, df_AD, pal[1], 'Additive, Deterministic')+ noaxislab,  
  f_plot_d(fit_SD, df_SD, pal[2], 'Substitutive, Deterministic')+ noaxislab, 
  f_plot_d(fit_AS, df_AS, pal[3], 'Additive, Stochastic')+ noaxislab, 
  f_plot_d(fit_SS, df_SS, pal[4], 'Substitutive, Stochastic', 2)+ noaxislab, 
  scale = .9, labels = "auto", hjust = -2
) +
  draw_label('Richness', x = .5, y = 0, vjust = -.5, size = 11) + 
  draw_label('Community disease density', x = .01, y = .5, angle = 90, size = 11)
ggsave('figures/four_trts/comm_dens.pdf', device = 'pdf', units = 'mm', width = 181, height = 100)











