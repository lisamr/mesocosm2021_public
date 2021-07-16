#plot out species-level model results from all trays

library(tidyverse)

merged <- read_csv('Outputs/merged_splevelanalysis.csv')
fit <- read_rds('Outputs/splevel_alltrays_fit3.rds')

#===============================================================================
#posterior predictions, plotting prevalence for only 1 species (radish)
#===============================================================================
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
  n_others = unname(quantile(merged_filtered$n_others_z, probs = c(.1,.5, .9))),
  n_radish = seq(min(merged_filtered$n_radish_zsq), max(merged_filtered$n_radish_zsq), length.out = 30)
)


predictions <- matrix(NA, ncol = nrow(newdat), nrow = nrow(B))
for(i in 1:nrow(B)){
  predictions[i,] <- t(as.matrix(newdat) %*% B[i,])
}
predictions <- inv_logit_scaled(predictions)

#summarize predictions
predictionsdf <- cbind.data.frame(
  newdat, 
  data.frame(median = apply(predictions, 2, median), 
             lower = apply(predictions, 2, rethinking::HPDI)[1,],
             upper = apply(predictions, 2, rethinking::HPDI)[2,] ))


#plot predictions against observations. Do for a single species.
pdf('Figures/radish_prevalence.pdf', width = 3, height = 2.5)
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

