df <- read_csv('data/plant_level_data.csv')
head(df)

df %>% 
  filter(tray_type %in% c('regular', 'augmented')) %>% 
  group_by(trayID, tray_type, richness) %>% 
  summarise(rad_prev = sum(spID == 1) / n()) %>% 
  ggplot(., aes(richness, rad_prev, group = tray_type)) +
  geom_jitter(height= 0,width = .1, alpha = .7, aes(color = tray_type)) +
  labs(x = 'Richness', y = 'Proportion radish', color = 'Tray type') +
  theme(legend.position = 'bottom')

ggsave('Figures/prop_rad_supplementary.pdf', width = 4, height = 3)
