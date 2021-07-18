df <- read_csv('data/plant_level_data.csv')
head(df)

tmp <- df %>% 
  filter(tray_type %in% c('regular', 'augmented')) %>% 
  group_by(trayID, tray_type, richness) %>% 
  summarise(rad_prev = sum(spID == 1) / n(), 
            rad_dens = sum(spID == 1)) 

#overall correlation??
cor(tmp$richness, tmp$rad_prev) # -0.3255051
cor(tmp$richness, tmp$rad_dens) # -0.09505544

tmp2 <- tmp %>% filter(tray_type == 'regular')
#barely made a dent, but whatever.
cor(tmp2$richness, tmp2$rad_prev) #-0.3594578. 
cor(tmp2$richness, tmp2$rad_dens) # -0.1589976



ggplot(tmp, aes(richness, rad_prev, group = tray_type)) +
  geom_jitter(height= 0,width = .1, alpha = .7, aes(color = tray_type)) +
  labs(x = 'Richness', y = 'Proportion radish', color = 'Tray type') +
  theme(legend.position = 'bottom')

ggsave('Figures/prop_rad_supplementary.pdf', width = 4, height = 3)


ggplot(tmp, aes(richness, rad_dens, group = tray_type)) +
  geom_jitter(height= 0,width = .1, alpha = .7, aes(color = tray_type)) +
  labs(x = 'Richness', y = 'Proportion radish', color = 'Tray type') +
  theme(legend.position = 'bottom')
