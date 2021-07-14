plantdf <- read_csv('output/plant_level_data.csv')
traydf <- read_csv('output/tray_level_data.csv')


#merge data sets
monos <- traydf %>% 
  filter(trayID < 300) %>% #remove uninoculated trays
  filter(richness ==1, n_allP == 238) %>% 
  pull(trayID)

plantdf %>% 
  filter(trayID %in% monos) %>% 
  group_by(spID, trayID) %>% 
  filter(state0_v2 != 'C', !is.na(state_final)) %>% 
  summarise(prev = sum(state_final == 'I') / n()) %>% 
  group_by(spID) %>% 
  summarise(mean = mean(prev), sd = sd(prev))
