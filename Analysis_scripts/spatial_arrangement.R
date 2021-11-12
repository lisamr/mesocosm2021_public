# how are species distributed in the tray with respect to challenged individuals?

rm(list = ls())

plantdf <- read_csv('Data/plant_level_data.csv')
dist_NN <- function(DF){
  # function to find distance to nearest challenged neighbor
  D <- fields::rdist(cbind(DF$x, DF$y)) # distance matrix
  D <- D[,DF$state0_v2 == 'C'] # only keep columns of challenged individuals
  apply(D, 1, min) # find minimum distance. 
}

#remove unanalyzed trays: uninoculated, unwatered, extras
plantdata <- plantdf %>% 
  filter(trayID < 300, trayID != 132, tray_type != 'augmented') 

# for all plants, get distances to nearest challenged individual
# split dataframe by tray, put into list.
plantlist <- split(plantdata, plantdata$trayID)

# get distances to nearest challenged neighbor
for(i in 1:length(plantlist)){
  plantlist[[i]]$nearestC <- dist_NN(plantlist[[i]]) 
} 

plantdf2 <- bind_rows(plantlist) %>% 
  filter(state0_v2 == 'S', !is.na(state_final)) 


# look at spatial patterns?? ==================================================

# distance to nearest challenged host for all non-challenged individuals

p1 <- ggplot(plantdf2, aes(species, nearestC)) +
  geom_jitter(height = 0, width = .2, alpha = .02, color = 'slateblue') +
  geom_boxplot(outlier.colour = NA, fill = NA, notch = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = 'Distance to nearest \nchallenged plant (cm)', 
       x = 'Species',
       title = 'All non-challenged plants') + 
  annotate(geom="text", x=5, y=15, label="n = 32,926 plants",
           color="black", size = 3)
  


# mean distance to nearest challenged host, grouped by tray and species
p2 <- plantdf2 %>% 
  group_by(trayID, species) %>% 
  summarise(mean = mean(nearestC)) %>% 
  ggplot(., aes(species, mean)) +
  geom_jitter(height = 0, width = .2, alpha = .5, color = 'slateblue') +
  geom_boxplot(outlier.colour = NA, fill = NA, notch = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = 'Mean distance to nearest \nchallenged plant (cm)', 
       x = 'Species',
       title = 'Grouped by tray and species')+   
  annotate(geom="text", x=5, y=6.2, size = 3,
           label="n = 487 \ntray-species groups", color="black")


cowplot::plot_grid(p1, p2, labels = 'auto')  
ggsave('Figures/distance_nearestC.pdf', width = 6.5, height = 3.5)
