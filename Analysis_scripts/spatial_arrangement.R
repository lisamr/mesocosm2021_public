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
