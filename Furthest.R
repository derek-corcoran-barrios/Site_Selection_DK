
# if you already made this layerrs skip to line 34
library(fasterize)
library(janitor)
library(raster)
library(sf)
library(tidyverse)
library(vegan)

Wetlands <- raster("O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/predictors/wetlands/hdr.adf")
Wetlands[Wetlands == 0] <- 1

Layers <- list.files("CHELSA/", full.names = T) %>% 
  purrr::map(raster) %>% 
  purrr::reduce(stack) %>% 
  projectRaster(crs = "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs")
Layers <- Layers[[c(1,4,12,15)]] %>% crop(Wetlands) %>% 
  resample(Wetlands, method ="ngb")
names(Layers) <- c("Temp", "TempSeas", "Prec", "PrecSeas")

Layers <- Layers*Wetlands

Wetness_lidar_files <- list.files(path = "O:/Nat_Ecoinformatics-tmp/au634851/dk_lidar_backup_2021-03-09/outputs/twi/", full.names = T, pattern = ".vrt")
Wetness_lidar_files <- raster(Wetness_lidar_files)%>% 
  projectRaster(crs = "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs") %>% 
  resample(Wetlands, method ="ngb")

Wetness_lidar_files <- Wetness_lidar_files*Wetlands

Layers <- addLayer(Layers, Wetness_lidar_files)

names(Layers)[5] <- "TWI"

# Start here
Layers <- readRDS("Layers_For_Distance.rds")

Alldata <- read_rds("AllData.rds")
Alldata <- Alldata[!is.na(st_dimension(Alldata)),]

RawVars <- raster::extract(Layers, y = Alldata) %>% as.data.frame()

WithVars <- RawVars %>% scale()  %>% as.data.frame() 

## Add a new variable

Canopy_Cover <- list.files(path = "O:/Nat_Ecoinformatics-tmp/au634851/dk_lidar_backup_2021-03-09/outputs/canopy_height/", full.names = T, pattern = ".vrt")

Canopy_Cover <- raster::extract(raster(Canopy_Cover), y = Alldata) %>% as.data.frame() #%>% scale()  %>% as.data.frame()
colnames(Canopy_Cover) <- "Canopy_Height"

RawVars <- cbind(RawVars, Canopy_Cover)
# Bins <1,1-3, 3-10, >10
RawVars <- RawVars %>% 
  mutate(Less_than_1 = ifelse(Canopy_Cover < 100, 1, 0),
         Between_1_and_3 = ifelse(Canopy_Cover >= 100 & Canopy_Cover < 300, 1, 0),
         Between_3_and_10 = ifelse(Canopy_Cover >= 300 & Canopy_Cover< 1000, 1, 0),
         Over_10 = ifelse(Canopy_Cover >= 1000,1,0))

Vegetation_dens <- Wetness_lidar_files <- list.files(path = "O:/Nat_Ecoinformatics-tmp/au634851/dk_lidar_backup_2021-03-09/outputs/proportions/vegetation_density/", full.names = T, pattern = ".vrt")

Vegetation_dens <- raster::extract(raster(Vegetation_dens), y = Alldata) %>% as.data.frame() #%>% scale()  %>% as.data.frame()
colnames(Vegetation_dens) <- "Vegetation_density"
RawVars <- cbind(RawVars, Vegetation_dens)

#WithVars <- readRDS("withVars.rds")

SoilRasters <- list.files(path = "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology",pattern = ".tif", recursive = T, full.names = T)
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = ".aux", negate = T)]
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = ".xml", negate = T)]
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = ".ovr", negate = T)]

SoilRasters <- SoilRasters[c(78,79)]

Names <- c("PhSurface", "phDeep")

Soils <- list()
for(i in 1:length(SoilRasters)){ 
  Soils[[i]] <- raster::extract(raster(SoilRasters[i]), y = Alldata) %>% as.data.frame() 
  colnames(Soils[[i]]) <- Names[i]
  RawVars <- cbind(RawVars, Soils[[i]])
  message(paste(i, "of", length(SoilRasters)))
}

saveRDS(RawVars, "RawVars.rds")

AllData <- Alldata %>% cbind(RawVars) %>% 
  dplyr::filter(!is.na(Temp), !is.na(Canopy_Height), !is.na(phDeep), !is.na(PhSurface))  %>% 
  mutate(Rank = case_when(Dataset %in% c("Biowide", "Agriculture", "Microflora Danica") ~ 0,
                          Dataset %in% c("Novana_Stable", "Novana_Increase", "Novana_Decrease")~ 1)) %>% 
  tibble::rowid_to_column() 

saveRDS(AllData, "AllData2.rds")

WithVars <- RawVars  %>% 
  dplyr::select("Temp", "TempSeas", "Prec", "PrecSeas", "TWI", 
                "Vegetation_density", "PhSurface", "phDeep") %>% 
  scale() %>% 
  as.data.frame() %>% 
  bind_cols(dplyr::select(RawVars, "Less_than_1", "Between_1_and_3", "Between_3_and_10", "Over_10")) %>% 
  dplyr::filter_all(~!is.na(.x))


saveRDS(AllData, "AllData2.rds")
saveRDS(WithVars, "WithVars.rds")

#WithVars <- raster::extract(Bios, Points) %>% as.data.frame() %>% scale()

Dist <- vegan::vegdist(x = as.matrix(WithVars), method = "euclidean") %>% as.matrix() %>% as.data.frame() 
#saveRDS(Dist, "Dist.rds")

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)

ToRank <- 2000

for(i in 1:ToRank){
  Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
  Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)

  Temp <- Dist[Used, Unused]
  rownames(Temp) <- Used
  colnames(Temp) <- Unused
  
  dmax <- max(apply(Temp,2,min,na.rm=TRUE))
  
  Cond <- which(Temp == dmax, arr.ind = TRUE)[1,] %>% as.numeric()
  
  AllData$Rank <- ifelse(AllData$rowid == Unused[Cond[2]], (i + 1), AllData$Rank)
  AllData$Dataset <- ifelse(AllData$rowid == Unused[Cond[2]], "Ranked", AllData$Dataset)
  print(paste(i, "of", ToRank, ", distance =", dmax))
}


OnlyPoints <- AllData %>% dplyr::filter(!is.na(Rank))

Prior <- OnlyPoints %>% dplyr::filter(Rank == 0)

New <- OnlyPoints  %>% dplyr::filter(Rank > 0)

New %>% dplyr::select(Rank) %>% write_sf("Ranked.shp")

Denmark <- read_rds("DK_Shape.rds")

ggplot() + 
  geom_sf(data = Denmark) +
  geom_sf(data = Prior) + 
  geom_sf(data = New, aes(color = Rank)) + 
#  scale_color_gradient(low = "#fff5f0", high = "#cb181d") + 
  scale_colour_viridis_b(option = "C") +
  theme_bw()

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)

RawDist <- vegan::vegdist(x = as.matrix(WithVars[Used,]), method = "euclidean")

saveRDS(RawDist, "RawDist.rds")

nmds = monoMDS(RawDist)# 

nmds_DF <- nmds$points %>% as.data.frame()

OnlyPoints <- cbind(OnlyPoints, nmds_DF)

ForGraph <- OnlyPoints %>% dplyr::filter(MDS1 < 20)

Prior <- ForGraph %>% dplyr::filter(Rank == 0)

New <- ForGraph  %>% dplyr::filter(Rank > 0)

ggplot(Prior, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset)) +
  theme_bw() +
  geom_point(data = New)

Animation <- ggplot(ForGraph, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset)) +
  theme_bw() +
  transition_reveal(along = Rank)

animate(Animation, width = 1100, height = 1100, nframes = 150, renderer = gifski_renderer(loop = F), end_pause = 30, res = 150, fps = 8)
anim_save("Rank.gif")

# plot(Points["Rank"])
# AllData <- Alldata %>% cbind(WithVars) %>% tibble::rowid_to_column() %>% 
#   dplyr::filter_all(~!is.na(.x))
# 
# Soils <- readRDS("Soils.rds")
# names(Soils) <- c("Humus", "TopPH", "DeepPH", "Yield wheat")
# 
# Layers <- readRDS("Layers.rds")
# 
# Height <- readRDS("Height.rds")
# names(Height) <- "Canopy height"
# 
# TWI <- readRDS("Twi.rds")
# names(TWI) <- "TWI"
# 
# Wetlands <- raster("O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/predictors/wetlands/hdr.adf")
# Wetlands[Wetlands == 0] <- NA
# 
# Valldepth <- raster("O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/predictors/valldepth/hdr.adf")
# Valldepth <- Valldepth*Wetlands
# names(Valldepth) <- "Valley depth"
# 
# Vars <- stack(Layers, Soils, Height, TWI, Valldepth)
# 
# 

Dist <- vegan::vegdist(x = as.matrix(WithVars[,-6]), method = "euclidean") %>% as.matrix() %>% as.data.frame() 

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)

ToRank <- 2000

for(i in 1:ToRank){
  Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
  Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)
  
  Temp <- Dist[Used, Unused]
  rownames(Temp) <- Used
  colnames(Temp) <- Unused
  
  dmax <- max(apply(Temp,2,min,na.rm=TRUE))
  
  Cond <- which(Temp == dmax, arr.ind = TRUE)[1,] %>% as.numeric()
  
  AllData$Rank <- ifelse(AllData$rowid == Unused[Cond[2]], (i + 1), AllData$Rank)
  AllData$Dataset <- ifelse(AllData$rowid == Unused[Cond[2]], "Ranked", AllData$Dataset)
  print(paste(i, "of", ToRank, ", distance =", dmax))
}


OnlyPoints <- AllData %>% dplyr::filter(!is.na(Rank))

Prior <- OnlyPoints %>% dplyr::filter(Rank == 0)

New <- OnlyPoints  %>% dplyr::filter(Rank > 0)

New %>% dplyr::select(Rank) %>% write_sf("Ranked2.shp")

Denmark <- read_rds("DK_Shape.rds")

ggplot() + 
  geom_sf(data = Denmark) +
  geom_sf(data = Prior) + 
  geom_sf(data = New, aes(color = Rank)) + 
  #  scale_color_gradient(low = "#fff5f0", high = "#cb181d") + 
  scale_colour_viridis_b(option = "C") +
  theme_bw()

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)

RawDist <- vegan::vegdist(x = as.matrix(WithVars[Used,-6]), method = "euclidean")

saveRDS(RawDist, "RawDist2.rds")

nmds = monoMDS(RawDist)# 

nmds_DF <- nmds$points %>% as.data.frame()

OnlyPoints <- cbind(OnlyPoints, nmds_DF)

saveRDS(OnlyPoints, "Onlypoints.rds")

OnlyPoints <- readRDS("Onlypoints.rds")

ForGraph <- OnlyPoints 

Prior <- ForGraph %>% dplyr::filter(Rank == 0)

New <- ForGraph  %>% dplyr::filter(Rank > 0)

ggplot(Prior, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset)) +
  theme_bw() +
  geom_point(data = New)

First_100 <- Onlypoints %>% 
  dplyr::filter(Rank <= 100)

library(gifski)
library(gganimate)

Animation <- ggplot(Onlypoints, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset, group = seq_along(Rank))) +
  theme_bw() +
  transition_reveal(along = Rank, range = c(1,50)) +
  enter_grow(size = 0.5)

library(gifski)

animate(Animation, width = 1100, height = 1100, nframes = 200, renderer = gifski_renderer(loop = F), end_pause = 30, res = 150, fps = 8)
anim_save("Test.gif")

## option 3

Dist <- vegan::vegdist(x = as.matrix(WithVars), method = "euclidean") %>% as.matrix() %>% as.data.frame() 

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)

ToRank <- 1000

for(i in 1:ToRank){
  
  Temp <- Dist[Used, Unused]
  rownames(Temp) <- Used
  colnames(Temp) <- Unused
  
  dmax <- max(apply(Temp,2,min,na.rm=TRUE))
  
  Cond <- which(Temp == dmax, arr.ind = TRUE)[1,] %>% as.numeric()
  
  AllData$Rank <- ifelse(AllData$rowid == Unused[Cond[2]], (i + 1), AllData$Rank)
  AllData$Dataset <- ifelse(AllData$rowid == Unused[Cond[2]], "Ranked", AllData$Dataset)
  Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
  Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)
  
  saveRDS(Used, "Used.rds")
  saveRDS(Unused, "Unused.rds")
  print(paste(i, "of", ToRank, ", distance =", dmax))
}


OnlyPoints <- AllData %>% dplyr::filter(!is.na(Rank))

Prior <- OnlyPoints %>% dplyr::filter(Rank == 0)

New <- OnlyPoints  %>% dplyr::filter(Rank > 0)

New %>% dplyr::select(Rank) %>% write_sf("Ranked3.shp")

Denmark <- read_rds("DK_Shape.rds")

ggplot() + 
  geom_sf(data = Denmark) +
  geom_sf(data = Prior) + 
  geom_sf(data = New, aes(color = Rank)) + 
  #  scale_color_gradient(low = "#fff5f0", high = "#cb181d") + 
  scale_colour_viridis_b(option = "C") +
  theme_bw()

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)

RawDist <- vegan::vegdist(x = as.matrix(WithVars[Used,]), method = "euclidean")

saveRDS(RawDist, "RawDist3.rds")

nmds = monoMDS(RawDist)# 

nmds_DF <- nmds$points %>% as.data.frame()

OnlyPoints <- cbind(OnlyPoints, nmds_DF)

saveRDS(OnlyPoints, "Onlypoints2.rds")

OnlyPoints <- readRDS("Onlypoints2.rds")

ForGraph <- OnlyPoints 

Prior <- ForGraph %>% dplyr::filter(Rank == 0)

New <- ForGraph  %>% dplyr::filter(Rank > 0)

ggplot(Prior, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset)) +
  theme_bw() +
  geom_point(data = New)

First_100 <- Onlypoints %>% 
  dplyr::filter(Rank <= 100)

library(gifski)
library(gganimate)

Animation <- ggplot(OnlyPoints, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset, group = seq_along(Rank))) +
  theme_bw() +
  transition_reveal(along = Rank, range = c(1,50))

library(gifski)

animate(Animation, width = 1100, height = 1100, nframes = 200, renderer = gifski_renderer(loop = F), end_pause = 30, res = 150, fps = 8)
anim_save("Test.gif")

## Add soil vars and slope and heat index

SoilRasters <- list.files(path = "O:/AUIT_Geodata/Denmark/Natural_ressources/Soil_geology/Texture3D_2014/geotiffs/",pattern = ".tif", recursive = T, full.names = T)
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = ".aux", negate = T)]
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = ".xml", negate = T)]
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = ".ovr", negate = T)]
SoilRasters <- SoilRasters[str_detect(SoilRasters,pattern = "/a", negate = F)]


Soil3D <- SoilRasters[c(5:7,10)] %>% purrr::map(raster) %>% purrr::reduce(stack)

Soil <- raster::extract(Soil3D, y = Alldata) %>% as.data.frame()

RawVars <- cbind(RawVars, Soil)


Slope <- list.files(path = "O:/Nat_Ecoinformatics-tmp/au634851/dk_lidar_backup_2021-03-09/outputs/slope/", full.names = T, pattern = ".vrt")

Slope <- raster::extract(raster(Slope), y = Alldata) %>% as.data.frame() #%>% scale()  %>% as.data.frame()
colnames(Slope) <- "Slope"

RawVars <- cbind(RawVars, Slope)

saveRDS(RawVars, "RawVars.rds")


Heat <- list.files(path = "O:/Nat_Ecoinformatics-tmp/au634851/dk_lidar_backup_2021-03-09/outputs/heat_load_index/", full.names = T, pattern = ".vrt")

Heat <- raster::extract(raster(Heat), y = Alldata) %>% as.data.frame() #%>% scale()  %>% as.data.frame()
colnames(Heat) <- "Heat_Load_Index"

RawVars <- cbind(RawVars, Heat)

saveRDS(RawVars, "RawVars.rds")

AllData <- Alldata %>% cbind(RawVars) %>% 
  dplyr::filter(!is.na(Temp), !is.na(Canopy_Height), !is.na(phDeep), !is.na(PhSurface), !is.na(afsandno), !is.na(Slope))  %>% 
  mutate(Rank = case_when(Dataset %in% c("Biowide", "Agriculture", "Microflora Danica") ~ 0,
                          Dataset %in% c("Novana_Stable", "Novana_Increase", "Novana_Decrease")~ 1)) %>% 
  tibble::rowid_to_column() 

saveRDS(AllData, "AllData3.rds")

WithVars <- RawVars  %>% 
  dplyr::select("Temp", "TempSeas", "Prec", "PrecSeas", "TWI", "Canopy_Height", 
                "Vegetation_density", "PhSurface", "phDeep", "aclaynor", "afsandno", 
                "agsandno", "asiltnor", "Slope", "Heat_Load_Index") %>% 
  scale() %>% 
  as.data.frame() %>% 
  bind_cols(dplyr::select(RawVars, "Less_than_1", "Between_1_and_3", "Between_3_and_10", "Over_10")) %>% 
  dplyr::filter_all(~!is.na(.x))


saveRDS(AllData, "AllData3.rds")
saveRDS(WithVars, "WithVars.rds")

## option 4

Dist <- vegan::vegdist(x = as.matrix(WithVars), method = "euclidean") %>% as.matrix() %>% as.data.frame() 

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)

ToRank <- 3000

for(i in 1:ToRank){
  
  Temp <- Dist[Used, Unused]
  rownames(Temp) <- Used
  colnames(Temp) <- Unused
  
  dmax <- max(apply(Temp,2,min,na.rm=TRUE))
  
  Cond <- which(Temp == dmax, arr.ind = TRUE)[1,] %>% as.numeric()
  
  AllData$Rank <- ifelse(AllData$rowid == Unused[Cond[2]], (i + 1), AllData$Rank)
  AllData$Dataset <- ifelse(AllData$rowid == Unused[Cond[2]], "Ranked", AllData$Dataset)
  Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)
  Unused <- AllData %>% dplyr::filter(is.na(Rank)) %>% pull(rowid)
  
  saveRDS(Used, "Used.rds")
  saveRDS(Unused, "Unused.rds")
  print(paste(i, "of", ToRank, ", distance =", dmax))
}


OnlyPoints <- AllData %>% dplyr::filter(!is.na(Rank))

Prior <- OnlyPoints %>% dplyr::filter(Rank == 0)

New <- OnlyPoints  %>% dplyr::filter(Rank > 0)

New  %>% write_sf("Ranked4.shp")

Denmark <- read_rds("DK_Shape.rds")

ggplot() + 
  geom_sf(data = Denmark) +
  geom_sf(data = Prior) + 
  geom_sf(data = New, aes(color = Rank)) + 
  #  scale_color_gradient(low = "#fff5f0", high = "#cb181d") + 
  scale_colour_viridis_b(option = "C") +
  theme_bw()

Used <- AllData %>% dplyr::filter(!is.na(Rank)) %>% pull(rowid)

RawDist <- vegan::vegdist(x = as.matrix(WithVars[Used,]), method = "euclidean")

saveRDS(RawDist, "RawDist4.rds")

nmds = monoMDS(RawDist)# 

nmds_DF <- nmds$points %>% as.data.frame()

OnlyPoints <- cbind(OnlyPoints, nmds_DF)

saveRDS(OnlyPoints, "Onlypoints4.rds")

OnlyPoints <- readRDS("Onlypoints4.rds")

ForGraph <- OnlyPoints 

Prior <- ForGraph %>% dplyr::filter(Rank == 0)

New <- ForGraph  %>% dplyr::filter(Rank > 0)

ggplot(Prior, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset)) +
  theme_bw() +
  geom_point(data = New)

First_100 <- OnlyPoints %>% 
  dplyr::filter(Rank <= 100)

library(gifski)
library(gganimate)

Animation <- ggplot(OnlyPoints, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Dataset, group = seq_along(Rank))) +
  theme_bw() +
  transition_reveal(along = Rank, range = c(1,50))

library(gifski)

animate(Animation, width = 1100, height = 1100, nframes = 200, renderer = gifski_renderer(loop = F), end_pause = 30, res = 150, fps = 8)
anim_save("Test4.gif")
