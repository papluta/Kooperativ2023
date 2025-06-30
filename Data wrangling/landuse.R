library(dplyr)
library(tidyverse)


c1000.raw = read.csv('Data/gis/Kooperativ_1000_raw_final.csv')

c1000 = c1000.raw %>% group_by(Landscape, jj_landuse4) %>%
  summarise(area = sum(area)) %>% mutate(area_h = area/10000) %>%
  mutate(proportion = area/(pi*1000^2), landuse = as.factor(gsub(' ','_', jj_landuse4))) %>% 
  dplyr::select(-jj_landuse4)

unique(c1000$landuse)

c1000.t = c1000 %>% ungroup() %>% dplyr::select(-c(area, area_h)) %>% 
  pivot_wider(names_from = landuse, values_from = proportion, values_fill = 0) 

sum = data.frame( sum = rowSums(c1000.t[,-1]))

str(c1000.t)

landuse1000 = c1000.t %>% mutate(SNH_ann = Kooperativ_flower_field + annual_flower_field + 0.5 * structural_flower_field,
                                 SNH_per = AESother + SNHopen + perennial_flower_field + 0.5 * structural_flower_field,
                                 Grassland = grassland,
                                 Grassland_organic = grassland_AKUM + grassland_organic,
                                 Unclear = unclear + `NA`
                                 ) %>%
  dplyr::select(Landscape, SNH_ann, SNH_per, Grassland, Grassland_organic, Crop = crop, Crop_organic = crop_organic, 
                Forest = forest, SNH_wood = SNHwoody, Settlement = settlement, Water = water, Unclear)

data.frame( sum = rowSums(landuse1000[,-1]))

library(vegan)
rownames(landuse1000) = landuse1000$Landscape

landuse1000$Shannon = diversity(landuse1000[,-1], index = 'shannon') # no "unklar"

write.csv(landuse1000, file = "Data/Kooperativ_landuse_1000_250613.csv")

######## OTHER RADII
h600.raw = read.csv('Data/Kooperativ2023_laduse_hexagons.csv')
c500.raw = read.csv('Data/Kooperativ2023_laduse_500m_circle_raw.csv')


c500 = c500.raw %>% group_by(Landscape, landuse) %>%
  summarise(area = sum(AREA)) %>% mutate(area_h = area/10000) %>%
  mutate(proportion = area/(pi*500^2), landuse = as.factor(gsub(' ','_', landuse)))

unique(c500$landuse)

c500.t = c500 %>% ungroup() %>% select(-c(area, area_h)) %>% 
  pivot_wider(names_from = landuse, values_from = proportion, values_fill = 0)

h600 = h600.raw %>% group_by(Landscape, landuse) %>%
  summarise(area = sum(AREA)) %>% mutate(area_h = area/10000) %>%
  mutate(landuse = as.factor(gsub(' ','_', landuse)))

h600.t = h600 %>% ungroup() %>% select(-c(area))  %>%
  pivot_wider(names_from = landuse, values_from = area_h, values_fill = 0) %>%
  mutate(structural_flower_field = structural_flower_field + `flower_field?`) %>% select(-`flower_field?`)



library(vegan)
rownames(h600.t) = h600.t$Landscape

h600.t$Shannon = diversity(h600.t[,-c(1,18)], index = 'shannon') # no "unklar"

h600.t2 = h600.t %>%
  mutate(SNH_combined = AESothers + SNHopen + SNHwoody,
         Flowers_combined = Kooperativ_flower_field + annual_flower_field + perennial_flower_field + structural_flower_field,
         #Other_combined = others + unklar + settlement + water,
         Grassland_combined = grassland + grassland_organic + grassland_AKUM) %>%
  rename(Crop_organic = organic_crop, Crop_conventional = crops, Forest = forest)

write.csv(h600.t2, 'Data/Kooperativ2023_laduse_600m_hexagons.csv', row.names = F)
