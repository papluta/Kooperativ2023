library(dplyr)
library(tidyverse)

c500.raw = read.csv('Data/Kooperativ2023_laduse_500m_circle.csv')
h600.raw = read.csv('Data/Kooperativ2023_laduse_hexagons.csv')

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
