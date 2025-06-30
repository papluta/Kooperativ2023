library(terra)
library(sf)
library(exactextractr)

lst = rast("C://Users/patry/Downloads/Kooperativ_MODIS_LST_Mean_2000_2023_DayNight.tif")
lst_dayonly = rast("C://Users/patry/Downloads/Kooperativ_MODIS_LST_Mean_2000_2023(1).tif")
coord  = read.csv('Data/land_rda.csv')

coord$LST_hist <- exactextractr::exact_extract(
  lst, st_buffer(st_as_sf(coord, coords = c("x", "y"),
                          crs = 4326, agr='constant'), 1000), fun = "mean")

coord$LST_hist_day <- exactextractr::exact_extract(
  lst_dayonly, st_buffer(st_as_sf(coord, coords = c("x", "y"),
                          crs = 4326, agr='constant'), 1000), fun = "mean")
plot(lst)

write.csv(coord, file = 'Data/land_env.csv')

plot(coord$LST ~ coord$LST_hist_day)

library(ggplot2)

ggplot(coord, aes(x,y))+
  geom_point(aes(col = Landscape))+
  geom_text(aes(label = Landscape, col = Landscape))

site_order = c(1, 9, 8,4,3,5,6,7,
               10,40,
               13,12,14,17,18,19,20,21,22,15,
               39,36,38,37,35,34,29,
               25,
               33,32,42,31,30,
               28,
               27)

coord$group = NULL
for (i in 1:nrow(coord)) {
  if (coord$Landscape[i] %in% c('NOM09', 'NOM08','NOM04','NOM03','NOM05','NOM06','NOM07')) {
    coord$group[i] = "g1"
  } else if (coord$Landscape[i] %in% c('NOM10','NOM40')) {
    coord$group[i] = 'g2'
  } else if (coord$Landscape[i] %in% c('NOM13','NOM12','NOM14','NOM17','NOM16','NOM18','NOM19','NOM20','NOM21','NOM22','NOM15')) {
    coord$group[i] = 'g3'
  } else if (coord$Landscape[i] %in% c('NOM39','NOM36','NOM38','NOM37','NOM35','NOM34','NOM29')) {
    coord$group[i] = 'g4'
  } else if (coord$Landscape[i] %in% c('NOM25')) {
    coord$group[i] = 'g5'
  } else if (coord$Landscape[i] %in% c('NOM33','NOM32','NOM42','NOM31','NOM30')) {
    coord$group[i] = 'g6'
  } else if (coord$Landscape[i] %in% c('NOM28')) {
    coord$group[i] = 'g7'
  } else if (coord$Landscape[i] %in% c('NOM01')) {
    coord$group[i] = 'g8'
  } else {
    coord$group[i] = 'g9'
  }
}

library(paletteer)
library(dplyr)
col_pal = paletteer_d("ggsci::category10_d3")[c(1:5,7,9,10,6)]
grouping = data.frame(col = col_pal, group = unique(coord$group))

coord = coord %>% left_join(grouping, by = 'group')


ggplot(coord, aes(x,y))+
  geom_point(aes(col = group))+
  geom_text(aes(label = Landscape, col = group))+
  scale_color_manual(values = col_pal)


write.csv(coord, file = 'coord_colors.csv', row.names = F)
