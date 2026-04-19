library(terra)
library(sf)
library(exactextractr)

lst = rast("C://Users/patry/Downloads/Kooperativ_MODIS_LST_Mean_2000_2023.tif")
coord  = read.csv('Data/coordinates.csv')

coord$LST_historical_daytime <- exactextractr::exact_extract(
  lst, st_buffer(st_as_sf(coord, coords = c("x", "y"),
                          crs = 4326, agr='constant'), 1000), fun = "mean")


write.csv(coord, file = 'Data/lst.csv')

