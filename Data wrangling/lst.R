library(terra)
library(sf)
library(exactextractr)

lst = rast("C://Users/patry/Downloads/Kooperativ_LST(1).tif")
coord  = read.csv('Data/land_rda.csv')

#buffer=(seq(200,1000,by=200))
buffer = 1000
# Extract lst with different buffers directly into different columns
for (distance in buffer) {
  col_name <- paste0("LSTC", distance)
  landuse[[col_name]] <- exactextractr::exact_extract(
    lst, st_buffer(st_as_sf(landuse, coords = c("Longitude", "Latitude"),
                            crs = 4326,agr='constant'), distance), fun = "mean") - 273.15
}


coord$LST <- exactextractr::exact_extract(
  lst, st_buffer(st_as_sf(coord, coords = c("x", "y"),
                          crs = 4326, agr='constant'), 1000), fun = "mean") - 273.15
plot(lst)

write.csv(coord, file = 'Data/land_rda.csv')
