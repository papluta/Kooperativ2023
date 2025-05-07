library(dplyr)
library(tidyverse)
library(geosphere)

#colony.ab = read.csv('Results/capwire_colony_ab.csv')
land = read.csv('Data/gis/Kooperativ_landuse1000.csv')
coor = read.csv('Data/Kooperativ2023_land.csv')

land.rda = land %>% select(Landscape, crop) %>% left_join(coor %>% select(Landscape = Landscape_ID, x, y), by = 'Landscape')

write.csv(land.rda, file = 'Data/land_rda.csv', row.names = F)

he = read.csv('Data/Goe_filtered095_Het60_pruned_heterozygosity.csv')
fst = read.csv('Data/Goe_filtered095_Het60_pruned_genetic_distances_bootstraps.csv')[, c(2,3,107)] %>% 
  filter(Population1 != 'NOM16' & Population2 != 'NOM16') 
pop = read.table('Data/populations.txt')
# samples_rel = read.csv('Data/my_project.BestConfig', sep = '') %>% mutate(id = sub('_.*', '', OffspringID)) %>% 
#   left_join(pop, by = join_by('id' == 'V1')) 
# samples_not_rel = read.csv('Data/samples_NR.txt', sep = '', header = F) %>% mutate(id = sub('_.*', '', V1)) %>% 
#   left_join(pop, by = join_by('id' == 'V1')) 
# n.r = samples_rel %>%
#   group_by(V2) %>% summarise(n.related = n()) %>% rename(Landscape = V2)
# n.nr = samples_not_rel %>%
#   group_by(V2) %>% summarise(n.not.related = n()) %>% rename(Landscape = V2)
coord = read.csv('Data/Kooperativ2023_land.csv', header = T) %>% 
  filter(Landscape_ID != 'NOM16' &  Landscape_ID != 'NOM11' &  Landscape_ID != 'NOM24')

n_caught = read.csv('Data/Kooperativ2023_samples.csv', header = T) %>% 
  filter(Landscape_ID != 'NOM11' &  Landscape_ID != 'NOM24') %>% group_by(Landscape_ID) %>% summarise(n = n())
mean(n_caught$n)
# Compute pairwise geographic distances (great-circle distances in km)
geo_dist_matrix <- distm(coord[, c("x", "y")]) / 1000 
rownames(geo_dist_matrix) <- coord$Landscape_ID
colnames(geo_dist_matrix) <- coord$Landscape_ID

fst2 = fst %>% 
  pivot_longer(Population1:Population2, names_to = 'pop', values_to = 'Landscape') %>%
  group_by(Landscape) %>%
  summarise(Fst.m = mean(Fst))

mean(fst$Fst)

# data = colony.ab %>% dplyr::select(pop, ml.pop.size) %>% rename(Landscape = pop, Pop.size = ml.pop.size) %>%
#   mutate(Colony.dens = Pop.size/((pi*(600+775)^2)/1000000)) %>%
#   mutate(Colony.dens1km = Pop.size/((pi*(1000)^2)/1000000)) %>%
#   left_join(n.r, by = 'Landscape') %>%
#   left_join(n.nr, by = 'Landscape') %>%
#   left_join(he %>% dplyr::select(pop, uHe, uHeSD), by = join_by('Landscape' == 'pop')) %>%
#   left_join(fst2, by = 'Landscape') %>% left_join(land, by = 'Landscape') %>% 
#   mutate(Good_habitat = Crop_organic + AESothers + Kooperativ_flower_field + SNHopen + SNHwoody + 
#            grassland_organic + annual_flower_field + grassland_AKUM + perennial_flower_field + 
#            structural_flower_field) %>% 
#   filter(Landscape != 'NOM16') 


data = he %>% dplyr::select(pop, uHe, uHeSD) %>% rename(Landscape = pop) %>%
  left_join(fst2, by = 'Landscape') %>% left_join(land, by = 'Landscape') %>% 
  mutate(Good_habitat = crop_organic + AESother + Kooperativ_flower_field + SNHopen + SNHwoody + 
           grassland_organic + annual_flower_field + grassland_AKUM + perennial_flower_field + 
           structural_flower_field) %>% 
  filter(Landscape != 'NOM16') 

colnames(data)
data.clean = data %>% drop_na(Pop.size)  

mean(data.clean$Colony.dens)
sd(data.clean$Colony.dens)

write.csv(data[,1:8], file = 'Results/gen_div_table.csv', row.names = F)
library(nlme)
library(DHARMa)
library(effects)
library(lme4)
library(car)
library(MASS)
library(sp)
library(spdep)
library(ape)

plot(density(data.clean$Pop.size))
hist(data.clean$Pop.size)
dotchart(data.clean$Pop.size)

# data.clean.no = data.clean %>% filter(Landscape != 'NOM37')

#colony.mod = glm.nb(Pop.size ~ Shannon + uHe + Kooperativ_flower_field, data = data.clean)
# colony.mod = lm(Colony.dens ~ poly(Shannon,2) + uHe, data = data.clean.no)
# colony.mod = lm(Colony.dens ~ Forest, data = data.clean.no)
# colony.mod = lm(Colony.dens ~ Good_habitat + Forest, data = data.clean)
# colony.mod = lm(Colony.dens ~ log1p(SNH_combined) + log1p(Flowers_combined) + log1p(Grassland_combined), data = data.clean.no)
# vif(colony.mod)
# summary(colony.mod)
# 
# #cor(data.clean[,-1])
# 
# 1500 km
# edge density
# no of good patches / area of good patches
# temp?
# 
# testDispersion(colony.mod)
# plot(simulateResiduals(colony.mod))
# 
# plot(allEffects(colony.mod))
# 
# data.pred = data.clean
# data.pred$colony.pred = predict(colony.mod, newdata = expand.grid(Good_habitat = seq(min(data.pred$Good_habitat), max(data.pred$Good_habitat), length.out = 31),
#                                                                  Forest = mean(data.pred$Forest)))
# 
# 
# ggplot(data.pred, aes(Good_habitat, colony.pred))+
#   geom_smooth(method = 'lm')
# 
# ggplot(data.pred, aes(Shannon, Pop.size))+
#   geom_smooth()
# 
# ggplot(data.clean, aes(n, Pop.size))+
#   geom_smooth(method = 'lm')


#### He

hist(data$uHe)
hist(log(data$uHe))

he.mod = lm(uHe ~ forest + Good_habitat, data)
he.mod = lm(uHe ~ Kooperativ_flower_field, data)
vif(he.mod)
summary(he.mod)
testDispersion(he.mod)
plot(simulateResiduals(he.mod))
plot(allEffects(he.mod))




### autocorrelation

coord = coord %>% filter(Landscape_ID %in% data.clean$Landscape)
matrix.dist = as.matrix(dist(cbind(coord$x, coord$y)))
matrix.dist.inv <- 1/matrix.dist
diag(matrix.dist.inv) <- 0

# decide what distance you want to use to define your neighborhood structure, e.g.
myDist = 1000

# run non-spatial model
nonsp_glm <- glm.nb(Pop.size ~ Shannon + uHe, data = data.clean)
# caluculate global Moran's I of nonsp model residuals
nonsp_moranI = Moran.I(resid(nonsp_glm), matrix.dist.inv)
# calculate residuals autocovariate (RAC)
rac <- autocov_dist(resid(nonsp_glm), as.matrix(coord[,c(2,3)]), nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T) 
# run RAC-model
sp_glm <- glm.nb(Pop.size ~ Shannon + uHe + rac, data.clean)
summary(sp_glm)
# caluculate global Moran's I of nonsp model residuals
sp_moranI = Moran.I(resid(sp_glm), matrix.dist.inv)
# calculate % Deviance explained

# plotting ea correlogram
correlogram <- correlog(coord$x,coord$y,resid(nonsp_glm),increment=500, resamp = 100, latlon=TRUE)
plot(correlogram)
