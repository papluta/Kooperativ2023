library(dplyr)
library(tidyverse)
library(geosphere)
library(ggeffects)
library(compositions)
library(effects)
library(car)
library(DHARMa)
library(performance)
library(betareg)
library(cowplot)


he = read.csv('Data/heterozygosity.csv')
lst = read.csv('Data/lst.csv')
land = read.csv("Data/landuse_1000m.csv")
coord = read.csv("Data/coordinates.csv")
landuse_edge = read.csv("Data/landuse_edge.csv")

data <- he %>% left_join(coord, by = "pop") %>%
  left_join(land, by = join_by("pop" == "Landscape")) %>%
  left_join(lst, by = "pop") %>%
  left_join(landuse_edge, by = join_by("pop" == "Landscape")) %>%
  mutate(forest_area = forest*(pi*1000^2),
         crop_conven_area = crop_conven*(pi*1000^2))

data_alr <- data %>%
  mutate(urban = urban + road,
         snh_per = snh_per + flower_fields) %>%
  mutate(across(crop_conven:water, ~ ifelse(.x == 0, 0.001, .x))) %>%
  mutate(across(flower_fields:water, ~ log(.x/crop_conven), .names = "{.col}_alr")) %>%
  mutate(forest_par = scale(forest_edge/(forest*(pi*1000^2)))[,1],
         forest_ed = scale(forest_edge/(pi*1000^2))[,1],
         crop_par = scale(crop_conven_edge/(forest*(pi*1000^2)))[,1],
         crop_ed = scale(crop_conven_edge/(pi*1000^2))[,1],
         crop_forest_ed = scale((crop_conven_edge+forest_edge)/(pi*1000^2))[,1]) %>%
  mutate(across(ends_with("alr"), ~scale(.x)[,1]),
         LST_historical_daytime_scaled = scale(LST_historical_daytime)[,1])


## LST ~ land use
te.mod = lm(LST_historical_daytime ~ crop_organic_alr + forest_alr + urban_alr +  
            + water_alr + snh_per_alr + x + y,  data = data_alr)
summary(te.mod)
vif(te.mod)

performance::check_model(te.mod, residual_type = "normal")
testDispersion(te.mod)
plot(simulateResiduals(te.mod))
plot(te.mod)

sum <- summary(te.mod)

summary_tbl <- as.data.frame(sum$coefficients) %>% rownames_to_column(var = "Predictor") %>%
  mutate(across(where(is.numeric), ~ round(.x,3)))

write.csv(summary_tbl, file = "Results/lst_model_summary.csv", row.names = F)


## adhHe ~ land use + LST

he.mod = betareg(adjHe  ~ LST_historical_daytime_scaled + forest_alr +  crop_organic_alr + urban_alr +
                 water_alr + snh_per_alr + crop_forest_ed,  data = data_alr, link = 'loglog')

he.mod.no.lst = betareg(adjHe  ~ forest_alr +  crop_organic_alr + urban_alr +
                   water_alr + snh_per_alr + crop_forest_ed,  data = data_alr, link = 'loglog')


vif(he.mod)
summary(he.mod)
summary(he.mod.no.lst)

performance::check_model(he.mod, residual_type = "normal")
performance::check_model(he.mod.no.lst, residual_type = "normal")

sum <- summary(he.mod)
summary_tbl <- as.data.frame(sum$coefficients$mean) %>% rownames_to_column(var = "Predictor") %>%
  mutate(across(where(is.numeric), ~ round(.x,3)))
write.csv(summary_tbl, file = "Results/he_model_summary.csv", row.names = F)


sum <- summary(he.mod.no.lst)
summary_tbl <- as.data.frame(sum$coefficients$mean) %>% rownames_to_column(var = "Predictor") %>%
  mutate(across(where(is.numeric), ~ round(.x,3)))
write.csv(summary_tbl, file = "Results/he_model_no_lst_summary.csv", row.names = F)

## plotting
text = summary(te.mod)
text$coefficients[,4]

pred.table.forest = ggpredict(te.mod, terms = "forest_alr [all]")

forest.plot = ggplot(pred.table.forest, aes(x, predicted))+
  geom_smooth(color = 'black', fill='black')+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = '#6EA130', alpha = 0.4)+
  geom_point(data = data_alr, aes(forest_alr, LST_historical_daytime), shape = 21, col = 'black')+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0, face = "bold"),plot.title.position = "plot")+
  labs(x = "Forest area (relative to\nconventional agriculture, ALR)", y = "Historical LST", title = "A")+
  scale_x_continuous(n.breaks = 6)+
  annotate("text", x = -0.5, y = 21, label = paste0('t = ',round(text$coefficients["forest_alr",3], 3),
                                                    ',  p = ',round(text$coefficients["forest_alr",4], 3)),
           size = 3)


text = summary(he.mod)
text$coefficients$mean[,4]

pred.table.lst = as.data.frame(ggpredict(he.mod, terms = "LST_historical_daytime_scaled [all]")) %>%
  mutate(x_unscaled = x * sd(data_alr$LST_historical_daytime) + mean(data_alr$LST_historical_daytime))


lst.plot = ggplot(pred.table.lst, aes(x_unscaled, predicted))+
  geom_smooth(color = 'black', fill='black')+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = '#306EA1', alpha = 0.4)+
  geom_point(data = data_alr, aes(LST_historical_daytime, adjHe), shape = 21, col = 'black')+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0, face = "bold"),plot.title.position = "plot")+
  labs(x = "Historical LST", y = expression(italic(adj) * H[EXP]), title = "B")+
  scale_x_continuous(n.breaks =  5)+
  annotate("text", x = 18.8, y = 0.252, label = paste0('z = ',round(text$coefficients$mean["LST_historical_daytime_scaled",3], 3),
                                                     ',  p = ',round(text$coefficients$mean["LST_historical_daytime_scaled",4], 3)),
           size = 3)

pred.table.forest2 = ggpredict(he.mod.no.lst, terms = "forest_alr [all]")

text = summary(he.mod.no.lst)
text$coefficients$mean[,4]

forest2.plot = ggplot(pred.table.forest2, aes(x, predicted))+
  geom_smooth(color = 'black', fill='black')+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = '#BF9200', alpha = 0.4)+
  geom_point(data = data_alr, aes(forest_alr, adjHe), shape = 21, col = 'black')+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0, face = "bold"),plot.title.position = "plot")+
  labs(x = "Forest area (relative to\nconventional agriculture, ALR)", y = expression(italic(adj) * H[EXP]), title = "D")+
  scale_x_continuous(n.breaks = 6)+
  annotate("text", x = -0.4, y = 0.252, label = paste0('z = ',round(text$coefficients$mean["forest_alr",3], 3),
                                                     ',  p = ',round(text$coefficients$mean["forest_alr",4], 3)),
           size = 3)


pred.table.water = ggpredict(he.mod, terms = "water_alr [all]") 

water.plot = ggplot(pred.table.water, aes(x, predicted))+
  geom_smooth(color = 'black', fill='black')+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = '#306EA1', alpha = 0.4)+
  geom_point(data = data_alr, aes(water_alr, adjHe), shape = 21, col = 'black')+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0, face = "bold"),plot.title.position = "plot")+
  labs(x = "Water body area (relative to\nconventional agriculture, ALR)", y = expression(italic(adj) * H[EXP]), title = "C")+
  scale_x_continuous(n.breaks = 5)+
  annotate("text", x = -0.1, y = 0.252, label = paste0('z = ',round(text$coefficients$mean["water_alr",3], 3),
                                                       ',  p = ',round(text$coefficients$mean["water_alr",4], 3)),
           size = 3)

plot_grid(forest.plot, lst.plot,  water.plot,forest2.plot,   nrow = 2, align = 'hv')

ggsave(file = "Results/reg_plots.png", width = 8, height = 8)

