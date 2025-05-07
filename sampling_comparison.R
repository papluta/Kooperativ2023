library(dplyr)
library(ggplot2)

He_full = read.csv('Data/Goe_filtered095_Het60_pruned_heterozygosity.csv')
He_size3 = read.csv('Data/size3_results.csv')
He_size9 = read.csv('Data/size9_results.csv')
He_size15 = read.csv('Data/size15_results.csv')
He_size21 = read.csv('Data/size21_results.csv')
He3.m = He_size3 %>% group_by(pop) %>% summarise(uHeSD = mean(uHeSD), uHe = mean(uHe), HeSD = mean(HeSD), He = mean(He)) %>% mutate(size = '3')
He9.m = He_size9 %>% group_by(pop) %>% summarise(uHeSD = mean(uHeSD), uHe = mean(uHe), HeSD = mean(HeSD), He = mean(He)) %>% mutate(size = '9')
He15.m = He_size15 %>% group_by(pop) %>% summarise(uHeSD = mean(uHeSD), uHe = mean(uHe), HeSD = mean(HeSD), He = mean(He)) %>% mutate(size = '15')
He21.m = He_size21 %>% group_by(pop) %>% summarise(uHeSD = mean(uHeSD), uHe = mean(uHe), HeSD = mean(HeSD), He = mean(He)) %>% mutate(size = '21')
He.m = He_full %>% filter(pop %in% He3.m$pop) %>% select(pop, uHe, uHeSD, He, HeSD) %>% mutate(size = 'full')

comp = rbind(He.m, He21.m,  He15.m, He9.m, He3.m)
comp.w = comp %>% pivot_wider(names_from = size, values_from = c(uHe, uHeSD, He, HeSD)) %>% select(-contains("SD"))

comp %>% mutate(size = factor(size, levels = c('3', '9','15', '21'))) %>% mutate(name = paste0(size, pop)) %>%
  ggplot(aes(name, He))+
  geom_pointrange(aes(ymin = uHe - HeSD, ymax = uHe + HeSD))

comp$uHe.s = scale(comp$uHe)

write.csv(comp.w, 'Results/sampling_comaprison_He.csv', row.names = F)

