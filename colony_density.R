library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)

clusters = read.csv('Data/my_project.BestConfig', header = T, sep = '')  %>% mutate(OffspringID = sub('_.*','', OffspringID))
pop = read.table('Data/populations.txt', header = F, sep = '\t') %>% mutate(V1 = sub('_.*','', V1)) %>%
  rename(M = V1, pop = V2) 

clusters_all = clusters %>% group_by(ClusterIndex) %>% summarise(cluster = paste0(OffspringID, collapse = ',')) %>%
  separate_wider_delim(cluster, delim = ",", too_few = "align_start", names = c('M1', 'M2', "M3", "M4", "M5", "M6"))

clusters2 = clusters_all %>%
  drop_na(M2)


clusters3 = clusters2 %>%
  left_join(pop, by = join_by('M1' == 'M')) %>% rename(pop1 = pop) %>%
  left_join(pop, by = join_by('M2' == 'M')) %>% rename(pop2 = pop) %>%
  left_join(pop, by = join_by('M3' == 'M')) %>% rename(pop3 = pop) %>%
  left_join(pop, by = join_by('M4' == 'M')) %>% rename(pop4 = pop) %>%
  left_join(pop, by = join_by('M5' == 'M')) %>% rename(pop5 = pop) %>%
  left_join(pop, by = join_by('M6' == 'M')) %>% rename(pop6 = pop) %>%
  mutate(same_site = ifelse(pop1 == pop2, T, F)) %>% filter(same_site == T)

n_siblings = clusters3 %>% mutate(across(M1:M6, function(x) ifelse(is.na(x), 0, 1))) %>% select(!c(pop2:pop6))
n_siblings$n_siblings = rowSums(n_siblings[,c(2:7)])

n_sibships = n_siblings %>% group_by(pop1) %>% summarise(n_sibships = n())


sibship_comparison = pop %>%
  filter(M %in% clusters$OffspringID) %>% distinct(pop) %>% left_join(n_sibships, by = join_by('pop' == 'pop1')) %>%
  left_join(n_siblings %>% group_by(pop1) %>% summarise(n_siblings = sum(n_siblings)), by = join_by('pop' == 'pop1')) %>%
  mutate(across(n_sibships:n_siblings, function(x) ifelse(is.na(x), 0, x))) %>% 
  left_join(pop %>% filter(M %in% clusters$OffspringID) %>% group_by(pop) %>% summarise(n_caught = n()), by = join_by('pop' == 'pop'))


sibship_comparison %>%
  ggplot(aes(n_caught, n_sibships))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = 'No. of bees analysed per site', y = 'No. of clusters per site', title = '"Colony" results')

write.csv(sibship_comparison, file = 'Data/sibships_colony.csv', row.names = F)

ggsave(file = 'Data/colony_res.png', height = 4, width = 4)



########### INFERRING COLONY NUMBER PER POP ###############


capwire.dat = clusters_all %>%
  left_join(pop, by = join_by('M1' == 'M')) %>% rename(pop1 = pop) %>%
  left_join(pop, by = join_by('M2' == 'M')) %>% rename(pop2 = pop) %>%
  mutate(same_site = ifelse(pop1 == pop2, T, F)) %>% filter(same_site == T | is.na(same_site)) %>% 
  mutate(across(M1:M6, function(x) ifelse(is.na(x), 0, 1))) %>% select(!pop2) 

capwire.dat$n_siblings = rowSums(capwire.dat[,c(2:7)])

capwire.dat2 = capwire.dat %>% mutate(n_siblings = as.factor(n_siblings)) %>%
  group_by(pop1, n_siblings) %>% summarise(n = n())

capwire.list = split(capwire.dat2, capwire.dat2$pop1)

capwire.list2 = lapply(capwire.list, function(x) x %>% ungroup() %>% select(!pop1) %>% mutate(n_siblings = as.integer(as.character(n_siblings))))

library(capwire)

capwire.results = list()
for (i in 1:length(capwire.list2)) {
  capwire.results[[i]] = fitTirm(as.data.frame(capwire.list2[[i]]), max.pop = 200, max.iter = 20)
}

names(capwire.results) = names(capwire.list2)

capwire.results2 = lapply(capwire.results, function(x) as.data.frame(rbind(x)))
capwire.results3 = capwire.results2 %>% bind_rows(.id = 'pop')
rownames(capwire.results3) = NULL

capwire.results4 = capwire.results3 %>% mutate(ml.pop.size = ifelse(is.na(likelihood), NA, ml.pop.size))
